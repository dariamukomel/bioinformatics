from pathlib import Path
import subprocess


SAMPLE = "SRR001666"

PROJECT_ROOT = Path(__file__).resolve().parents[5]

REF = PROJECT_ROOT / "data/reference/ecoli_reference.fna"
READ1 = PROJECT_ROOT / "data/fastq/SRR001666_1.fastq.gz"
READ2 = PROJECT_ROOT / "data/fastq/SRR001666_2.fastq.gz"

FASTQC_DIR = PROJECT_ROOT / "results/fastqc"
ALIGN_DIR = PROJECT_ROOT / "results/alignment"
FLAGSTAT_DIR = PROJECT_ROOT / "results/flagstat"
VARIANT_DIR = PROJECT_ROOT / "results/variants"
LOG_DIR = PROJECT_ROOT / "logs"

SAM_FILE = ALIGN_DIR / f"{SAMPLE}.kedro.sam"
BAM_FILE = ALIGN_DIR / f"{SAMPLE}.kedro.bam"
SORTED_BAM_FILE = ALIGN_DIR / f"{SAMPLE}.kedro.sorted.bam"
FLAGSTAT_FILE = FLAGSTAT_DIR / f"{SAMPLE}_kedro_flagstat.txt"
MAPPING_RESULT_FILE = FLAGSTAT_DIR / f"{SAMPLE}_kedro_mapping_quality.txt"
VCF_FILE = VARIANT_DIR / f"{SAMPLE}.kedro.vcf"


def run_command(command: list[str], log_file: Path | None = None) -> None:
    if log_file is not None:
        with log_file.open("w") as log:
            subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True)
    else:
        subprocess.run(command, check=True)


def prepare_directories() -> str:
    for directory in [FASTQC_DIR, ALIGN_DIR, FLAGSTAT_DIR, VARIANT_DIR, LOG_DIR]:
        directory.mkdir(parents=True, exist_ok=True)

    return "directories_prepared"


def run_fastqc(_: str) -> str:
    run_command(
        [
            "fastqc",
            str(READ1),
            str(READ2),
            "-o",
            str(FASTQC_DIR),
        ],
        LOG_DIR / "kedro_fastqc.log",
    )

    return "fastqc_finished"


def index_reference(_: str) -> str:
    if not Path(f"{REF}.bwt").exists():
        run_command(
            ["bwa", "index", str(REF)],
            LOG_DIR / "kedro_bwa_index.log",
        )

    return "reference_indexed"


def run_bwa_mem(_: str) -> str:
    with SAM_FILE.open("w") as sam_output, (LOG_DIR / "kedro_bwa_mem.log").open("w") as log:
        subprocess.run(
            [
                "bwa",
                "mem",
                str(REF),
                str(READ1),
                str(READ2),
            ],
            stdout=sam_output,
            stderr=log,
            check=True,
        )

    return str(SAM_FILE)


def convert_sam_to_bam(sam_file: str) -> str:
    with BAM_FILE.open("w") as bam_output:
        subprocess.run(
            [
                "samtools",
                "view",
                "-S",
                "-b",
                sam_file,
            ],
            stdout=bam_output,
            check=True,
        )

    return str(BAM_FILE)


def run_flagstat(bam_file: str) -> str:
    with FLAGSTAT_FILE.open("w") as flagstat_output:
        subprocess.run(
            [
                "samtools",
                "flagstat",
                bam_file,
            ],
            stdout=flagstat_output,
            check=True,
        )

    return str(FLAGSTAT_FILE)


def parse_mapping_quality(flagstat_file: str) -> str:
    mapped_percent = None

    with open(flagstat_file, "r") as file:
        for line in file:
            if " mapped (" in line:
                mapped_percent = line.split("(")[1].split("%")[0]
                break

    if mapped_percent is None:
        raise ValueError("Mapped reads percentage was not found in flagstat output")

    status = "OK" if float(mapped_percent) > 90 else "not OK"

    with MAPPING_RESULT_FILE.open("w") as output:
        output.write(f"Mapped reads: {mapped_percent}%\n")
        output.write(f"{status}\n")

    print(f"Mapped reads: {mapped_percent}%")
    print(status)

    return status


def sort_bam(status: str) -> str:
    if status != "OK":
        print("Mapping quality is not OK. Sorting was skipped.")
        return "sorting_skipped"

    subprocess.run(
        [
            "samtools",
            "sort",
            str(BAM_FILE),
            "-o",
            str(SORTED_BAM_FILE),
        ],
        check=True,
    )

    return str(SORTED_BAM_FILE)


def index_sorted_bam(sorted_bam_file: str) -> str:
    if sorted_bam_file == "sorting_skipped":
        print("BAM indexing was skipped.")
        return "bam_indexing_skipped"

    subprocess.run(
        [
            "samtools",
            "index",
            sorted_bam_file,
        ],
        check=True,
    )

    return f"{sorted_bam_file}.bai"


def run_freebayes(bam_index_file: str) -> str:
    if bam_index_file == "bam_indexing_skipped":
        print("FreeBayes was skipped.")
        return "freebayes_skipped"

    if not Path(f"{REF}.fai").exists():
        subprocess.run(
            [
                "samtools",
                "faidx",
                str(REF),
            ],
            check=True,
        )

    with VCF_FILE.open("w") as vcf_output:
        subprocess.run(
            [
                "freebayes",
                "-f",
                str(REF),
                str(SORTED_BAM_FILE),
            ],
            stdout=vcf_output,
            check=True,
        )

    print(f"VCF file: {VCF_FILE}")
    return str(VCF_FILE)