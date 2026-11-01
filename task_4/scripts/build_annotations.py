from __future__ import annotations

import argparse
import csv
import gzip
from collections import defaultdict
from pathlib import Path
from urllib.parse import unquote


PROJECT_DIR = Path(__file__).resolve().parent.parent
DEFAULT_VCF = PROJECT_DIR / "data" / "student_20.vcf"
DEFAULT_VEP = PROJECT_DIR / "results" / "vep.tsv.gz"
DEFAULT_RESULTS = PROJECT_DIR / "results"

IMPACT_RANK = {"HIGH": 0, "MODERATE": 1, "LOW": 2, "MODIFIER": 3}

CONSEQUENCE_ORDER = [
    "transcript_ablation",
    "splice_acceptor_variant",
    "splice_donor_variant",
    "stop_gained",
    "frameshift_variant",
    "stop_lost",
    "start_lost",
    "transcript_amplification",
    "inframe_insertion",
    "inframe_deletion",
    "missense_variant",
    "protein_altering_variant",
    "splice_donor_5th_base_variant",
    "splice_region_variant",
    "splice_donor_region_variant",
    "splice_polypyrimidine_tract_variant",
    "incomplete_terminal_codon_variant",
    "start_retained_variant",
    "stop_retained_variant",
    "synonymous_variant",
    "coding_sequence_variant",
    "mature_miRNA_variant",
    "5_prime_UTR_variant",
    "3_prime_UTR_variant",
    "non_coding_transcript_exon_variant",
    "intron_variant",
    "NMD_transcript_variant",
    "non_coding_transcript_variant",
    "upstream_gene_variant",
    "downstream_gene_variant",
    "TFBS_ablation",
    "TFBS_amplification",
    "TF_binding_site_variant",
    "regulatory_region_ablation",
    "regulatory_region_amplification",
    "feature_elongation",
    "regulatory_region_variant",
    "feature_truncation",
    "intergenic_variant",
]
CONSEQUENCE_RANK = {name: rank for rank, name in enumerate(CONSEQUENCE_ORDER)}

TRANSCRIPT_RULE = (
    "MANE Select > MANE Plus Clinical > canonical > protein-coding > other; "
    "within the selected class: VEP consequence severity"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--vcf", type=Path, default=DEFAULT_VCF)
    parser.add_argument("--vep", type=Path, default=DEFAULT_VEP)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_RESULTS)
    return parser.parse_args()


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return path.open("r", encoding="utf-8", newline="")


def read_vcf(path: Path) -> list[dict[str, str]]:
    variants = []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            chrom, pos, variant_id, ref, alt, *_ = line.rstrip("\n").split("\t")
            variants.append(
                {
                    "chrom": chrom.removeprefix("chr"),
                    "pos": pos,
                    "id": variant_id,
                    "ref": ref,
                    "alt": alt,
                }
            )
    return variants


def read_vep(path: Path) -> list[dict[str, str]]:
    with open_text(path) as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def consequence_rank(value: str) -> int:
    parts = value.split("&") if value else []
    return min(
        (CONSEQUENCE_RANK.get(part, len(CONSEQUENCE_RANK)) for part in parts),
        default=len(CONSEQUENCE_RANK),
    )


def transcript_priority(row: dict[str, str]) -> tuple[int, str]:
    if row.get("MANE_SELECT"):
        return 0, "MANE_Select"
    if row.get("MANE_PLUS_CLINICAL"):
        return 1, "MANE_Plus_Clinical"
    if row.get("CANONICAL") in {"YES", "1"}:
        return 2, "canonical"
    if row.get("BIOTYPE") == "protein_coding":
        return 3, "protein_coding_fallback"
    return 4, "other_fallback"


def row_sort_key(row: dict[str, str]) -> tuple[int, int, int, str, str]:
    priority, _ = transcript_priority(row)
    return (
        priority,
        IMPACT_RANK.get(row.get("IMPACT", ""), len(IMPACT_RANK)),
        consequence_rank(row.get("Consequence", "")),
        row.get("SYMBOL", ""),
        row.get("Feature", ""),
    )


def clean_list(value: str) -> str:
    return unquote(value).replace("+", "; ") if value else ""


def first_rsid(*values: str) -> str:
    for value in values:
        for token in value.replace(",", "&").split("&"):
            token = token.strip()
            if token.startswith("rs") and token[2:].isdigit():
                return token
    return ""


def write_tsv(path: Path, rows: list[dict[str, str]], fields: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    args = parse_args()
    variants = read_vcf(args.vcf)
    vep_rows = read_vep(args.vep)

    if len(variants) != 50:
        raise ValueError(f"Expected 50 VCF variants, found {len(variants)}")

    by_id: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in vep_rows:
        by_id[row["Uploaded_variation"]].append(row)

    final_rows = []
    clinvar_rows = []
    gnomad_rows = []

    for variant in variants:
        variant_id = variant["id"]
        candidates = by_id.get(variant_id, [])
        if not candidates:
            coordinate = f'{variant["chrom"]}:{variant["pos"]}'
            candidates = [
                row
                for row in vep_rows
                if row.get("Location") == coordinate
                and row.get("Ref") == variant["ref"]
                and row.get("Allele") == variant["alt"]
            ]
        if not candidates:
            raise ValueError(f"No VEP annotation found for {variant_id}")

        selected = min(candidates, key=row_sort_key)
        _, selected_by = transcript_priority(selected)
        rsid = first_rsid(selected.get("Existing_variation", ""), variant_id)
        label = rsid or (
            f'{variant["chrom"]}:{variant["pos"]} '
            f'{variant["ref"]}>{variant["alt"]}'
        )

        exomes_af = selected.get("gnomAD_exomes_AF", "")
        genomes_af = selected.get("gnomAD_genomes_AF", "")
        if genomes_af:
            selected_af, af_source = genomes_af, "gnomAD_genomes_AF"
        elif exomes_af:
            selected_af, af_source = exomes_af, "gnomAD_exomes_AF"
        else:
            selected_af, af_source = "", "not_available"

        clinvar_id = selected.get("ClinVar", "")
        common = {
            "variant": label,
            "chrom": variant["chrom"],
            "pos": variant["pos"],
            "ref": variant["ref"],
            "alt": variant["alt"],
            "rsid": rsid,
        }

        final_rows.append(
            {
                **common,
                "gene": selected.get("SYMBOL", ""),
                "selected_transcript": selected.get("Feature", ""),
                "selected_by": selected_by,
                "transcript_rule": TRANSCRIPT_RULE,
                "vep_consequence": selected.get("Consequence", ""),
                "vep_impact": selected.get("IMPACT", ""),
                "hgvsc": unquote(selected.get("HGVSc", "")),
                "hgvsp": unquote(selected.get("HGVSp", "")),
                "gnomad_allele_frequency": selected_af,
                "gnomad_af_source": af_source,
                "clinvar_present": "yes" if clinvar_id else "no",
                "clinvar_variation_id": clinvar_id,
                "clinvar_significance": clean_list(
                    selected.get("ClinVar_CLNSIG", "")
                ),
                "clinvar_condition": clean_list(selected.get("ClinVar_CLNDN", "")),
                "clinvar_review_status": clean_list(
                    selected.get("ClinVar_CLNREVSTAT", "")
                ),
            }
        )

        clinvar_rows.append(
            {
                **common,
                "clinvar_present": "yes" if clinvar_id else "no",
                "clinvar_variation_id": clinvar_id,
                "clinical_significance": clean_list(
                    selected.get("ClinVar_CLNSIG", "")
                ),
                "conflicting_interpretations": clean_list(
                    selected.get("ClinVar_CLNSIGCONF", "")
                ),
                "review_status": clean_list(
                    selected.get("ClinVar_CLNREVSTAT", "")
                ),
                "conditions": clean_list(selected.get("ClinVar_CLNDN", "")),
                "rcv_records": clean_list(selected.get("ClinVar_CLNRCV", "")),
            }
        )

        gnomad_rows.append(
            {
                **common,
                "gnomad_exomes_variant": selected.get("gnomAD_exomes", ""),
                "gnomad_exomes_af": exomes_af,
                "gnomad_genomes_variant": selected.get("gnomAD_genomes", ""),
                "gnomad_genomes_af": genomes_af,
                "selected_af": selected_af,
                "selected_af_source": af_source,
            }
        )

    args.output_dir.mkdir(parents=True, exist_ok=True)
    final_fields = [
        "variant", "chrom", "pos", "ref", "alt", "gene",
        "selected_transcript", "selected_by", "transcript_rule",
        "vep_consequence", "vep_impact", "hgvsc", "hgvsp", "rsid",
        "gnomad_allele_frequency", "gnomad_af_source", "clinvar_present",
        "clinvar_variation_id", "clinvar_significance", "clinvar_condition",
        "clinvar_review_status",
    ]
    clinvar_fields = [
        "variant", "chrom", "pos", "ref", "alt", "rsid",
        "clinvar_present", "clinvar_variation_id", "clinical_significance",
        "conflicting_interpretations", "review_status", "conditions",
        "rcv_records",
    ]
    gnomad_fields = [
        "variant", "chrom", "pos", "ref", "alt", "rsid",
        "gnomad_exomes_variant", "gnomad_exomes_af",
        "gnomad_genomes_variant", "gnomad_genomes_af", "selected_af",
        "selected_af_source",
    ]

    write_tsv(args.output_dir / "final_annotations.tsv", final_rows, final_fields)
    write_tsv(
        args.output_dir / "clinvar_annotations.tsv", clinvar_rows, clinvar_fields
    )
    write_tsv(
        args.output_dir / "gnomad_annotations.tsv", gnomad_rows, gnomad_fields
    )

    print(f"VEP annotation rows: {len(vep_rows)}")
    print(f"Final variants: {len(final_rows)}")
    print(
        "ClinVar records: "
        f'{sum(row["clinvar_present"] == "yes" for row in clinvar_rows)}'
    )
    print(
        "gnomAD selected AF available: "
        f'{sum(bool(row["selected_af"]) for row in gnomad_rows)}'
    )


if __name__ == "__main__":
    main()
