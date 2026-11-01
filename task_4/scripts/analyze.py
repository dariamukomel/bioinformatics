from pathlib import Path
from collections import Counter

path = Path(__file__).resolve().parent.parent / "data" / "student_20.vcf"

counts = Counter()
transitions = {"AG", "GA", "CT", "TC"}

with path.open() as f:
    for line in f:
        if line.startswith("#"):
            continue

        chrom, pos, rsid, ref, alt, *_ = line.strip().split("\t")

        counts["total"] += 1
        counts[f"chr{chrom.removeprefix('chr')}"] += 1

        pair = ref + alt

        if pair in transitions:
            counts["transitions"] += 1
        elif len(ref) == len(alt) == 1:
            counts["transversions"] += 1

for key in ["total", "chr21", "chr22",
            "transitions", "transversions"]:
    print(f"{key}: {counts[key]}")
