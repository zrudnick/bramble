#!/usr/bin/env python3
import subprocess
from pathlib import Path

ref = Path('/Users/rob/software_src/bramble/test_data/nanosim.chr21.ref_projected_long.bam')
rust = Path('/Users/rob/software_src/bramble/test_data/nanosim.chr21.rust_projected_long.bam')

exclude = {
    "NM_001321704_6_aligned_6105603_F_2_324_1",
    "NM_003720_892_aligned_9505676_R_1_561_19",
    "NM_004965_389_aligned_5668657_R_2_806_3",
    "NM_004965_773_aligned_2459988_R_47_457_87",
    "NM_006052_1560_aligned_6543549_R_4_1363_1",
    "NM_006936_345_aligned_1900844_R_3_644_6",
    "NM_013240_831_aligned_7250963_R_6_2401_2",
    "NM_017446_761_aligned_10160222_R_5_219_2",
    "NM_021254_1171_aligned_1575360_R_6_1208_7",
    "NR_024027_104_aligned_1129087_R_5_131_2",
    "NR_024027_1184_aligned_1645391_R_2_241_9",
}

remove_tags = {"HI"}
mask_secondary = 0x100


def normalize(line: str) -> str:
    fields = line.rstrip("\n").split("\t")
    qname = fields[0]
    if qname in exclude:
        return ""
    flag = int(fields[1])
    flag &= ~mask_secondary
    fields[1] = str(flag)
    fixed = fields[:11]
    tags = [t for t in fields[11:] if t[:2] not in remove_tags]
    tags.sort()
    return "\t".join(fixed + tags)


def read_sam(path: Path):
    proc = subprocess.Popen(["samtools", "view", str(path)], stdout=subprocess.PIPE, text=True)
    for line in proc.stdout:
        norm = normalize(line)
        if norm:
            yield norm
    proc.stdout.close()
    proc.wait()


def main() -> int:
    ref_lines = sorted(read_sam(ref))
    rust_lines = sorted(read_sam(rust))

    s_ref = set(ref_lines)
    s_rust = set(rust_lines)
    only_ref = s_ref - s_rust
    only_rust = s_rust - s_ref

    print("ref_only", len(only_ref))
    print("rust_only", len(only_rust))
    if only_ref or only_rust:
        print("EXAMPLE REF")
        for line in sorted(only_ref)[:3]:
            print(line)
        print("EXAMPLE RUST")
        for line in sorted(only_rust)[:3]:
            print(line)
        return 1
    print("MATCH_EXCLUDING_HI_AND_SECONDARY_FLAG")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
