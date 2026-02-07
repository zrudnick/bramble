# Long-read mismatches (11 reads / 60 diffs)

We have 11 reads (60 alignment diffs) where the Rust output disagrees with the C++ output in long-read mode.

Observed behavior:
- These reads are single-exon long reads with no XS/ts tags in the input.
- The transcript strand is `-` but the input read is forward (`orig_rev=false`).
- C++ reverses them (flags 16/272, reversed CIGAR), while Rust matches C++ only if we disable the
  single-exon strand override for these cases.
- However, disabling the override broadly causes large regressions (tens of thousands of diffs), so we
  are keeping the override and accepting these 11 mismatches for now.

Read names:
- NM_001321704_6_aligned_6105603_F_2_324_1
- NM_003720_892_aligned_9505676_R_1_561_19
- NM_004965_389_aligned_5668657_R_2_806_3
- NM_004965_773_aligned_2459988_R_47_457_87
- NM_006052_1560_aligned_6543549_R_4_1363_1
- NM_006936_345_aligned_1900844_R_3_644_6
- NM_013240_831_aligned_7250963_R_6_2401_2
- NM_017446_761_aligned_10160222_R_5_219_2
- NM_021254_1171_aligned_1575360_R_6_1208_7
- NR_024027_104_aligned_1129087_R_5_131_2
- NR_024027_1184_aligned_1645391_R_2_241_9

Next steps (when revisiting):
- Identify a precise condition that matches these reads (e.g., input CIGAR pattern or tags)
  without affecting the broader dataset.
- Confirm with Zoe whether C++ behavior here is intended.
