"""
What end spread did TAMA actually absorb?  -> pick -a and -z from it.
 
tama_trans_report.txt reports per-exon wobble. The FIRST entry of
start_wobble_list and the LAST entry of end_wobble_list are the transcript's
two terminal wobbles; strand decides which is 5' and which is 3'.
 
Run this on a GENEROUS merge (large -a/-z): the distribution is then only
right-censored at those thresholds, so where it tapers tells you how small
you can go without losing anything.
 
  python3 wobble.py tama_trans_report.txt merged.bed
"""
import sys
from collections import Counter

report, bed = sys.argv[1], sys.argv[2]
report, bed = '/projects/splitorfs/work/PacBio/merged_bam_files/merge_mando_stringtie_isoquant_rescue_up_1000_down_500_25_august_2026/HUVEC/HUVEC_merged_tama_trans_report.txt', '/projects/splitorfs/work/PacBio/merged_bam_files/merge_mando_stringtie_isoquant_rescue_up_1000_down_500_25_august_2026/HUVEC/HUVEC_merged_tama.bed'
report, bed = '/projects/splitorfs/work/PacBio/merged_bam_files/merge_mando_stringtie_isoquant_rescue_up_2000_down_500_longest_ends_05_09_2026/HUVEC/HUVEC_merged_tama_trans_report.txt', '/projects/splitorfs/work/PacBio/merged_bam_files/merge_mando_stringtie_isoquant_rescue_up_2000_down_500_longest_ends_05_09_2026/HUVEC/HUVEC_merged_tama.bed'
report, bed = '/projects/splitorfs/work/PacBio/merged_bam_files/merge_mando_stringtie_isoquant_rescue_up_2000_down_1000_longest_ends_05_09_2026/HUVEC/HUVEC_merged_tama_trans_report.txt', '/projects/splitorfs/work/PacBio/merged_bam_files/merge_mando_stringtie_isoquant_rescue_up_2000_down_1000_longest_ends_05_09_2026/HUVEC/HUVEC_merged_tama.bed'
report, bed = '/projects/splitorfs/work/PacBio/merged_bam_files/merge_mando_stringtie_isoquant_rescue_up_3000_down_2000_longest_ends_05_09_2026/HUVEC/HUVEC_merged_tama_trans_report.txt', '/projects/splitorfs/work/PacBio/merged_bam_files/merge_mando_stringtie_isoquant_rescue_up_3000_down_2000_longest_ends_05_09_2026/HUVEC/HUVEC_merged_tama.bed'


strand = {}
for line in open(bed):
    f = line.split()
    if len(f) >= 6:
        if len(f[3].split(";")) > 1:
            strand[f[3].split(";")[1]] = f[5]
        else:
            print(f[3])

five, three, skipped = [], [], 0
for i, line in enumerate(open(report)):
    if i == 0 and line.startswith("transcript_id"):
        continue
    f = line.rstrip("\n").split("\t")
    tid, n = f[0], int(f[1])
    if n < 2:                       # single contributor: nothing was merged
        continue
    s = strand.get(tid.split(";")[-1])
    if s is None:
        skipped += 1
        continue
    sw = [int(x) for x in f[3].split(",") if x.strip() != ""]
    ew = [int(x) for x in f[4].split(",") if x.strip() != ""]
    if not sw or not ew:
        continue
    left, right = sw[0], ew[-1]     # genomic-left start, genomic-right end
    if s == "+":
        five.append(left)
        three.append(right)
    else:
        five.append(right)
        three.append(left)


def summarise(name, v, thresholds):
    v = sorted(v)
    n = len(v)
    if n == 0:
        print(f"{name}: no multi-contributor models")
        return

    def q(p): return v[int(p * (n - 1))]
    print(f"\n{name}  (n = {n})")
    print(
        f"  median {q(.50)}   75th {q(.75)}   90th {q(.90)}   99th {q(.99)}   max {v[-1]}")
    for t in thresholds:
        k = sum(1 for x in v if x > t)
        print(f"  wobble > {t:>5} nt : {k:>7} ({100*k/n:5.1f}%)")


summarise("5' end wobble", five,  [100, 250, 500, 1000, 2000, 3000])
summarise("3' end wobble", three, [50, 100, 250, 500, 1000, 2000])
if skipped:
    print(f"\n({skipped} report rows had no matching BED id — check id formats)")
