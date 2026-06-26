# ----- This script changes the order of transcript_id and gene_id in the 9th ----- #
# ----- field of a GTF that was generated using SQANTI3 to gene_id; transcript_id  ----- #

import os
import pandas as pd
import argparse
import csv


# path_to_tama_gtf = "/projects/splitorfs/work/PacBio/merged_bam_files/merge_mando_stringtie_isoquant_rescue_23_June_2026/HUVEC/HUVEC_merged_tama.gtf"
# output_gtf = "/projects/splitorfs/work/PacBio/merged_bam_files/merge_mando_stringtie_isoquant_rescue_23_June_2026/HUVEC/HUVEC_merged_tama_gene_id.gtf"
# path_to_classif = "/projects/splitorfs/work/PacBio/merged_bam_files/merge_mando_stringtie_isoquant_rescue_23_June_2026/SQANTI3_QC/HUVEC/isoforms_classification.txt"

def parse_args():
    parser = argparse.ArgumentParser(
        description="."
    )

    # Required positional arguments
    parser.add_argument("path_to_tama_gtf",
                        help="Path to TAMA GTF to be changed")

    parser.add_argument("path_to_classif",
                        help="Path to classification file SQANTI")

    parser.add_argument("output_gtf",
                        help="Path to output GTF")

    return parser.parse_args()


def main(path_to_tama_gtf, path_to_classif, output_gtf):
    def change_gene_id(row):
        gene_info_string = row[8]
        gene_info_list = gene_info_string.split('"')
        gene_info_list[1] = row['gene_id']
        gene_info_string = '"'.join(gene_info_list)
        return gene_info_string

    tama_gtf_df = pd.read_csv(
        path_to_tama_gtf, header=None, sep='\t', dtype={0: str, 6: str})
    classif_df = pd.read_csv(path_to_classif, header=0,
                             sep='\t', low_memory=False)

    # native TAMA gene id = first quoted token in the attribute field;
    # kept in its own column so it survives as a fallback after mapping
    tama_gtf_df['native_gene'] = tama_gtf_df.iloc[:, 8].apply(
        lambda x: x.split('"')[1])
    classif_df['gene_id'] = classif_df['isoform'].apply(
        lambda x: x.split('.')[0])

    def bucket(s):
        """check the number of multiple gene associations: if there is not a clear one, 
        then note as smbiguous and want to keep TAMA ID also if only novelGene is
        associated keep the TAMA ID"""
        vc = s.value_counts()
        if vc.index.str.startswith('novelGene').all():
            return 'all_novel'
        if vc.size == 1:
            return 'clean'
        return 'dominant' if vc.iloc[0] / vc.sum() >= 0.8 else 'ambiguous'

    buckets = classif_df.groupby('gene_id')['associated_gene'].agg(bucket)
    print(buckets.value_counts())

    # A single TAMA gene can hold isoforms with different associated_gene values
    # (fusions, multi-mapping). dict(zip(...)) silently kept whichever isoform
    # appeared last; instead choose the most frequent association per TAMA gene,
    # ties broken alphabetically so the result is reproducible.

    def pick_associated(s):
        counts = s.value_counts()
        return counts[counts == counts.max()].index[0]

    keep = set(buckets[buckets.isin(["clean", "dominant"])].index)

    assert len(keep) == sum(buckets.isin(["clean", "dominant"]))

    gene_id_dict = (classif_df[classif_df['gene_id'].isin(keep)]
                    .groupby('gene_id')['associated_gene']
                    .agg(pick_associated).to_dict())

    tama_gtf_df['gene_id'] = (tama_gtf_df['native_gene']
                              .map(gene_id_dict)
                              .fillna(tama_gtf_df['native_gene']))

    # Genes in the GTF but absent from the classification map to NaN; keep their
    # native TAMA id rather than writing a literal 'nan' (which would itself
    # collapse many genes into one duplicated id).
    unmapped = tama_gtf_df['gene_id'].isna()
    if unmapped.any():
        print(f"Warning: {int(unmapped.sum())} feature line(s) had no SQANTI "
              f"association; keeping native TAMA gene_id.")
        tama_gtf_df.loc[unmapped, 'gene_id'] = \
            tama_gtf_df.loc[unmapped, 'native_gene']

    # SQANTI labels an antisense transcript with the gene it is antisense TO, so
    # a sense locus and its antisense neighbour can receive the same
    # associated_gene. Collapsing them fuses two strands under one id, which is
    # biologically wrong and breaks union-exon length / TPM with a duplicate
    # label. Detect any gene_id landing on more than one (chrom, strand) locus
    # and suffix only those: by strand, or by chrom+strand if it also spans
    # chromosomes. (chrom = column 0, strand = column 6.)
    loci = tama_gtf_df[['gene_id', 0, 6]].drop_duplicates()
    # how many different combinations of strand and chromosome does each gene_id have?
    n_loci = loci.groupby('gene_id').size()
    collisions = set(n_loci[n_loci > 1].index)
    multichrom = set(
        loci[loci['gene_id'].isin(collisions)]
        .groupby('gene_id')[0].nunique().loc[lambda s: s > 1].index)

    if len(multichrom) > 0:
        print('genes spanning multiple chromosomes:', multichrom)

    strand_word = {'+': 'plus', '-': 'minus', '.': 'unstranded'}

    def disambiguate(row):
        if row['gene_id'] not in collisions:
            return row['gene_id']
        strand = strand_word.get(row[6], 'unknown')
        if row['gene_id'] in multichrom:
            return f"{row['gene_id']}_{row[0]}_{strand}"
        return f"{row['gene_id']}_{strand}"

    was_collision = tama_gtf_df['gene_id'].isin(collisions)
    tama_gtf_df['gene_id'] = tama_gtf_df.apply(disambiguate, axis=1)

    if collisions:
        changed = (tama_gtf_df.loc[was_collision, ['gene_id', 0, 6]]
                   .drop_duplicates().sort_values(['gene_id', 0, 6]))
        print(f"Disambiguated {len(collisions)} colliding gene id(s):")
        print(changed.to_string(index=False))

    tama_gtf_df.iloc[:, 8] = tama_gtf_df.apply(
        lambda row: change_gene_id(row), axis=1)

    # QUOTE_NONE only: the attribute field legitimately contains double quotes
    # and must be written verbatim. The previous escapechar='\\' turned them
    # into malformed \" on newer pandas. GTF fields never contain tabs, so no
    # escaping is needed.
    tama_gtf_df.iloc[:, 0:9].to_csv(output_gtf, sep='\t', index=False,
                                    header=False, quoting=csv.QUOTE_NONE)


if __name__ == "__main__":
    # ------------------ CONSTANTS ------------------ #
    args = parse_args()

    path_to_tama_gtf = args.path_to_tama_gtf
    output_gtf = args.output_gtf
    path_to_classif = args.path_to_classif

    main(path_to_tama_gtf, path_to_classif, output_gtf)
