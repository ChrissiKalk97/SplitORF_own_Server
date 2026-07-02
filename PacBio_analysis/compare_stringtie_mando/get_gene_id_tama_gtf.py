# ----- This script changes the order of transcript_id and gene_id in the 9th ----- #
# ----- field of a GTF that was generated using SQANTI3 to gene_id; transcript_id  ----- #

import os
import pandas as pd
import argparse
import csv


# path_to_tama_gtf = "/projects/splitorfs/work/PacBio/merged_bam_files/merge_mando_stringtie_isoquant_rescue_23_June_2026/HUVEC/HUVEC_merged_tama.gtf"
# output_gtf = "/projects/splitorfs/work/PacBio/merged_bam_files/merge_mando_stringtie_isoquant_rescue_23_June_2026/HUVEC/HUVEC_merged_tama_gene_id_02_07_26.gtf"
# path_to_classif = "/projects/splitorfs/work/PacBio/merged_bam_files/merge_mando_stringtie_isoquant_rescue_23_June_2026/SQANTI3_QC/HUVEC/isoforms_classification.txt"
# path_to_reference_gtf = "/projects/splitorfs/work/reference_files/filtered_Ens_reference_correct_29_09_25/Ensembl_110_filtered_equality_and_tsl1_2_correct_29_09_25.gtf"

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

    # parser.add_argument("path_to_reference_gtf",
    #                     help="Path to output GTF")

    return parser.parse_args()


def main(path_to_tama_gtf, path_to_classif, output_gtf):  # , path_to_reference_gtf
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
        then note as ambiguous and want to keep TAMA ID also if only novelGene is
        associated keep the TAMA ID"""
        vc = s.value_counts()
        if vc.index.str.startswith('novelGene').all():
            return 'all_novel'
        if vc.size == 1:
            return 'clean'
        return 'dominant' if vc.iloc[0] / vc.sum() >= 0.8 else 'ambiguous'

    buckets = classif_df.groupby('gene_id')['associated_gene'].agg(bucket)
    print(buckets.value_counts())

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

    unmapped = tama_gtf_df['gene_id'].isna()
    if unmapped.any():
        print(f"Warning: {int(unmapped.sum())} feature line(s) had no SQANTI "
              f"association; keeping native TAMA gene_id.")
        tama_gtf_df.loc[unmapped, 'gene_id'] = \
            tama_gtf_df.loc[unmapped, 'native_gene']

    # split gene_ids that landed on more than one (chrom, strand) locus
    loci = tama_gtf_df[['gene_id', 0, 6]].drop_duplicates()
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

    # --- collapse duplicate 'gene' records that share a gene_id on one locus ---
    # After disambiguate, a shared gene_id means one locus (multi-strand/chrom
    # were already suffixed). Give every gene line the full min..max span of its
    # gene_id, then keep only one gene line per id. transcript/exon lines are
    # untouched, so counts are unaffected. MUST run after disambiguate.
    is_gene = tama_gtf_df[2] == 'gene'
    span = tama_gtf_df[is_gene].groupby('gene_id').agg(start=(3, 'min'),
                                                       end=(4, 'max'))
    tama_gtf_df.loc[is_gene, 3] = \
        tama_gtf_df.loc[is_gene, 'gene_id'].map(span['start']).values
    tama_gtf_df.loc[is_gene, 4] = \
        tama_gtf_df.loc[is_gene, 'gene_id'].map(span['end']).values
    dup_gene = is_gene & tama_gtf_df.duplicated(subset=[2, 'gene_id'],
                                                keep='first')
    n_dup = int(dup_gene.sum())
    if n_dup:
        print(
            f"Collapsed {n_dup} duplicate gene record(s) into spanning lines.")
    tama_gtf_df = tama_gtf_df[~dup_gene]

    tama_gtf_df.iloc[:, 8] = tama_gtf_df.apply(
        lambda row: change_gene_id(row), axis=1)

    # QUOTE_NONE only: the attribute field legitimately contains double quotes
    # and must be written verbatim (escapechar='\\' would mangle them).
    tama_gtf_df.iloc[:, 0:9].to_csv(output_gtf, sep='\t', index=False,
                                    header=False, quoting=csv.QUOTE_NONE)


if __name__ == "__main__":
    args = parse_args()
    main(args.path_to_tama_gtf, args.path_to_classif,
         args.output_gtf)  # , args.path_to_reference_gtf
