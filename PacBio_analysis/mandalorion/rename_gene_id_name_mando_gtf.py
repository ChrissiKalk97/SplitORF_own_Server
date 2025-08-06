# Rename the Mandolorion gene_id and gene_name attributes
# these are currently gene_id_gene_name, want to have them seperately
# so gene_id for "gene_id" and gene_name for "gene_name"

# usage: python rename_gene_id_name_mando_gtf.py Mando.gtf Output.gtf


import re
import argparse


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Rename the Mandolorion gene_id and gene_name attributes.")

    parser.add_argument(
        "mando_gtf_file",
        type=str,
        help="Path to the input Mando GTF file"
    )

    parser.add_argument(
        "output_gtf_file",
        type=str,
        help="Path to the output GTF file"
    )

    return parser.parse_args()


def change_MSTRG_gene_name(mando_gtf_file, output_gtf_file):
    """
    process GTF file from Mandalorion and write to a new GTF file
    change the gene_id to gene_id only and gene_name to gene_name only
    """
    mando_gtf = open(mando_gtf_file, 'r')
    output_gtf = open(output_gtf_file, 'w')
    for l in mando_gtf:
        if not l.startswith('#'):
            # info = re.split(r'\s+|;', l)

            l = re.sub(r'(gene_id "ENSG[0-9]*)(_[^;]*)"', r'\1"', l)
            l = re.sub(r'(gene_name )("ENSG[0-9]*_)([^;]*)"', r'\1"\3"', l)
        # Write the feature to the output GTF file
        output_gtf.write(l)


def main(mando_gtf_file, output_gtf_file):
    change_MSTRG_gene_name(mando_gtf_file, output_gtf_file)


if __name__ == "__main__":
    args = parse_arguments()
    main(args.mando_gtf_file, args.output_gtf_file)
