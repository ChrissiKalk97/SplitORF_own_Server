import os
import pandas as pd
import seaborn as sbn
import argparse
import matplotlib.pyplot as plt


def parse_args():
    parser = argparse.ArgumentParser(
        description='.'
    )

    # Required positional arguments
    parser.add_argument('path_to_mando_output',
                        help='Path to Mando Output')

    parser.add_argument('nr_samples',
                        help='Nr of samples used for assembly')

    parser.add_argument('x_max',
                        help='maximum x value')

    parser.add_argument('--quant_file_name',
                        default="Isoforms.filtered.clean.quant",
                        help='Name of the Mando isoform qunatifiation file')

    return parser.parse_args()


def main(path_to_mando_output, nr_samples, x_max, quant_file_name):

    sample_name = os.path.basename(path_to_mando_output)

    isoform_quant_path = os.path.join(path_to_mando_output, quant_file_name)
    mando_quant_df = pd.read_csv(isoform_quant_path, header=0, sep='\t')
    samples = [str(nr) for nr in range(1, int(nr_samples) + 1)]

    mando_quant_df['expression_sum'] = mando_quant_df[samples].sum(axis=1)

    # plot histogram of the summed expression values
    sbn.histplot(mando_quant_df['expression_sum'], bins=100000,
                 kde=False, color='skyblue')

    plt.xlim(0, x_max)
    plt.xlabel('expression of LRs per isoform')
    plt.ylabel('Count')
    plt.title(f'Expression distribution of LRs per {sample_name} isoform')
    plt.show()

    plt.savefig(os.path.join(path_to_mando_output,
                             f'{sample_name}_LR_expression_quant_mando_{x_max}.png'),
                dpi=300, bbox_inches='tight')

    plt.close()


if __name__ == '__main__':
    # ------------------ CONSTANTS ------------------ #
    args = parse_args()

    path_to_mando_output = args.path_to_mando_output
    nr_samples = args.nr_samples
    x_max = int(args.x_max)
    quant_file_name = args.quant_file_name

    main(path_to_mando_output, nr_samples, x_max, quant_file_name)
