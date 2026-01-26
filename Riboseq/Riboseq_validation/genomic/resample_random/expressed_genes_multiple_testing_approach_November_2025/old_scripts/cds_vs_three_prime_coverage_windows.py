import os
import pandas as pd
import argparse
import seaborn as sbn
import matplotlib.pyplot as plt
import numpy as np
# set_config(transform_output='pandas')


def parse_args():
    parser = argparse.ArgumentParser(
        description='.'
    )

    # Required positional arguments
    parser.add_argument('three_prime_coverage_file',
                        help='Path to three prime coverage file')

    parser.add_argument('cds_coverage_file',
                        help='Path to three cds coverage file')

    parser.add_argument('sample',
                        help='sample')

    parser.add_argument('three_prime_original_for_strand_info',
                        help='Path to 3 prime bed with strand information')

    parser.add_argument('alpha',
                        help='alpha')

    return parser.parse_args()


def main(three_prime_coverage_file, sample, cds_coverage_file, three_prime_original_for_strand_info, alpha):
    outdir = os.path.dirname(three_prime_coverage_file)

    three_prime_coverage_file_df = pd.read_csv(
        three_prime_coverage_file, sep='\t', header=None,
        names=['chr', 'start', 'stop', '3_prime_name', 'nr_overlap_reads',
               'nr_bases_covered', 'length_3_prime_UTR', 'covered_fraction'],
        low_memory=False)

    # length normalization of reads
    three_prime_coverage_file_df['RPK'] = three_prime_coverage_file_df['nr_overlap_reads'] / \
        (three_prime_coverage_file_df['length_3_prime_UTR']/1000)

    cds_coverage_file_df = pd.read_csv(
        cds_coverage_file, sep='\t', header=None,
        names=['chr', 'start', 'stop', '3_prime_name', 'nr_overlap_reads',
               'nr_bases_covered', 'length_cds', 'covered_fraction'],
        low_memory=False)

    # length normalization of reads
    cds_coverage_file_df['RPK'] = cds_coverage_file_df['nr_overlap_reads'] / \
        (cds_coverage_file_df['length_cds']/1000)

    # plot RPK distribution CDS hist
    ax = sbn.histplot(
        cds_coverage_file_df, x='RPK', bins=35000)
    ax.set_xlim(0, 2000)
    ax.set_title(f'CDS RPK for {sample}')
    plt.show()

    fig = ax.get_figure()
    fig.savefig(os.path.join(outdir,
                f'{sample}_CDS_RPK_hist_0_2000.png'), dpi=300, bbox_inches='tight')
    plt.close(fig)

    cds_coverage_file_df['log_RPK'] = np.log10(cds_coverage_file_df['RPK'] + 1)

    # plot RPK distribution CDS box
    ax = sbn.boxplot(
        cds_coverage_file_df, x='log_RPK')
    ax.set_title(f'log CDS RPK for {sample}')
    plt.show()

    fig = ax.get_figure()
    fig.savefig(os.path.join(outdir,
                f'{sample}_CDS_log_RPK_box.png'), dpi=300, bbox_inches='tight')
    plt.close(fig)

    # plot RPK distribution CDS hist
    ax = sbn.histplot(
        cds_coverage_file_df, x='log_RPK', bins=150)
    ax.set_title(f'Log CDS RPK for {sample}')
    plt.show()

    fig = ax.get_figure()
    fig.savefig(os.path.join(outdir,
                f'{sample}_CDS_log_RPK_hist.png'), dpi=300, bbox_inches='tight')
    plt.close(fig)

    # plot RPK distribution 3primes hist
    ax = sbn.histplot(
        three_prime_coverage_file_df[three_prime_coverage_file_df['RPK'] < 2000],
        x='RPK',
        bins=15000)
    ax.set_xlim(0, 200)
    ax.set_title(f'3 prime RPK for {sample}')
    plt.show()

    fig = ax.get_figure()
    fig.savefig(os.path.join(outdir,
                f'{sample}_three_primes_RPK_hist_0_200.png'), dpi=300, bbox_inches='tight')
    plt.close(fig)

    ax = sbn.histplot(
        three_prime_coverage_file_df[(three_prime_coverage_file_df['RPK'] > 0) &
                                     (three_prime_coverage_file_df['RPK'] < 2000)],
        x='RPK',
        bins=15000)
    ax.set_xlim(0, 200)
    ax.set_title(f'filtered 3 prime RPK for {sample}')
    plt.show()

    fig = ax.get_figure()
    fig.savefig(os.path.join(outdir,
                f'{sample}_three_primes_RPK_hist_greater_0_0_200.png'),
                dpi=300, bbox_inches='tight')
    plt.close(fig)

    three_prime_coverage_file_df['log_RPK'] = np.log10(
        three_prime_coverage_file_df['RPK'] + 1)

    # plot RPK distribution
    ax = sbn.boxplot(
        three_prime_coverage_file_df, x='log_RPK')
    ax.set_title(f'log 3 prime RPK for {sample}')
    plt.show()

    fig = ax.get_figure()
    fig.savefig(os.path.join(outdir,
                f'{sample}_three_primes_log_RPK_box.png'), dpi=300, bbox_inches='tight')
    plt.close(fig)

    # print('CDS 25th quantile', cds_coverage_file_df['RPK'].quantile(0.25))
    # print('CDS median', cds_coverage_file_df['RPK'].quantile(0.5))

    # print('three prime 25th quantile',
    #       three_prime_coverage_file_df['RPK'].quantile(0.25))
    # print('three prime median',
    #       three_prime_coverage_file_df['RPK'].quantile(0.5))
    # print('three prime 75th quantile',
    #       three_prime_coverage_file_df['RPK'].quantile(0.75))

    # FILTER OUT VALUES THAT ARE 0
    from scipy import stats
    cds_coverage_file_df_filtered = cds_coverage_file_df[cds_coverage_file_df['RPK'] > 0]
    distributions = [stats.skewnorm, stats.weibull_min,
                     stats.gamma, stats.norm, stats.t]
    fits = {}

    for dist in distributions:
        params = dist.fit(cds_coverage_file_df_filtered['log_RPK'])
        loglik = np.sum(dist.logpdf(
            cds_coverage_file_df_filtered['log_RPK'], *params))
        k = len(params)
        aic = 2*k - 2*loglik
        bic = np.log(
            len(cds_coverage_file_df_filtered['log_RPK']))*k - 2*loglik
        fits[dist.name] = {"params": params,
                           "AIC": aic, "BIC": bic, "loglik": loglik}

        fig = plt.figure()
        plt.hist(cds_coverage_file_df_filtered['log_RPK'],
                 bins=100, density=True, alpha=0.5, label='Data')

        x = np.linspace(0, max(cds_coverage_file_df_filtered['log_RPK']), 1000)
        pdf_fitted = dist.pdf(x, *params)
        plt.plot(x, pdf_fitted, 'r-', lw=2, label=f'{dist.name} fit')

        plt.xlabel('Value')
        plt.ylabel('Density')
        plt.title(f'Histogram with {dist.name} fit - {sample}')
        plt.legend()
        plt.savefig(os.path.join(outdir, f'{sample}_{dist.name}_fit.png'))
        # plt.show()
        plt.close(fig)

        fig = plt.figure()
        plt.hist(cds_coverage_file_df_filtered['log_RPK'],
                 bins=50, density=True, alpha=0.5, label='Data')

        x = np.linspace(0, max(cds_coverage_file_df_filtered['log_RPK']), 1000)
        pdf_fitted = dist.pdf(x, *params)
        plt.plot(x, pdf_fitted, 'r-', lw=2, label=f'{dist.name} fit')

        plt.xlabel('Value')
        plt.ylabel('Density')
        plt.title(f'Histogram with {dist.name} fit - {sample}')
        plt.legend()
        plt.savefig(os.path.join(
            outdir, f'{sample}_{dist.name}_50_bins_fit.png'))
        # plt.show()
        plt.close(fig)
    print(fits)

    fig = plt.figure()
    stats.probplot(cds_coverage_file_df_filtered['log_RPK'], dist="gamma",
                   sparams=fits['gamma']['params'], plot=plt)
    plt.title(f"{sample} – Gamma QQ plot")
    plt.savefig(os.path.join(
        outdir, f'{sample}_gamma_QQ_plot.png'), dpi=300, bbox_inches='tight')
    plt.close(fig)

    fig = plt.figure()
    stats.probplot(cds_coverage_file_df_filtered['log_RPK'], dist="skewnorm",
                   sparams=fits['skewnorm']['params'], plot=plt)
    plt.title(f"{sample} – Skew-normal QQ plot")
    plt.savefig(os.path.join(
        outdir, f'{sample}_skewnorm_QQ_plot.png'), dpi=300, bbox_inches='tight')
    plt.close(fig)

    fig = plt.figure()
    stats.probplot(cds_coverage_file_df_filtered['log_RPK'], dist="t",
                   sparams=fits['t']['params'], plot=plt)
    plt.title(f"{sample} – t QQ plot")
    plt.savefig(os.path.join(
        outdir, f'{sample}_t_QQ_plot.png'), dpi=300, bbox_inches='tight')
    plt.close(fig)

    fig = plt.figure()
    stats.probplot(cds_coverage_file_df_filtered['log_RPK'], dist="weibull_min",
                   sparams=fits['weibull_min']['params'], plot=plt)
    plt.title(f"{sample} – Weibull-min QQ plot")
    plt.savefig(os.path.join(
        outdir, f'{sample}_weibull_min_QQ_plot.png'), dpi=300, bbox_inches='tight')
    plt.close(fig)

    lowest_BIC = np.inf
    for dist, values in fits.items():
        if values['BIC'] < lowest_BIC:
            best_dist_bic = dist
            lowest_BIC = values['BIC']

    print(f'best distribution in terms of BIC is {best_dist_bic}')
    print(fits[best_dist_bic])

    lowest_AIC = np.inf
    for dist, values in fits.items():
        if values['AIC'] < lowest_AIC:
            best_dist_aic = dist
            lowest_AIC = values['AIC']

    print(f'best distribution in terms of AIC is {best_dist_aic}')
    print(fits[best_dist_aic])

    # n = len(cds_coverage_file_df_filtered['log_RPK'])
    # loc = fits[best_dist_bic]['params'][1]
    # df = fits[best_dist_bic]['params'][0]
    # scale = fits[best_dist_bic]['params'][2]

    # t quantile for two-sided PI
    # t_quantile = stats.t.ppf(alpha, df)

    # Prediction interval
    # only interested in lower
    # PI_lower = loc + t_quantile * scale * np.sqrt(1 + 1/n)
    dist = getattr(stats, best_dist_bic)
    PI_lower = dist.ppf(alpha, *fits[best_dist_bic]['params'])

    # here we can keep the 0-counts as we anyway want to have them
    # select all the values below the PI_lower
    three_prime_coverage_file_df_filtered = three_prime_coverage_file_df[
        three_prime_coverage_file_df['log_RPK'] < PI_lower].copy()

    # how does the three prime RPK distribution look like after the filtering?
    print('Fitlering threshold for log RPK:', PI_lower)

    ax = sbn.histplot(
        three_prime_coverage_file_df_filtered[
            (three_prime_coverage_file_df_filtered['RPK'] > 0) &
            (three_prime_coverage_file_df_filtered['RPK'] < 2000)],
        x='RPK',
        bins=15000)
    ax.set_xlim(0, 200)
    ax.set_title(f'CDS filtered 3 prime RPK for {sample}')
    plt.show()

    fig = ax.get_figure()
    fig.savefig(os.path.join(outdir,
                f'{sample}_three_primes_RPK_after_CDS_filter_hist_greater_0_0_200.png'), dpi=300, bbox_inches='tight')
    plt.close(fig)

    # need to recover the strand information, based on gene
    strand_info_df = pd.read_csv(three_prime_original_for_strand_info, sep='\t', names=[
                                 'chr', 'start', 'stop', 'name', 'score', 'strand'])
    strand_info_df['Gene'] = strand_info_df['name'].apply(
        lambda x: x.split('|')[0])
    strand_dict = dict(zip(strand_info_df['Gene'], strand_info_df['strand']))

    three_prime_coverage_file_df_filtered['Gene'] = three_prime_coverage_file_df_filtered['3_prime_name'].apply(
        lambda x: x.split('|')[0])
    three_prime_coverage_file_df_filtered['strand'] = three_prime_coverage_file_df_filtered['Gene'].map(
        strand_dict)
    three_prime_coverage_file_df_filtered['score'] = 0

    three_prime_coverage_file_df_filtered.loc[:, ['chr', 'start', 'stop', '3_prime_name', 'score', 'strand']].to_csv(os.path.join(
        outdir, f'3_primes_filtered_for_CDS_distribution_{sample}.bed'), sep='\t', header=False, index=False)


if __name__ == '__main__':
    # ------------------ CONSTANTS ------------------ #
    args = parse_args()

    three_prime_coverage_file = args.three_prime_coverage_file
    sample = args.sample
    cds_coverage_file = args.cds_coverage_file
    three_prime_original_for_strand_info = args.three_prime_original_for_strand_info
    alpha = float(args.alpha)

    # three_prime_coverage_file = '/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/resample_q10_expression_filter/NMD_genome/ERR3367798/3_primes_genomic_merged_numbered_ERR3367797_windows_coverage.tsv'
    # cds_coverage_file = '/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/resample_q10_expression_filter/NMD_genome/ERR3367798/Ens_110_CDS_coordinates_genomic_protein_coding_tsl_refseq_filtered_ERR3367797_windows_coverage.tsv'
    # sample = 'ERR3367798'
    # three_prime_original_for_strand_info = '/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/resample_q10_expression_filter/NMD_genome/ERR3367798/3_primes_genomic_merged_numbered_ERR3367797.bed'
    # alpha = 0.05

    # check ENSG00000144136|ENST00000272542|11291 for ERR8, jhas very high coverage

    main(three_prime_coverage_file, sample, cds_coverage_file,
         three_prime_original_for_strand_info, alpha)
