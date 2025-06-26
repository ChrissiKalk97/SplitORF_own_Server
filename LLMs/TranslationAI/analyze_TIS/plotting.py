import os
import seaborn as sbn
import matplotlib.pyplot as plt
import pandas as pd


def plot_so_background(df_merged, background_transcripts_df, nr_start, datatype, outdir, so='SO'):
    # HISTOGRAM ################################################################
    plt.figure(figsize=(8, 6))
    sbn.histplot(df_merged[f'{nr_start}_start_SO'], color="blue",
                 label=f"{so} {nr_start} start", bins=10, alpha=0.6)
    sbn.histplot(background_transcripts_df[f'{nr_start}_start'],
                 color="red", label=f"Non-{so} {nr_start} start", bins=10, alpha=0.6)
    # Labels and legend
    plt.xlabel(f"Probability of the {nr_start} start site")
    plt.ylabel("Frequency")
    plt.title(
        f"Distribution of {nr_start} Start Probabilities for {so} and non-{so} {datatype} transcripts")
    plt.legend()
    # Save plot
    plt.savefig(f"{outdir}/{datatype}_{nr_start}_start_{so}_non_{so}_historgam.svg",
                dpi=300, bbox_inches="tight")
    plt.close()

    # DENSITY PLOT #############################################################
    plt.figure(figsize=(8, 6))
    sbn.kdeplot(df_merged[f'{nr_start}_start_SO'], color="blue",
                label=f"{so} {nr_start} start", cut=0, fill=True)
    sbn.kdeplot(background_transcripts_df[f'{nr_start}_start'],
                color="red", label=f"Non-{so} {nr_start} start", cut=0, fill=True)

    # Labels and legend
    plt.xlabel(f"Probability of the {nr_start} start site")
    plt.ylabel("Density")
    plt.title(
        f"Distribution of {nr_start} Start Probabilities for {so} and non-{so} {datatype} transcripts")
    plt.legend()

    # Save plot
    plt.savefig(f"{outdir}/{datatype}_{nr_start}_start_{so}_non_{so}_density.svg",
                dpi=300, bbox_inches="tight")
    plt.close()


def plot_validated_so_probs_TransAI(Ribo_df_merged_first, Ribo_df_merged_second, outdir, datatype, sample):
    plt.figure(figsize=(8, 6))
    sbn.histplot(Ribo_df_merged_first['OrfProb'], color="blue",
                 label="First val ORFs", bins=10, alpha=0.6)
    sbn.histplot(Ribo_df_merged_second['OrfProb'], color="red",
                 label="Second val ORFs", bins=10, alpha=0.6)

    # Labels and legend
    plt.xlabel("Probability of first and second Riboseq val ORFs")
    plt.ylabel("Frequency")
    plt.title(
        f"Distribution of first and second Riboseq val ORF probs - {datatype}")
    plt.legend()

    plt.savefig(f"{outdir}/{datatype}_{sample}_Riboseq_empirical_first_and_second_ORF.svg",
                dpi=300, bbox_inches="tight")
    plt.close()


def plot_three_category_pie(cat1, cat2, cat3, nr_total_so, names_list, title, out_dir, figname, region_type, color_order):
    # plot pie chart with numbers: first, middle, last ORF
    fig, ax = plt.subplots()
    ax.pie([cat1, cat2, cat3],
           labels=names_list,
           autopct=lambda p: '{:.0f}'.format(p * nr_total_so / 100),
           colors=color_order)
    plt.title(f'{title} - {region_type}')

    plt.savefig(
        os.path.join(out_dir, f'{region_type}_{figname}.png'))
    plt.close()


def plot_val_so_sets(nr_orfs_with_UR, nr_validated_so, total_nr_so, outdir, region_type):
    nr_so_not_validated = nr_orfs_with_UR - nr_validated_so
    nr_so_no_UR = total_nr_so - nr_orfs_with_UR

    assert nr_validated_so + nr_so_not_validated + nr_so_no_UR == total_nr_so

    plot_three_category_pie(nr_validated_so,
                            nr_so_not_validated,
                            nr_so_no_UR,
                            total_nr_so,
                            ['# validated SO', '# SO with UR not validated',
                                '# SO without UR'],
                            'Split-ORF validation pie chart',
                            outdir,
                            'SO_validation_pie_chart',
                            region_type,
                            ['#CC79A7', '#FFC500', '#75C1C5']
                            )


def plot_start_prob_by_orf_position(so_df, title, region_type, out_dir, hue='OrfPosition', palette=['#75C1C5', '#CC79A7', '#FFC500']):
    plt.figure(figsize=(8, 6))
    sbn.histplot(data=so_df, x="OrfProb",
                 hue=hue, palette=palette, multiple='stack')

    # Labels and legend
    plt.xlabel("Probability of SO start sites")
    plt.ylabel("Frequency")
    plt.title(
        f"{title} - {region_type}")

    plt.savefig(f"{out_dir}/{region_type}_Riboseq_empirical_SO_start_probs_by_position.svg",
                dpi=300, bbox_inches="tight")
    plt.close()


def plot_so_start_probs_by_position_box(start_probs_position_info_df, title, region_type, out_dir, palette=['#CC79A7', '#FFC500', '#75C1C5']):
    plt.figure(figsize=(8, 6))
    sbn.boxplot(data=start_probs_position_info_df, x='OrfPosition', y="OrfProb",
                hue='SO_class', palette=palette, order=['first', 'middle', 'last'])

    # Labels and legend
    plt.xlabel("ORF position in transcript")
    plt.ylabel("Probability of SO start sites")
    plt.title(
        f"{title} - {region_type}")

    plt.legend(
        bbox_to_anchor=(1.05, 1),  # shift to the right
        loc='upper left',
        borderaxespad=0.
    )

    plt.tight_layout()

    plt.savefig(f"{out_dir}/{region_type}_boxplot_SO_set_start_probs_by_ORF_position.png",
                dpi=300, bbox_inches="tight")
    plt.show()
    plt.close()

#################################################################################
# ------------------ RiboTISH plotting                       ------------------ #
#################################################################################


def plot_RiboTISH_TransAI(URs_RiboTISH_TransAI_df_first, URs_RiboTISH_TransAI_df_second, datatype, outdir, sample):
    plt.figure(figsize=(8, 6))
    sbn.histplot(URs_RiboTISH_TransAI_df_first['OrfProb'], color="blue",
                 label="First val ORFs", bins=10, alpha=0.6)
    sbn.histplot(URs_RiboTISH_TransAI_df_second['OrfProb'], color="red",
                 label="Second val ORFs", bins=10, alpha=0.6)

    # Labels and legend
    plt.xlabel("Probability of first and second RiboTISH val ORFs")
    plt.ylabel("Frequency")
    plt.title(
        f"Distribution of first and second RiboTISH val ORF probs - {datatype}")
    plt.legend()

    plt.savefig(f"{outdir}/RiboTISH/{datatype}_{sample}_RiboTISH_first_and_second_ORF.svg",
                dpi=300, bbox_inches="tight")
    plt.close()


def plot_RiboTISH_inframecount_vs_probs(RiboTISH_TransAI_df, datatype, outdir, sample):
    RiboTISH_TransAI_df['RiboPValue'] = RiboTISH_TransAI_df['rt_name'].apply(
        lambda x: float(x.split('|')[4]))
    RiboTISH_TransAI_df['InFrameCount'] = RiboTISH_TransAI_df['rt_name'].apply(
        lambda x: int(x.split('|')[5]))

    plt.figure(figsize=(8, 6))
    sbn.scatterplot(x=RiboTISH_TransAI_df['InFrameCount'],
                    y=RiboTISH_TransAI_df['OrfProb'], label='RiboTISH val ORFs')

    # Labels and legend
    plt.xlabel("Inframecount")
    plt.ylabel("Probability of RiboTISH val ORFs")
    plt.title(
        f"RiboTISH predicted probs vs Inframecounts - {datatype}")
    plt.legend()

    plt.savefig(f"{outdir}/RiboTISH/{datatype}_{sample}_RiboTISH_inframecount_vs_prob.svg",
                dpi=300, bbox_inches="tight")
    plt.close()

    plt.figure(figsize=(8, 6))
    sbn.scatterplot(x=RiboTISH_TransAI_df['RiboPValue'],
                    y=RiboTISH_TransAI_df['OrfProb'], label='RiboTISH val ORFs')

    # Labels and legend
    plt.xlabel("RiboPValue")
    plt.ylabel("Probability of RiboTISH val ORFs")
    plt.title(
        f"RiboTISH predicted probs vs RiboPvalue - {datatype}")
    plt.legend()

    plt.savefig(f"{outdir}/RiboTISH/{datatype}_{sample}_RiboTISH_RiboPValue_vs_prob.svg",
                dpi=300, bbox_inches="tight")
    plt.close()
