import seaborn as sbn
import matplotlib.pyplot as plt
import pandas as pd


def plot_SO_background(df_merged, background_transcripts_df, nr_start, datatype, outdir, SO = 'SO'):
    # HISTOGRAM ################################################################
    plt.figure(figsize=(8, 6))
    sbn.histplot(df_merged[f'{nr_start}_start_SO'], color="blue",
                 label=f"{SO} {nr_start} start", bins=10, alpha=0.6)
    sbn.histplot(background_transcripts_df[f'{nr_start}_start'],
                 color="red", label=f"Non-{SO} {nr_start} start", bins=10, alpha=0.6)
    # Labels and legend
    plt.xlabel(f"Probability of the {nr_start} start site")
    plt.ylabel("Frequency")
    plt.title(
        f"Distribution of {nr_start} Start Probabilities for {SO} and non-{SO} {datatype} transcripts")
    plt.legend()
    # Save plot
    plt.savefig(f"{outdir}/{datatype}_{nr_start}_start_{SO}_non_{SO}_historgam.svg",
                dpi=300, bbox_inches="tight")
    plt.close()

    # DENSITY PLOT #############################################################
    plt.figure(figsize=(8, 6))
    sbn.kdeplot(df_merged[f'{nr_start}_start_SO'], color="blue",
                 label=f"{SO} {nr_start} start", cut = 0, fill = True)
    sbn.kdeplot(background_transcripts_df[f'{nr_start}_start'],
                 color="red", label=f"Non-{SO} {nr_start} start", cut = 0, fill = True)

    # Labels and legend
    plt.xlabel(f"Probability of the {nr_start} start site")
    plt.ylabel("Density")
    plt.title(
        f"Distribution of {nr_start} Start Probabilities for {SO} and non-{SO} {datatype} transcripts")
    plt.legend()

    # Save plot
    plt.savefig(f"{outdir}/{datatype}_{nr_start}_start_{SO}_non_{SO}_density.svg",
                dpi=300, bbox_inches="tight")
    plt.close()




def plot_emp_background_TransAI(Ribo_df_merged_first, Ribo_df_merged_second, outdir, datatype, sample):
    plt.figure(figsize=(8, 6))
    sbn.histplot(Ribo_df_merged_first['ORF_prob'], color="blue",
                 label="First val ORFs", bins=10, alpha=0.6)
    sbn.histplot(Ribo_df_merged_second['ORF_prob'], color="red",
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


def plot_RiboTISH_TransAI(URs_RiboTISH_TransAI_df_first, URs_RiboTISH_TransAI_df_second, datatype, outdir, sample):
    plt.figure(figsize=(8, 6))
    sbn.histplot(URs_RiboTISH_TransAI_df_first['ORF_prob'], color="blue",
                 label="First val ORFs", bins=10, alpha=0.6)
    sbn.histplot(URs_RiboTISH_TransAI_df_second['ORF_prob'], color="red",
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
    RiboTISH_TransAI_df['RiboPValue'] = RiboTISH_TransAI_df['rt_name'].apply(lambda x: float(x.split('|')[4]))
    RiboTISH_TransAI_df['InFrameCount'] = RiboTISH_TransAI_df['rt_name'].apply(lambda x: int(x.split('|')[5]))

    plt.figure(figsize=(8, 6))
    sbn.scatterplot(x = RiboTISH_TransAI_df['InFrameCount'], y= RiboTISH_TransAI_df['ORF_prob'], label = 'RiboTISH val ORFs')

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
    sbn.scatterplot(x = RiboTISH_TransAI_df['RiboPValue'], y= RiboTISH_TransAI_df['ORF_prob'], label = 'RiboTISH val ORFs')

    # Labels and legend
    plt.xlabel("RiboPValue")
    plt.ylabel("Probability of RiboTISH val ORFs")
    plt.title(
        f"RiboTISH predicted probs vs RiboPvalue - {datatype}")
    plt.legend()

    plt.savefig(f"{outdir}/RiboTISH/{datatype}_{sample}_RiboTISH_RiboPValue_vs_prob.svg",
                dpi=300, bbox_inches="tight")
    plt.close()


