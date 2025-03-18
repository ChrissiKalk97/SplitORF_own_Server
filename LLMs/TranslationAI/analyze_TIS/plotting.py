import seaborn as sbn
import matplotlib.pyplot as plt
import pandas as pd


def plot_SO_background(df_merged, background_transcripts_df, nr_start, datatype, outdir):
    plt.figure(figsize=(8, 6))
    sbn.histplot(df_merged[f'{nr_start}_start_SO'], color="blue",
                 label=f"SO {nr_start} start", bins=10, alpha=0.6)
    sbn.histplot(background_transcripts_df[f'{nr_start}_start'],
                 color="red", label=f"Non-SO {nr_start} start", bins=10, alpha=0.6)

    # Labels and legend
    plt.xlabel(f"Probability of the {nr_start} start site")
    plt.ylabel("Frequency")
    plt.title(
        f"Distribution of {nr_start} Start Probabilities for SO and non-SO {datatype} transcripts")
    plt.legend()

    # Save plot
    plt.savefig(f"{outdir}/{datatype}_{nr_start}_start_SO_non_SO_historgam.svg",
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

    plt.savefig(f"{outdir}/{datatype}_{sample}_RiboTISH_first_and_second_ORF.svg",
                dpi=300, bbox_inches="tight")
    plt.close()
