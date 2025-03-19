import sys
import glob
import re
import os


import pandas as pd
import upsetplot
from upsetplot import from_contents
from pybedtools import BedTool
import matplotlib.pyplot as plt


TransAI_scripts_path = "/home/ckalk/scripts/SplitORFs/LLMs/TranslationAI/analyze_TIS/"
sys.path.append(os.path.dirname(TransAI_scripts_path))


RiboTISH_path = sys.argv[1]
datatype = sys.argv[2]


# RiboTISH_path = '/projects/splitorfs/work/Riboseq/Output/RiboTISH_NMD_custom'
# datatype = 'NMD'


# here comes the import
from analysis import create_RiboTISH_BedTool, obtain_correct_ORF_RiboTISH  # noqa: E402


def main():
    UR_BedTool = BedTool(
        f'/projects/splitorfs/work/Riboseq/data/region_input/genomic/Unique_DNA_regions_genomic_{datatype}_16_12_24_chrom_sorted.bed')

    RiboTISH_hits_dict = {}
    for RiboTISH_file in glob.glob(f"{RiboTISH_path}/*.csv"):
        sample = os.path.basename(
            RiboTISH_file).rsplit('_', 3)[0]
        print(sample)
        RiboTISH_SO_df = pd.read_csv(RiboTISH_file, header=0)
        if len(RiboTISH_SO_df.index) > 0:
            RiboTISH_BedTool, RiboTISH_SO_df = create_RiboTISH_BedTool(
                RiboTISH_SO_df)

            URs_found_in_RiboTISH_df = obtain_correct_ORF_RiboTISH(
                UR_BedTool, RiboTISH_BedTool)

            RiboTISH_hits_dict[sample] = URs_found_in_RiboTISH_df['name'].to_list(
            )

    Upsetplot_data = from_contents(RiboTISH_hits_dict)

    custom_order = ["OHMX20220060_001", "OHMX20220060_002", "OHMX20220060_003", "OHMX20220060_004",
                    "OHMX20220060_005", "OHMX20220060_006", "ERR3367797", "ERR3367798"]  # Define your preferred order

    # Convert index to categorical with custom order
    Upsetplot_data = Upsetplot_data.reorder_levels(custom_order)

    fig = plt.figure(figsize=(8, 6))  # Create a figure

    upsetplot.plot(Upsetplot_data, sort_categories_by='-input',
                   sort_by='-degree')

    # Save the plot as an SVG file
    plt.savefig(f"{RiboTISH_path}/plots/{datatype}_upset_plot.svg",
                format="svg", bbox_inches="tight")


if __name__ == "__main__":
    main()
