import sys
import pandas as pd
from pygtftk.gtf_interface import GTF


# map tids to their respective gids
def main():

    def obtain_gid_tid_dict(gtf_path):
        gtf = GTF(gtf_path, check_ensembl_format=False)
        gid_tid_dict = gtf.select_by_key("feature", "transcript").extract_data(
            keys="gene_id,transcript_id", as_dict_of_values=True)
        assert len(gid_tid_dict.keys()) == len(set(gid_tid_dict.keys()))
        gid_tid_dict = {k.split('.')[0]: v for (k, v) in gid_tid_dict.items()}
        return gid_tid_dict

    def map_tid_txt_to_gid_txt(gid_txt_file, gid_tid_dict, outfile_path, map_csv):
        tid_series = pd.read_csv(gid_txt_file, sep=";", header=0)
        tid_series = tid_series.iloc[1:, :]
        tid_series.iloc[:, 0] = tid_series.iloc[:, 0].apply(
            lambda x: x.split('.')[0])
        tid_series['tid'] = tid_series.iloc[:, 0].map(gid_tid_dict)
        tid_series.to_csv(map_csv)
        with open(outfile_path, 'w') as outfile:
            for tid in tid_series['tid'].dropna():
                outfile.write(tid + '\n')

    gtf_path = sys.argv[1]
    gid_txt_file = sys.argv[2]
    outfile_path = sys.argv[3]
    map_csv = sys.argv[4]
    gid_tid_dict = obtain_gid_tid_dict(gtf_path)
    map_tid_txt_to_gid_txt(gid_txt_file, gid_tid_dict, outfile_path, map_csv)


if __name__ == "__main__":
    main()
