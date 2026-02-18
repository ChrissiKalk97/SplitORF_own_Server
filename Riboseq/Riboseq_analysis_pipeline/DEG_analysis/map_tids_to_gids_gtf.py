import sys
import pandas as pd
from pygtftk.gtf_interface import GTF


# map tids to their respective gids
def main():

    def obtain_tid_gid_dict(gtf_path):
        gtf = GTF(gtf_path, check_ensembl_format=False)
        tid_gid_dict = gtf.select_by_key("feature", "transcript").extract_data(
            keys="transcript_id,gene_id", as_dict_of_values=True)
        assert len(tid_gid_dict.keys()) == len(set(tid_gid_dict.keys()))
        return tid_gid_dict

    def map_tid_txt_to_gid_txt(tid_txt_file, tid_gid_dict, outfile_path):
        tid_series = pd.read_csv(tid_txt_file, header=None)
        tid_series['gid'] = tid_series.loc[:, 0].map(tid_gid_dict)
        with open(outfile_path, 'w') as outfile:
            for gid in tid_series['gid']:
                if '.' in gid:
                    gid = gid.split('.')[0]
                outfile.write(gid + '\n')

    gtf_path = sys.argv[1]
    tid_txt_file = sys.argv[2]
    outfile_path = sys.argv[3]
    tid_gid_dict = obtain_tid_gid_dict(gtf_path)
    map_tid_txt_to_gid_txt(tid_txt_file, tid_gid_dict, outfile_path)


if __name__ == "__main__":
    main()
