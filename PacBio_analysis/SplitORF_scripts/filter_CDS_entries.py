import sys
from typing import Dict, List
from pygtftk.gtf_interface import GTF

def get_transcript_string(transcript_ids: List[str]) -> str:
    '''get string of transcript or other ids for filtering 
    of a GTF class object from pygtftk'''
    transcript_string = ''
    if len(transcript_ids) > 1:
        for transcript_id in transcript_ids[:-1]:
            transcript_string += transcript_id+','
    transcript_string += transcript_ids[-1]
    return transcript_string

orfanage_gtf = GTF(sys.argv[1], check_ensembl_format=False)
orf_cds = orfanage_gtf.select_by_key('feature', 'CDS')
tids = orf_cds.get_tx_ids(nr = True)
tid_string = get_transcript_string(tids)

orfanage_gtf.select_by_key('transcript_id', tid_string).write(sys.argv[2])