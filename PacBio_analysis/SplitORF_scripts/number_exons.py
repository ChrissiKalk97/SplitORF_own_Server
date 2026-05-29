#number_exons.py in.gtf out.gtf numbers the exons in in.gtf and wirtes
#the gtf with numbered exons to out.gtf

import sys
from pygtftk.gtf_interface import GTF
orfanage_gtf = GTF(sys.argv[1], check_ensembl_format=False)
orfanage_gtf.add_exon_number(key = "exon_number").write(sys.argv[2])