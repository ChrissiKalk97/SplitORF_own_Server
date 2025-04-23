import sys
from Bio import Entrez

NCBI_cont_fasta = sys.argv[1]

Entrez.email = "kalk@med.uni-frankfurt.de"

accessions = ["NR_023379.1", "NR_146151.1", "NR_146144.1",
              "NR_145819.1", "NR_146117.1", "NR_003285.3", "NR_003286.4", "X12811.1", "NR_003287.4"]


handle = Entrez.efetch(
    db="nucleotide",


    id=",".join(accessions),
    rettype="fasta",
    retmode="text"
)
fasta_data = handle.read()
handle.close()


with open(f"{NCBI_cont_fasta}", "w") as f:
    f.write(fasta_data)
