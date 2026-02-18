import sys
import os
import os.path
from Bio import SeqIO
import gzip
import pandas as pd

fastq_dir = sys.argv[1]
output_file = sys.argv[2]
adapter_dimer_df = pd.DataFrame(
    columns=['nr dimers', 'nr reads', 'dimer percentage'])

for fastq in os.listdir(fastq_dir):
    if fastq.endswith("R1.fastq.gz") and fastq.startswith('uf'):
        sample = fastq[:-12]
        fastq_filepath = os.path.join(fastq_dir, fastq)
        dimer_counter = 0
        nr_reads = 0
        with gzip.open(fastq_filepath, 'rt') as handle:
            for record in SeqIO.parse(handle, "fastq"):
                nr_reads += 1
                if record.seq[9:].startswith('AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC'):
                    dimer_counter += 1
            dimer_perc = dimer_counter/nr_reads
            adapter_dimer_df.loc[sample] = [
                dimer_counter, nr_reads, dimer_perc]


adapter_dimer_df.to_csv(output_file)
