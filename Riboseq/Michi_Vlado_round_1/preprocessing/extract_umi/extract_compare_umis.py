import os
import sys

from Bio import SeqIO
import gzip

# fastq_path = '/projects/splitorfs/work/Riboseq/Output/Michi_Vlado_round_1/preprocess/fastp'
# umi_trimmed_path = '/projects/splitorfs/work/Riboseq/Output/Michi_Vlado_round_1/preprocess/fastp/UMI_trimmed_custom'

fastq_path = sys.argv[1]
umi_trimmed_path = sys.argv[2]

files = os.listdir(fastq_path)
# Filter out only files (if needed)
files = [f for f in files if os.path.isfile(os.path.join(
    fastq_path, f)) and 'R1' in f and f.endswith('fastq.gz')]
print(files, flush=True)


for file in files:
    file = os.path.join(fastq_path, file)
    print(file, flush=True)
    read2_file = file.replace('R1', 'R2')
    # Open the FASTQ files
    with gzip.open(file, "rt") as fwd_file, \
            gzip.open(read2_file, 'rt') as rev_file:
        fwd_reads = SeqIO.parse(fwd_file, "fastq")
        rev_reads = SeqIO.parse(rev_file, "fastq")

        nr_paired_reads = 0
        nr_no_umi_match = 0
        nr_umi_match = 0
        umi_trimmed_fwd_reads = []
        umi_trimmed_rev_reads = []
        nr_wrong_umi_rev = 0
        nr_wrong_umi_forward = 0
        # Iterate through the paired reads
        for fwd, rev in zip(fwd_reads, rev_reads):
            nr_paired_reads = nr_paired_reads + 1
            # Extract UMI sequences (first 9 bases for forward, last 9 for reverse)
            fwd_umi = fwd.seq[:9]
            rev_umi = rev.seq[-9:]

            # assert fwd_umi[-1] != 'A'
            # assert rev_umi[0] != 'T'
            if fwd_umi[-1] == 'A':
                nr_wrong_umi_forward += 1
            if rev_umi[0] == 'T':
                nr_wrong_umi_rev += 1

            # create new header name with only fwd UMi attached: compliability with umi-tools for deduplicaiton
            fwd_id = fwd.id + '_' + fwd_umi  # + '_' + rev_umi
            rev_id = rev.id + '_' + fwd_umi  # + '_' + rev_umi

            # subset RPFs to trim UMIs
            fwd_RPF = fwd[9:]
            rev_RPF = rev[:-9]

            # Change read name to append the new header with UMIs
            fwd_RPF.id = fwd_id
            rev_RPF.id = rev_id

            fwd_RPF.description = fwd.description.split(' ')[1]
            rev_RPF.description = rev.description.split(' ')[1]

            umi_trimmed_fwd_reads.append(fwd_RPF)
            umi_trimmed_rev_reads.append(rev_RPF)

            if fwd_umi == rev_umi.reverse_complement():
                nr_umi_match += 1
            else:
                nr_no_umi_match += 1
        print('nr paired reads', nr_paired_reads, flush=True)
        print('nr times UMI matches', nr_umi_match, flush=True)
        print('nr of times UMI does not match', nr_no_umi_match, flush=True)

        print('Nr of UMI with A in 9th position, not supposed to:',
              nr_wrong_umi_forward)
        print('Nr of rev UMI with T in first position, not supposed to:',
              nr_wrong_umi_rev)

    forward_trimmed_file_name = umi_trimmed_path + '/' + \
        os.path.basename(file).split(
            '.')[0] + '.' + 'R1.UMI_adapter_trimmed.fastq.gz'
    rev_trimmed_file_name = forward_trimmed_file_name.replace('R1', 'R2')
    with gzip.open(forward_trimmed_file_name, "wt") as out_handle:
        SeqIO.write(umi_trimmed_fwd_reads, out_handle, "fastq")

    with gzip.open(rev_trimmed_file_name, "wt") as out_handle:
        SeqIO.write(umi_trimmed_rev_reads, out_handle, "fastq")
