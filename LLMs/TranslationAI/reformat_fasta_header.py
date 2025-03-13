import sys
from Bio import SeqIO
from Bio.Seq import Seq

fasta_file = sys.argv[1]
outfasta = sys.argv[2]

records_to_keep = []
with open(fasta_file) as input_handle, open (outfasta, 'w') as output_handle:
        for record in SeqIO.parse(input_handle, "fasta"):
                record_items = record.id.split('|')
                if record_items[3] == '1':
                    record_items[3] = '+'
                elif record_items[3] == '-1':
                    record_items[3] = '-'
                record.id = record_items[0] + ':' + record_items[1] + '-' + record_items[2] + '('+  record_items[3] + ')' + '(' + record_items[4] + '|' + record_items[5] + ')' + '(10, 100,)'
                record.description = ''
                record.seq = Seq(str(record.seq).replace('\n', ''))
                seq_string = str(record.seq).replace('\n', '')
                output_handle.write(f">{record.id}\n{seq_string}\n")

                


