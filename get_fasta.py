from Bio import SeqIO

import sys

file = sys.argv[1]

primer = str(list(SeqIO.parse(file,'fasta'))[0].seq)

print(primer)