import os
import sys
from Bio import SeqIO

input1 = sys.argv[1]
output1 = sys.argv[2]
f_write = open(output1,'w')
for record in SeqIO.parse(input1,'fasta'):
    f_write.write ( str(record.id) + '\t0\t' + str(len(str(record.seq))) + '\n')
f_write.close()
 

