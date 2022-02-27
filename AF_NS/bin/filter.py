import os
import sys
from Bio import SeqIO

folder = sys.argv[1]
fasta_name = sys.argv[2]
bed_name = sys.argv[3]
out =  sys.argv[4]
len_cutoff = int(sys.argv[5])

seq = {}
for record in SeqIO.parse(folder + fasta_name,'fasta'):
    seq[str(record.id)] = str(record.seq)

f= open(folder + bed_name,'r')
f_write1 = open(folder + out + '.bed','w')
f_write = open(folder + out + '.fa','w')
for line in f.readlines():
    line = line.strip().split('\t')
    length = int(line[2]) - int(line[1])
    if length >= len_cutoff and line[0] in seq:
       f_write1.write(line[0]+ ':' + line[1] + '-'+line[2] + '\t0\t' + str(length) + '\n')
       f_write.write('>' +line[0]+ ':' + line[1] + '-'+line[2] + '\n')
       f_write.write(seq[line[0]][int(line[1]):int(line[2])] + '\n')
f.close()
f_write1.close()
f_write.close()




