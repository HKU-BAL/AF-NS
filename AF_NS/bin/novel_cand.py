import os
import sys
from Bio import SeqIO

path = sys.argv[1]
len_cutoff = int(sys.argv[2])


names = set()
f = open(path + '3rd.paf','r')
for line in f.readlines():
    line = line.strip().split('\t')
    names.add(line[0])
f.close()

f_write = open(path + '/novel_tmp.fa','w')
    
num = 0
for record in SeqIO.parse(path + '/clean.fa','fasta'):
    if str(record.id) not in names and len(str(record.seq)) >= len_cutoff:
       f_write.write('>' + str(num)  + '\n')
       f_write.write(str(record.seq) + '\n')
       num += 1
f_write.close()
   



