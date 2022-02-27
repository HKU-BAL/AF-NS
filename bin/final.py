import os
import sys
from Bio import SeqIO
import collections

path = sys.argv[1]
repeat = collections.defaultdict(int)
repeat_percent = float(sys.argv[2])

f = open(path + '/tmp/repeat.bed','r')
for line in f.readlines():
    line = line.strip().split('\t')
    repeat[line[0]] += int(line[2]) - int(line[1])
f.close()

f_write = open(path + '/novel.fa','w')
for record in SeqIO.parse(path + '/tmp/cluster.fa','fasta'):
    if 1.0*repeat[str(record.id)]/ len(str(record.seq))<repeat_percent:
       f_write.write( ">"+str(record.id) + '\n')
       f_write.write( str(record.seq )+ '\n')
f_write.close()

