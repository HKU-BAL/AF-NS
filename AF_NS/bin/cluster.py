import os
import sys
from Bio import SeqIO
import networkx as nx

path = sys.argv[1]
cov = float(sys.argv[2])

sequence = {}
for record in SeqIO.parse(path + '/novel_tmp.fa','fasta'):
    sequence[str(record.id)] = str(record.seq)

reads = []
f = open(path + 'overlap.paf','r')
for line in f.readlines():
    line = line.strip().split('\t')
    if 1.0*int(line[9]) / int(line[1]) > cov and line[0]!=line[5] and [line[0], line[5]] not in reads:
       reads.append([line[0], line[5]])
f.close()

G = nx.Graph()
G.add_nodes_from(sum(reads, []))
q = [[(s[i],s[i+1]) for i in range(len(s)-1)] for s in reads]
f_write = open(path + 'cluster.fa','w')
for i in q:
    G.add_edges_from(i)
reads1 = set()
for info in nx.connected_components(G):
    length = 0
    for i in info:
        reads1.add(i)
        if len(sequence[i]) > length:
           rep = i
           length = len(sequence[i])
    f_write.write('>' + rep + '\n')
    f_write.write(sequence[rep] + '\n')
for seq in sequence:
    if seq not in reads1:
       f_write.write('>' + seq + '\n')
       f_write.write(sequence[seq] + '\n')
f_write.close()


