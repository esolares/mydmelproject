#!/usr/bin/env python
from collections import defaultdict
mygenedict = defaultdict(set)
for line in open('../../data/processed/dmel-all-r6.13_exonlengths.txt','r'):
    myline = line.replace("\n","").replace('"',"").replace(";","").split(" ")
    mygenedict[myline[3]].add(myline[7])

fout = open('../../output/reports/dmel-all-r6.13_transcripts_per_gene.txt','w')

for key in mygenedict:
    fout.write(key + '\t' + str(len(mygenedict[key])) + '\n')

fout.close()
quit()
