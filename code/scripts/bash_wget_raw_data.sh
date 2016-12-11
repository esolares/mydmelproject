#!/bin/bash
DATADIR="../../data/raw"

wget -P $DATADIR ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r6.13_FB2016_05/fasta/dmel-all-chromosome-r6.13.fasta.gz
wget -P $DATADIR ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r6.13_FB2016_05/fasta/dmel-all-transcript-r6.13.fasta.gz
wget -P $DATADIR ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r6.13_FB2016_05/gtf/dmel-all-r6.13.gtf.gz
wget -P $DATADIR ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r6.13_FB2016_05/fasta/dmel-all-CDS-r6.13.fasta.gz

for i in $(ls ${DATADIR}/*.gz)
	do gunzip $i
done
