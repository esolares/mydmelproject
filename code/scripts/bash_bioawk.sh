#!/bin/bash
DATADIR="../../data/raw"
INPUT="dmel-all-chromosome-r6.13.fasta"
PROCDIR="../../data/processed"
OUTDIR="../../output/reports"
FIGDIR="../../output/figures"
GTFFILE="dmel-all-r6.13.gtf"

for i in $(ls $DATADIR/*.fasta)
	do echo -e "Length\tContig\tAssembly" > ${PROCDIR}/$(basename ${i} .fasta).lengths; bioawk -c fastx '{ print length($seq), $name, "dmel_chrom" }' ${i} | sort -rnk 1,1 >> ${PROCDIR}/$(basename ${i} .fasta).lengths; plotCDF2 ${PROCDIR}/$(basename ${i} .fasta).lengths ${FIGDIR}/$(basename ${i} .fasta).png
done

for i in $(ls $PROCDIR/*.fasta)
        do echo -e "Length\tContig\tAssembly" > ${PROCDIR}/$(basename ${i} .fasta).lengths; bioawk -c fastx '{ print length($seq), $name, "dmel_chrom" }' ${i} | sort -rnk 1,1 >> ${PROCDIR}/$(basename ${i} .fasta).lengths; plotCDF2 ${PROCDIR}/$(basename ${i} .fasta).lengths ${FIGDIR}/$(basename ${i} .fasta).png
done

###get summary of nucleotide make up
grep -v '^>' $DATADIR/$INPUT | awk 'BEGIN{a=0; c=0; g=0; t=0; n=0;} {IGNORECASE=1; a+=gsub("A",""); c+=gsub("C",""); g+=gsub("G",""); t+=gsub("T",""); n+=gsub("N","")} END{print "Total # of Nucleotides: " a+c+g+t " bases", "Total # of Ns: " n, "Total Bases: " a+c+g+t+n " bases"}' > ${OUTDIR}/$(basename $INPUT .fasta)_summary.txt

###get gc content per contig
bioawk -cfastx '{print $name, gc($seq)}' ${DATADIR}/${INPUT} > ${PROCDIR}_gc_content.txt

###
#Filter reads by size
faFilter -maxSize=100000 $DATADIR/$INPUT $PROCDIR/$(basename $INPUT .fasta)_lt100kb_filtered.fasta
faFilter -minSize=100000 $DATADIR/$INPUT $PROCDIR/$(basename $INPUT .fasta)_gt100kb_filtered.fasta

###get summary of nucleotide make up for filtered reads
#<100kb
grep -v '^>' $PROCDIR/$(basename $INPUT .fasta)_lt100kb_filtered.fasta | awk 'BEGIN{a=0; c=0; g=0; t=0; n=0;} {IGNORECASE=1; a+=gsub("A",""); c+=gsub("C",""); g+=gsub("G",""); t+=gsub("T",""); n+=gsub("N","")} END{print "Total # of Nucleotides: " a+c+g+t " bases", "Total # of Ns: " n, "Total Bases: " a+c+g+t+n " bases"}' > ${OUTDIR}/$(basename $INPUT .fasta)_lt100kb_filtered_summary.txt
#>100kb
grep -v '^>' $PROCDIR/$(basename $INPUT .fasta)_gt100kb_filtered.fasta | awk 'BEGIN{a=0; c=0; g=0; t=0; n=0;} {IGNORECASE=1; a+=gsub("A",""); c+=gsub("C",""); g+=gsub("G",""); t+=gsub("T",""); n+=gsub("N","")} END{print "Total # of Nucleotides: " a+c+g+t " bases", "Total # of Ns: " n, "Total Bases: " a+c+g+t+n " bases"}' > ${OUTDIR}/$(basename $INPUT .fasta)_gt100kb_filtered_summary.txt

###get gc content per contig
bioawk -cfastx '{print $name, gc($seq)}' $PROCDIR/$(basename $INPUT .fasta)_lt100kb_filtered.fasta > ${PROCDIR}_lt100kb_filtered_gc_content.txt
bioawk -cfastx '{print $name, gc($seq)}' $PROCDIR/$(basename $INPUT .fasta)_gt100kb_filtered.fasta > ${PROCDIR}_gt100kb_filtered_gc_content.txt

###count 
cat $DATADIR/$GTFFILE | cut -f 1 | sort -k1,1 | uniq -c | sort -nr >> $OUTDIR/$(basename $GTFFILE .gtf)_count_per_xsome.txt

