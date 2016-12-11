#!/bin/bash
DATADIR="../../data/raw"
INPUT="dmel-all-chromosome-r6.13.fasta"
PROCDIR="../../data/processed"
OUTDIR="../../output/reports"
FIGDIR="../../output/figures"
GTFFILE="dmel-all-r6.13.gtf"

#./bash_wget_raw_data.sh

for i in $(ls $DATADIR/*.fasta)
	do echo -e "Length\tContig\tAssembly" > ${PROCDIR}/$(basename ${i} .fasta).lengths; bioawk -c fastx '{ print length($seq), $name, "dmel_chrom" }' ${i} | sort -rnk 1,1 >> ${PROCDIR}/$(basename ${i} .fasta).lengths; plotCDF2 ${PROCDIR}/$(basename ${i} .fasta).lengths ${FIGDIR}/$(basename ${i} .fasta).png
done

for i in $(ls $PROCDIR/*.fasta)
        do echo -e "Length\tContig\tAssembly" > ${PROCDIR}/$(basename ${i} .fasta).lengths; bioawk -c fastx '{ print length($seq), $name, "dmel_chrom" }' ${i} | sort -rnk 1,1 >> ${PROCDIR}/$(basename ${i} .fasta).lengths; plotCDF2 ${PROCDIR}/$(basename ${i} .fasta).lengths ${FIGDIR}/$(basename ${i} .fasta).png
done

###get summary of nucleotide make up
grep -v '^>' $DATADIR/$INPUT | awk 'BEGIN{a=0; c=0; g=0; t=0; n=0;} {IGNORECASE=1; a+=gsub("A",""); c+=gsub("C",""); g+=gsub("G",""); t+=gsub("T",""); n+=gsub("N","")} END{print "Total # of Nucleotides: " a+c+g+t " bases", "Total # of Ns: " n, "Total Bases: " a+c+g+t+n " bases"}' > ${OUTDIR}/$(basename $INPUT .fasta)_summary.txt

###get gc content per contig
bioawk -cfastx '{print $name, gc($seq)}' ${DATADIR}/${INPUT} > ${PROCDIR}/$(basename $INPUT .fasta)_gc_content.txt

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
#bioawk -cfastx '{print $name, gc($seq)}' $PROCDIR/$(basename $INPUT .fasta)_lt100kb_filtered.fasta > ${PROCDIR}/$(basename $INPUT .fasta)_lt100kb_filtered_gc_content.txt
#bioawk -cfastx '{print $name, gc($seq)}' $PROCDIR/$(basename $INPUT .fasta)_gt100kb_filtered.fasta > ${PROCDIR}/$(basename $INPUT .fasta)_gt100kb_filtered_gc_content.txt

###count 
cat $DATADIR/$GTFFILE | bioawk -cgff '$3=="gene" {print $1 "\t" $end-$start "\t" $9}' > ${PROCDIR}/$(basename $GTFFILE .gtf)_genelengths.txt
cat $DATADIR/$GTFFILE | bioawk -cgff '$3=="exon" {print $1 "\t" $end-$start "\t" $9}' > ${PROCDIR}/$(basename $GTFFILE .gtf)_exonlengths.txt
cat $DATADIR/$GTFFILE | cut -f 3 | sort -k1,1 | uniq -c | sort -nr > ${PROCDIR}/$(basename $GTFFILE .gtf)_sorted_feature_counts.txt
cat ${PROCDIR}/$(basename $GTFFILE .gtf)_genelengths.txt | cut -f 1 | sort -k1,1 | uniq -c | sort -nr > $OUTDIR/$(basename $GTFFILE .gtf)_gene_count_per_xsome.txt
#cat ${PROCDIR}/$(basename $GTFFILE .gtf)_exonlengths.txt | cut -f 1 | sort -k1,1 | uniq -c | sort -nr > $OUTDIR/$(basename $GTFFILE .gtf)_exon_count_per_xsome.txt

./bash_python_tr_per_gene.py
./bash_rscript_plots.r

cat $OUTDIR/$(basename $GTFFILE .gtf)_transcripts_per_gene.txt | sort -rnk 2,2 > $OUTDIR/$(basename $GTFFILE .gtf)_transcripts_per_gene_sorted.txt
