#!/bin/bash

SECONDS=0 #in seconds

#change working directory
#cd /mnt/f/bulkRNASeq/GM12878

#Download BAM file for GM12878 replicate 1: 
#wget -O alignments.bam http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeCaltechRnaSeq/wgEncodeCaltechRnaSeqGm12892R2x75Il200AlignsRep1V2.bam

#-----------MAPPING----------------#

#view alignments
#samtools view -h alignments.bam | less
#sort by chromosome coordinates - for genome browsers, analysis tools
#samtools sort --threads 6 alignments.bam -o alignments.sorted.bam 
#sort by read names - for expression quantitation tools
#samtools sort --threads 8 -n alignments.bam -o alignments.namesorted.bam 
#create a bam index for IGV
#samtools index -@ 6 alignments.sorted.bam
#echo "samtools finished sorting!"

#------------QUANTITATION--------------#

#download and unzip genome annotation file
#wget https://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/Homo_sapiens.GRCh38.110.gtf.gz

#sudo gunzip -k Homo_sapiens.GRCh38.110.gtf.gz
#mv Homo_sapiens.GRCh38.110.gtf GRCh38.gtf

#quantitating gene expression using featureCounts
#mkdir counts 
#featureCounts -p -T 6 -s 0 -Q 10 -a GRCh38.gtf -o counts/GM12878_Rep1_featureCounts.txt alignments.namesorted.bam
#Note: for paired-end reads, we should count read pairs (fragments) rather than reads 
#because counting fragments will give us more accurate counts, so we're using the option -p here
#featureCounts took almost 4 hrs to run :(
#echo "featureCounts finished running!"

#extract gene IDs and counts from the featureCounts output
#cut -f1,7 counts/GM12878_Rep1_featureCounts.txt > counts/GM12878_Rep1_Counts.txt 
#echo "counts file is ready!"

#create a matrix of counts from all replicates  
cd counts
paste * | awk 'BEGIN {OFS="\t"; FS="\t"}; {j=$1; for (i=2;i<=NF;i+=2) {j=j FS $i} print j}' > counts_matrix.csv

hours=$SECONDS/3600
mins=($SECONDS/60)%60
seconds=$SECONDS%60
echo "Time elapsed: $(($hours)) hours, $(($mins)) minutes and $(($seconds)) seconds"
