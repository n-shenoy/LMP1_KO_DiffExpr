#!/bin/bash

SECONDS=0 #in seconds

#change working directory
#cd /mnt/f/bulkRNASeq/GSE228167

#Download FASTQ files for GSE228167 
#replicate 3
#wget -O reads1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR239/012/SRR23957812/SRR23957812_1.fastq.gz
#wget -O reads2.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR239/012/SRR23957812/SRR23957812_2.fastq.gz

#------------QUALITY CONTROL-----------#

#generate quality report using fastqc for both files
#fastqc -t 8 reads1.fastq.gz
#fastqc -t 8 reads2.fastq.gz
#echo "fastqc finished running!"

#removing Illumina universal adapter from the reads
#cutadapt -a AGATCGGAAGAG -A AGATCGGAAGAG --cores=4 -o adapter.filtered.reads1.fastq.gz -p adapter.filtered.reads2.fastq.gz\
#	reads1.fastq.gz reads2.fastq.gz

#running fastqc on processed reads:
#fastqc -t 8 adapter.filtered.reads1.fastq.gz
#fastqc -t 8 adapter.filtered.reads2.fastq.gz
#echo "fastqc finished running again!"

#-----------MAPPING----------------#

#HISAT2
#mkdir HISAT2
#download genome indices
#wget https://genome-idx.s3.amazonaws.com/hisat/grch38_genome.tar.gz
#sudo tar -xvzf grch38_genome.tar.gz -C HISAT2/
#hisat2 --phred33 -p 8 -x HISAT2/grch38/genome -1 adapter.filtered.reads1.fastq.gz -2 adapter.filtered.reads2.fastq.gz | samtools view -bS -o alignments.bam
#echo "hisat2 finished running!"

#alignment summary
#samtools flagstat -@ 8 alignments.bam -O alignment_summary.tsv
#Overall alignment rate was 95.75% and 88.76% of the reads aligned concordantly exactly one time

#view alignments
#samtools view -h alignments.bam | less
#sort by chromosome coordinates - for genome browsers, analysis tools
#samtools sort --threads 8 alignments.bam -o alignments.sorted.bam 
#create a bam index for IGV
#samtools index -@ 8 alignments.sorted.bam
#sort by read names - for expression quantitation tools
#samtools sort --threads 8 -n alignments.bam -o alignments.namesorted.bam
#echo "samtools finished sorting!"

#------------QUANTITATION--------------#

#download and unzip genome annotation file
#wget https://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/Homo_sapiens.GRCh38.110.gtf.gz

#sudo gunzip -k Homo_sapiens.GRCh38.110.gtf.gz
#mv Homo_sapiens.GRCh38.110.gtf GRCh38.gtf

#quantitating gene expression using featureCounts
#mkdir counts 
featureCounts -p -T 8 -s 0 -Q 10 -a GRCh38.gtf -o counts/GSE228167_Rep3_featureCounts.txt alignments.namesorted.bam
#Note: for paired-end reads, we should count read pairs (fragments) rather than reads 
#because counting fragments will give us more accurate counts, so we're using the option -p here
echo "featureCounts finished running!"

#extract gene IDs and counts from the featureCounts output
cut -f1,7 counts/GSE228167_Rep3_featureCounts.txt > counts/GSE228167_Rep3_Counts.txt 
echo "counts file is ready!"

hours=$SECONDS/3600
mins=($SECONDS/60)%60
seconds=$SECONDS%60
echo "Time elapsed: $(($hours)) hours, $(($mins)) minutes and $(($seconds)) seconds"
