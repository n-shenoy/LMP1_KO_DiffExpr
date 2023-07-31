#!/bin/bash

SECONDS=0 #in seconds

#change working directory
#cd /mnt/f/bulkRNASeq/GM12878

#Download FASTQ files for GM12878 replicate 2:
#wget -O Gm12892Rd1Rep2.fastq.gz http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeCaltechRnaSeq/wgEncodeCaltechRnaSeqGm12892R2x75Il200FastqRd1Rep2V2.fastq.gz
#wget -O Gm12892Rd2Rep2.fastq.gz http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeCaltechRnaSeq/wgEncodeCaltechRnaSeqGm12892R2x75Il200FastqRd2Rep2V2.fastq.gz

#------------QUALITY CONTROL-----------#

#generate quality report using fastqc for both files
#fastqc Gm12892Rd1Rep2.fastq.gz
#fastqc Gm12892Rd2Rep2.fastq.gz
echo "fastqc finished running!"

#lots of reads entirely consisting of Ns in both reads 
#filtering out paired-end reads with more than 2 Ns using cutadapt
#cutadapt --max-n 2 -o nfiltered.Rd1Rep2.fastq.gz -p nfiltered.Rd2Rep2.fastq.gz Gm12892Rd1Rep2.fastq.gz Gm12892Rd2Rep2.fastq.gz
#685,175 (2.0%) reads had more than 2 instances of N bases
#About 98.0% (33,546,906) reads survived the filtering process

#In the pair file, bases at positions 50 and above have poor quality (<25)
#Plus, starting positions of the paired-end reads have poor sequence content per base
#So, we're gonna trim the reads with poor starting and trailing base qualities 
#and filter out reads that are too short

#trimming reads from the 3' end while maintaining a read length of 50 and strictless of 0.8:
#java -jar /usr/share/java/trimmomatic-0.39.jar PE -threads 4 -phred64 nfiltered.Rd1Rep2.fastq.gz nfiltered.Rd2Rep2.fastq.gz paired1.Rd1Rep2.fastq.gz unpaired1.Rd1Rep2.fastq.gz paired2.Rd2Rep2.fastq.gz unpaired2.Rd2Rep2.fastq.gz MAXINFO:50:0.8 MINLEN:50
#83.16% of paired-end reads survived trimming

#finally, filtering out reads with mean base quality below 20
#java -jar /usr/share/java/trimmomatic-0.39.jar PE -threads 4 -phred64 paired1.Rd1Rep2.fastq.gz paired2.Rd2Rep2.fastq.gz filt.paired1.Rd1Rep2.fastq.gz filt.unpaired1.Rd1Rep2.fastq.gz filt.paired2.Rd2Rep2.fastq.gz filt.unpaired2.Rd2Rep2.fastq.gz AVGQUAL:20
#99.93% pairs survived filtering!

#running fastqc on processed reads:
#fastqc filt.paired1.Rd1Rep2.fastq.gz
#fastqc filt.paired2.Rd2Rep2.fastq.gz
echo "fastqc finished running again!"
#-----------MAPPING----------------#

#HISAT2
#mkdir HISAT2
#download genome indices
#wget https://genome-idx.s3.amazonaws.com/hisat/grch38_genome.tar.gz
#sudo tar -xvzf grch38_genome.tar.gz -C HISAT2/
#hisat2 -q --phred64 -p 6 -x HISAT2/grch38/genome -1 filt.paired1.Rd1Rep2.fastq.gz -2 filt.paired2.Rd2Rep2.fastq.gz | samtools view -bS -o alignments.bam
#hisat2 took about 2.15 hours to run on one thread
echo "hisat2 finished running!"

#alignment summary
#samtools flagstat -@ 6 alignments.bam -O alignment_summary.txt
#Overall alignment rate was 96.14% and about 76% of the reads aligned concordantly exactly one time

#view alignments
#samtools view -h alignments.bam | less
#sort by chromosome coordinates - for genome browsers, analysis tools
#samtools sort --threads 6 alignments.bam -o alignments.sorted.bam 
#sort by read names - for expression quantitation tools
#samtools sort --threads 8 -n alignments.bam -o alignments.namesorted.bam 
#create a bam index for IGV
#samtools index -@ 6 alignments.sorted.bam
echo "samtools finished sorting!"

#------------QUANTITATION--------------#

#download and unzip genome annotation file
#wget https://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/Homo_sapiens.GRCh38.110.gtf.gz

#sudo gunzip -k Homo_sapiens.GRCh38.110.gtf.gz
#mv Homo_sapiens.GRCh38.110.gtf GRCh38.gtf

#quantitating gene expression using featureCounts
#mkdir counts 
#featureCounts -p -T 6 -s 0 -Q 10 -a GRCh38.gtf -o counts/GM12878_Rep2_featureCounts.txt alignments.sorted.bam
#Note: for paired-end reads, we should count read pairs (fragments) rather than reads 
#because counting fragments will give us more accurate counts, so we're using the option -p here
#featureCounts took almost 4 hrs to run :(
echo "featureCounts finished running!"

#extract gene IDs and counts from the featureCounts output
#cut -f1,7 counts/GM12878_Rep2_featureCounts.txt > counts/GM12878_Rep2_Counts.txt 
echo "counts file is ready!"

hours=$SECONDS/3600
mins=($SECONDS/60)%60
seconds=$SECONDS%60
echo "Time elapsed: $(($hours)) hours, $(($mins)) minutes and $(($seconds)) seconds"
