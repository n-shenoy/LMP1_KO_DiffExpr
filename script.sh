#!/bin/bash
#Author: Navami Shenoy
#pipeline for quantifying gene expression for control and LMP1 knockout LCL samples

SECONDS=0 

#change working directory
cd /mnt/f/bulkRNASeq/LMP1

#----------------DOWNLOADS----------------#

#make directories for replicates
mkdir Control_Rep{1..3} LMP1_KO_Rep{1..3}
mkdir counts #for saving counts

#Control
#BAM file for replicate 1
wget -O Control_Rep1/alignments.bam http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeCaltechRnaSeq/wgEncodeCaltechRnaSeqGm12892R2x75Il200AlignsRep1V2.bam
#Replicate 2 FASTQ
wget -O Control_Rep2/reads1.fastq.gz http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeCaltechRnaSeq/wgEncodeCaltechRnaSeqGm12892R2x75Il200FastqRd1Rep2V2.fastq.gz
wget -O Control_Rep2/reads2.fastq.gz http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeCaltechRnaSeq/wgEncodeCaltechRnaSeqGm12892R2x75Il200FastqRd2Rep2V2.fastq.gz
#Replicate 3 FASTQ
wget -O Control_Rep3/reads1.fastq.gz http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeCaltechRnaSeq/wgEncodeCaltechRnaSeqGm12892R2x75Il200FastqRd1Rep3V2.fastq.gz
wget -O Control_Rep3/reads2.fastq.gz http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeCaltechRnaSeq/wgEncodeCaltechRnaSeqGm12892R2x75Il200FastqRd2Rep3V2.fastq.gz

#LMP1 Knockout
#Replicate 1
wget -O LMP1_KO_Rep1/reads1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR239/010/SRR23957810/SRR23957810_1.fastq.gz
wget -O LMP1_KO_Rep1/reads2.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR239/010/SRR23957810/SRR23957810_2.fastq.gz
#Repilcate 2
wget -O LMP1_KO_Rep2/reads1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR239/011/SRR23957811/SRR23957811_1.fastq.gz
wget -O LMP1_KO_Rep2/reads2.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR239/011/SRR23957811/SRR23957811_2.fastq.gz
#Replicate 3
wget -O LMP1_KO_Rep2/reads1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR239/012/SRR23957812/SRR23957812_1.fastq.gz
wget -O LMP1_KO_Rep2/reads2.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR239/012/SRR23957812/SRR23957812_2.fastq.gz


#Genome indices
mkdir HISAT2/
wget https://genome-idx.s3.amazonaws.com/hisat/grch38_genome.tar.gz
sudo tar -xvzf grch38_genome.tar.gz -C HISAT2/ #unzip

#Genome annotation file
wget https://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/Homo_sapiens.GRCh38.110.gtf.gz
sudo gunzip -k Homo_sapiens.GRCh38.110.gtf.gz #unzip
mv Homo_sapiens.GRCh38.110.gtf GRCh38.gtf #rename


#------------QUALITY CONTROL-----------#

#CONTROL

#generate quality report using fastqc for paired-end reads
find -name '*.fastq.gz' | xargs fastqc -t 8
echo "fastqc finished running!"


for i in 2 3;
    do
        #lots of reads entirely consisting of Ns in both replicates
        #filter out paired-end reads with more than 2 Ns using cutadapt
        cutadapt --max-n 2 --cores=3 -o Control_Rep${i}/nfiltered.reads1.fastq.gz -p Control_Rep${i}/nfiltered.reads2.fastq.gz\
            Control_Rep${i}/reads1.fastq.gz Control_Rep${i}/reads2.fastq.gz

        #For replicates 2 and 3, bases towards the 3' ends have poor quality (<25)
        #trim reads from the 3' end while maintaining a read length of 50 and strictless of 0.8:
        java -jar /usr/share/java/trimmomatic-0.39.jar PE -threads 8 -phred64 Control_Rep${i}/nfiltered.reads1.fastq.gz Control_Rep${i}/nfiltered.reads2.fastq.gz\
            Control_Rep${i}/paired.reads1.fastq.gz Control_Rep${i}/unpaired.reads1.fastq.gz\
            Control_Rep${i}/paired.reads2.fastq.gz Control_Rep${i}/unpaired.reads2.fastq.gz MAXINFO:50:0.8 MINLEN:50
        #83.16% and 79% of Rep2 and Rep3 paired-end reads survived trimming

        #finally, filter out reads with mean base quality below 20
        java -jar /usr/share/java/trimmomatic-0.39.jar PE -threads 8 -phred64 Control_Rep${i}/paired.reads1.fastq.gz Control_Rep${i}/paired.reads2.fastq.gz\
            Control_Rep${i}/paired.reads1.filt.fastq.gz Control_Rep${i}/unpaired.reads1.filt.fastq.gz\
            Control_Rep${i}/paired.reads2.filt.fastq.gz Control_Rep${i}/unpaired.reads2.filt.fastq.gz AVGQUAL:20
        #99% pairs survived filtering
    done


#LMP1 KNOCKOUT

#these reads already had great base quality but have a lot of adapter content
#remove Illumina universal adapter from the reads using cutadapt
for i in 1 2 3;
    do
        cutadapt -a AGATCGGAAGAG -A AGATCGGAAGAG --cores=4 -o LMP1_KO_Rep${i}/paired.reads1.filt.fastq.gz -p LMP1_KO_Rep${i}/paired.reads2.filt.fastq.gz\
            LMP1_KO_Rep${i}/reads1.fastq.gz LMP1_KO_Rep${i}/reads2.fastq.gz
    done

#run fastqc on filtered
find -name '*.filt.fastq.gz' | xargs fastqc -t 8
echo "fastqc finished running again!"


#---------------MAPPING & QUANTITATION----------------#

#CONTROL
for i in 2 3; 
    do 
        #map reads to reference genome using HISAT2
        hisat2 -q --phred64 -p 8 -x HISAT2/grch38/genome -1 Control_Rep${i}/paired.reads1.filt.fastq.gz\
        -2 Control_Rep${i}paired.reads1.filt.fastq.gz | samtools view -bS -o Control_Rep${i}/alignments.bam
        echo "hisat2 finished aligning Control_Rep${i}!"

        #alignment summary
        samtools flagstat -@ 8 Control_Rep${i}/alignments.bam

        #sort by chromosome coordinates - for genome browsers, analysis tools
        samtools sort --threads 8 Control_Rep${i}/alignments.bam -o Control_Rep${i}/alignments.sorted.bam 

        #create a bam index for IGV
        samtools index -@ 8 Control_Rep${i}/alignments.sorted.bam

        #sort by read names - for expression quantitation tools
        samtools sort --threads 8 -n Control_Rep${i}/alignments.bam -o Control_Rep${i}/alignments.namesorted.bam
        echo "samtools finished sorting Control_Rep${i}!"

        #quantitate gene counts using featureCounts
        featureCounts -p -T 8 -s 0 -Q 10 -a GRCh38.gtf -o Control_Rep${i}/Control_Rep${i}_featureCounts.txt Control_Rep${i}/alignments.namesorted.bam
        #Note: for paired-end reads, we should count read pairs (fragments) rather than reads 
        #because counting fragments will give us more accurate counts, so we're using the option -p here
        echo "featureCounts finished quantitating Control_Rep${i}!"

        #extract gene IDs and counts from the featureCounts output
        cut -f1,7 Control_Rep${i}/Control_Rep${i}_featureCounts.txt > counts/Control_Rep${i}_Counts.txt 
        echo "Control_Rep${i} counts file is ready!"
    done 

#KNOCKOUT
for i in 1 2 3; 
do 
    #mapping
    hisat2 -q --phred64 -p 8 -x HISAT2/grch38/genome -1 LMP1_KO_Rep${i}/paired.reads1.filt.fastq.gz\
     -2 LMP1_KO_Rep${i}paired.reads1.filt.fastq.gz | samtools view -bS -o LMP1_KO_Rep${i}/alignments.bam
    echo "hisat2 finished aligning LMP1_KO_Rep${i}!" 

    #sorting
    samtools flagstat -@ 8 LMP1_KO_Rep${i}/alignments.bam
    samtools sort --threads 8 LMP1_KO_Rep${i}/alignments.bam -o LMP1_KO_Rep${i}/alignments.sorted.bam 
    samtools index -@ 8 LMP1_KO_Rep${i}/alignments.sorted.bam
    samtools sort --threads 8 -n LMP1_KO_Rep${i}/alignments.bam -o LMP1_KO_Rep${i}/alignments.namesorted.bam
    echo "samtools finished sorting LMP1_KO_Rep${i}!"

    #quantitation
    featureCounts -p -T 8 -s 0 -Q 10 -a GRCh38.gtf -o LMP1_KO_Rep${i}/LMP1_KO_Rep${i}_featureCounts.txt LMP1_KO_Rep${i}/alignments.namesorted.bam
    echo "featureCounts finished quantitating LMP1_KO_Rep${i}!"

    cut -f1,7 LMP1_KO_Rep${i}/LMP1_KO_Rep${i}_featureCounts.txt > counts/LMP1_KO_Rep${i}_Counts.txt 
    echo "LMP1_KO_Rep${i} counts file is ready!"
done 


#create a counts matrix from all replicates  
#note: one of the control samples was dropped so at the end we had n=5
cd counts
paste * | awk 'BEGIN {OFS="\t"; FS="\t"}; {j=$1; for (i=2;i<=NF;i+=2) {j=j FS $i} print j}' > counts_matrix_5_samples.csv


#time elapsed
hours=$SECONDS/3600
mins=($SECONDS/60)%60
seconds=$SECONDS%60
echo "Time elapsed: $(($hours)) hours, $(($mins)) minutes and $(($seconds)) seconds"