# Differential Expression Analysis of LMP1 Knockout Lymphoblastoid Cell Lines
Complete sequencing pipeline for differential gene expression analysis of bulk RNA from human latent membrane protein 1 (LMP1) knockout (KO) lymphoblastoid cell lines (LCL). 

### Control RNA
The control transcriptome sequence was obtained from GM12878, a well-established LCL produced from the blood of a female donor by Epstein-Barr virus (EBV) transformation as a part of the ENCODE project. The control consists of raw paired-end reads from GM12878 replicates 2 and 3, and can be downloaded from the command line:
```Bash
#Replicate 2
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeCaltechRnaSeq/wgEncodeCaltechRnaSeqGm12892R2x75Il200FastqRd1Rep2V2.fastq.gz
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeCaltechRnaSeq/wgEncodeCaltechRnaSeqGm12892R2x75Il200FastqRd2Rep2V2.fastq.gz
#Replicate 3
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeCaltechRnaSeq/wgEncodeCaltechRnaSeqGm12892R2x75Il200FastqRd1Rep3V2.fastq.gz
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeCaltechRnaSeq/wgEncodeCaltechRnaSeqGm12892R2x75Il200FastqRd2Rep3V2.fastq.gz
```
For replicate 1, the BAM file with aligned reads was used due to large FASTQ files and hardware limitations:
```Bash
http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeCaltechRnaSeq/wgEncodeCaltechRnaSeqGm12892R2x75Il200AlignsRep1V2.bam
```

### LMP1 KO RNA
The RNA sequence data for LMP1 knockout LCLs came from Mitra et al 2023, (GEO Accession ID GSE228167) generated from the GM12878 cell lines and consists of three replicates. The data can be obtained from:
```Bash
#Replicate 1
#wget -O reads1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR239/010/SRR23957810/SRR23957810_1.fastq.gz
#wget -O reads2.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR239/010/SRR23957810/SRR23957810_2.fastq.gz
#Repilcate 2
#wget -O reads1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR239/011/SRR23957811/SRR23957811_1.fastq.gz
#wget -O reads2.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR239/011/SRR23957811/SRR23957811_2.fastq.gz
#Replicate 3
#wget -O reads1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR239/012/SRR23957812/SRR23957812_1.fastq.gz
#wget -O reads2.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR239/012/SRR23957812/SRR23957812_2.fastq.gz
```


`fastqc` contains quality reports for raw and processed reads.

The bash script for the complete RNA-Seq pipeline can be found in `script.sh`.

References:
1. [Mitra et al. (2023).](https://www.biorxiv.org/content/10.1101/2023.04.10.536234v1.full) Characterization of Target Gene Regulation by the Two Epstein-Barr Virus Oncogene LMP1 Domains Essential for B-cell Transformation.
2. [RNA-seq Data Analysis: A Practical Approach](https://www.amazon.com/RNA-seq-Data-Analysis-Mathematical-Computational/dp/1466595000/ref=sr_1_1?crid=2WJG6RRFNQTWV&keywords=Rna-Seq+Data+Analysis%3A+A+Practical+Approach&qid=1689673621&s=books&sprefix=rna-seq+data+analysis+a+practical+approach%2Cstripbooks%2C316&sr=1-1&ufe=app_do%3Aamzn1.fos.006c50ae-5d4c-4777-9bc0-4513d670b6bc) 
