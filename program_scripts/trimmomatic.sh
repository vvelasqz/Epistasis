#!/bin/bash

module load trimmomatic

trimmomatic SE \
   -phred33 \
   samples/Sample_${1}/*.fastq.gz \
   qc_reads/${1}.qcreads.fastq.gz \
   LEADING:15 \
   TRAILING:15 \
   SLIDINGWINDOW:4:20 \
   MINLEN:30 \
   ILLUMINACLIP:adapters.fa:2:30:10

