#!/bin/bash

module load samtools

sample=`basename ${1%.*}`

samtools view \
	-@ 8 \
	-bS ${1} \
	| \
	samtools sort \
	- \
	-o sorted_alignments/${sample}.sorted.bam
	
samtools index \
	-@ 8 \
	-o sorted_alignments/${sample}.sorted.bam.bai \
	sorted_alignments/${sample}.sorted.bam
