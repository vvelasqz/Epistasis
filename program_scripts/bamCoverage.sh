#!/bin/bash

module load deeptools

sample=`basename ${1%.*}`

bamCoverage \
	-p 8 \
	-b ${1} \
	-o bw/${sample}.bw
