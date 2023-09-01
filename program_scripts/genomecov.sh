#!/bin/bash

module load bedtools2

sample=`basename ${1%.*}`

bedtools \
	genomecov \
	-ibam ${1} \
	-bga > genomecov/${sample}_genome_coverage.txt
