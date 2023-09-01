#!/bin/bash

module load trinity

sample=$(basename "${1}")
sample="${sample%%.*}"

Trinity \
	--genome_guided_bam ${1} \
	--genome_guided_max_intron 1000 \
	--max_memory 16G \
	--CPU 8 \
	--output trinity/${sample}_trinity_gg \
	--no_salmon
