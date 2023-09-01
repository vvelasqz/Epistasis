#!/bin/bash

module load salmon

sample=$(basename "$2")
sample="${sample%%.*}"

salmon quant \
	-t ${1} \
	--libType A \
	-a ${2} \
	-o quantifications/${sample} \
	-p 8
