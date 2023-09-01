#!/bin/bash

module load bwa

ref=$(basename "$1")
ref="${ref%.*}"

sample=$(basename "$2")
sample="${sample%%.*}"

bwa mem ${1} ${2} > alignments/${sample}_${ref}.sam
