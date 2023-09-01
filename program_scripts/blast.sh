#!/bin/bash

module load blast-plus

blastn \
	-query ${1} \
	-db reference/blast_reference_trinity \
	-out blast/${2}.txt \
	-outfmt "6 qseqid sseqid qcovhsp pident length mismatch qstart qend sstart send evalue"

sed -i "s/$/\t${2}/" blast/${2}.txt
