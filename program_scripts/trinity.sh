#!/bin/bash

module load trinity

sample=`basename ${1%.*}`

ends_with() {
    case "$1" in
        *"$2") return 0 ;; # Return 0 if the suffix is found
        *) return 1 ;;      # Return 1 otherwise
    esac
}

if ends_with ${1} "fa" || ends_with ${1} "fasta" || ends_with ${1} "fna"; then
    file_type="fa"
else
    file_type="fq"
fi

Trinity \
	--seqType fq \
	--single ${1} \
	--max_memory 128G \
	--CPU 8 \
	--output trinity/${sample}_trinity
