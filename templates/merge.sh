#!/bin/bash
set -euo pipefail

printf "** Start process 'merge_bams' at: \$(date)\n\n"
printf "    Process versions:
    \$(samtools --version | tr '\n' ' ')\n\n"
printf "    Process command:
    samtools merge ${key}.bam $sorted_bam\n\n"


samtools merge -@ $cores_merge ${key}.bam $sorted_bam


printf "${key}\t\$(samtools view -c ${key}.bam)" > ${key}.read_count.txt

printf "** End process 'merge_bams' at: \$(date)\n\n"