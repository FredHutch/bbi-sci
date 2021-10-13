#!/bin/bash
set -euo pipefail

cat ${logfile} > remove_dups.log
printf "** Start processes 'remove duplicates, assign_genes, merge_assignment' at: \$(date)\n\n" >> remove_dups.log
printf "    Process versions:
    \$(bedtools --version)
    \$(samtools --version | tr '\n' ' ')
    \$(bamtools --version | grep bamtools)
    \$(python --version)\n\n" >> remove_dups.log

printf "    Process stats:
    remove_dups starting reads: \$(samtools view -c $merged_bam)" >> remove_dups.log


mkdir split_bams
bamtools split -in $merged_bam -reference -stub split_bams/split
cd split_bams
if [[ \$(ls | grep "_[0-9A-Za-z\\.]\\{3,\\}.bam\$") ]]; then
    ls | grep "_[0-9A-Za-z\\.]\\{3,\\}.bam\$" | samtools merge split.REFnonstand.bam -b -
    ls | grep "_[0-9A-Za-z\\.]\\{3,\\}.bam\$" | xargs -d"\\n" rm
    mv split.REFnonstand.bam split.REF_nonstand.bam
fi