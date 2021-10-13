#!/bin/bash
set -euo pipefail

printf "** Start process 'sort_and_filter' for $aligned_bam at: \$(date)\n\n" > ${name}_sf.log
printf "    Process versions:
    \$(samtools --version | tr '\n' ' ')\n\n" >> ${name}_sf.log

samtools view -bh -q 30 -F 4 "$aligned_bam" \
    | samtools sort -@ $cores_sf - \
    > "${name}.bam"

printf "    Process stats:
    sort_and_filter starting reads: \$(samtools view -c $aligned_bam)
    sort_and_filter ending reads  : \$(samtools view -c ${name}.bam)\n\n" >> ${name}_sf.log
printf "** End process 'sort_and_filter' at: \$(date)\n\n" >> ${name}_sf.log
