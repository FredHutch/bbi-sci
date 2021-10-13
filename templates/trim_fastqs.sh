#!/usr/bin/env bash

set -euo pipefail

printf "** Start process 'trim_fastqs' for ${input_fastq} "
printf "Process versions: "
printf "  python      \$(python --version &>> piece.log) "
printf "  trim_galore \$(trim_galore -v | grep version | awk '{\$1=\$1;print}') "
printf "  cutadapt    \$(cutadapt --version) "

mkdir trim_out
trim_galore ${input_fastq} \
    -a AAAAAAAA \
    --three_prime_clip_R1 1 \
    --gzip \
    -o ./trim_out/

cat trim_out/*trimming_report.txt | sed '/Overview of/,/RUN/{//!d}' | sed 's/Overview of removed sequences//' >> ${name}_trim.txt
cat ${name}_trim.txt
