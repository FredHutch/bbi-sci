#!/bin/bash
set -euo pipefail

printf "** Start process 'align_reads' for $trimmed_fastq" >> align.log
printf "    Process versions:
          star   \$(STAR --version) " >> align.log
printf "
Reference genome information:
  \$(grep fastq_url $star_path/../*gsrc/record.out | awk '{\$1=\$2=""; print \$0}')
    FASTA download date: \$(grep fastq_download_date $star_path/../*gsrc/record.out | awk '{\$1=\$2=""; print \$0}')
    Non REF sequences removed.

  \$(grep gtf_url $star_path/../*gsrc/record.out | awk '{\$1=\$2=""; print \$0}')
    GTF download date: \$(grep gtf_download_date $star_path/../*gsrc/record.out | awk '{\$1=\$2=""; print \$0}')
    \$(grep gtf_include_biotypes $star_path/record.out | awk '{\$1=\$2=""; print \$0}')

Process output:\n" >> align.log

mkdir align_out
STAR \
    --runThreadN $cores_align \
    --genomeDir $star_path \
    --readFilesIn $trimmed_fastq \
    --readFilesCommand zcat \
    --outFileNamePrefix ./align_out/${name}  \
    --outSAMtype BAM Unsorted \
    --outSAMmultNmax 1 \
    --outSAMstrandField intronMotif


cat align_out/*Log.final.out >> align.log

printf "\n** End process 'align_reads'" >> align.log

cp align.log ${name}_align.txt