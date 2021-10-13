#!/bin/bash
set -euo pipefail

spec=$(sed 's/ *\$//g' sample_sheet.csv | awk 'BEGIN {FS=",";OFS=","}{split(\$2,a,"_fq_part");gsub("[_ /-]", ".", a[1]);print(\$1, a[1], \$3)}' | awk 'BEGIN {FS=","}; \$2=="$key" {print \$3}' | uniq)
star_mem=$(awk -v var="\$spec" '\$1==var {print \$3}' ${params.star_file} | uniq)
star_path=$(awk -v var="\$spec" '\$1==var {print \$2}' ${params.star_file} | uniq)
gtf_path=$(awk -v var="\$spec" '\$1==var {print \$2}' ${params.gene_file} | uniq)
