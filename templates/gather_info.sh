#!/bin/bash
set -euo pipefail

spec=\$(sed 's/ *\$//g' ${good_sample_sheet} | awk 'BEGIN {FS=",";OFS=","}{split(\$2,a,"_fq_part");gsub("[_ /-]", ".", a[1]);print(\$1, a[1], \$3)}' | awk 'BEGIN {FS=","}; \$2=="$key" {print \$3}' | uniq)
star_mem=\$(awk -v var="\$spec" '\$1==var {print \$3}' ${star_file} | uniq)
star_path=\$(awk -v var="\$spec" '\$1==var {print \$2}' ${star_file} | uniq)
gtf_path=\$(awk -v var="\$spec" '\$1==var {print \$2}' ${gene_file} | uniq)
