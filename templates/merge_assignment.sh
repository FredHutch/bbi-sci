
!#/bin/bash
set -euo pipeline

cat ${logfile} > merge_assignment.log
cat $split_bed > "${key}.bed"
sort -m -k1,1 -k2,2 $split_gene_assign > "${key}_ga.txt"

datamash -g 1,2 count 2 < "${key}_ga.txt" \
| gzip > "${key}.gz"


umi=`cat $split_umi_count | awk '{ sum += \$1 } END { print sum }'`
read=`cut -f2 $read_count`
perc=\$(echo "100.0 * (1 - \$umi/\$read)" | bc -l)
printf "%-18s   %10d    %10d    %7.1f\\n" $key \$read \$umi \$perc \
>"${key}.duplication_rate_stats.txt"

printf "
    remove_dups ending reads  : \$(wc -l ${key}.bed | awk '{print \$1;}')\n\n
    Read assignments:\n\$(awk '{count[\$3]++} END {for (word in count) { printf "            %-20s %10i\\n", word, count[word]}}' ${key}_ga.txt)\n\n" >> merge_assignment.log

printf "** End processes 'remove duplicates, assign_genes, merge_assignment' at: \$(date)\n\n" >> merge_assignment.log
