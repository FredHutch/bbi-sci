
/*************

Process: merge_assignment

 Inputs:
    key - sample id
    split_gene_assign - text file with 3 columns: sample|cell, gene, gene type (exonic, intronic)
    split_bed - deduplicated sorted bed file
    logfile - running log
    split_umi_count - file with count of umis in split bam
    read_count - file with the total reads listed

 Outputs:
    key - sample id
    gene_assign - text file with 3 columns: sample|cell, gene, gene type (exonic, intronic)
    cell_gene_count - gzipped text file with a count of cell, gene pairs
    logfile - running log
    dup_stats - file with duplication rate information for the sample

 Pass through:

 Summary:
    merge bed files by sample
    merge gene assignment files by sample
    make cell gene count file
    calculate duplication rate

 Downstream:
    count_umis_by_sample
    reformat_qc

 Published:

 Notes:

*************/


process merge_assignment {

    input:
        tuple val(key), file(split_bed), file(split_gene_assign), file(split_umi_count), file(logfile), file(read_count) from for_cat_dups

    output:
        tuple val(key), file("*.gz"), file("*_ga.txt"), file("merge_assignment.log") into merge_assignment_out
        tuple val(key), file("*duplication_rate_stats.txt") into duplication_rate_out
        file "*.bed"  into temp_bed

    """
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

    """
}
