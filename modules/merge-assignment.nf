
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
    container "${params.container__tools}"

    input:
        tuple val(key), file(split_bed), file(split_gene_assign), file(split_umi_count), file(logfile), file(read_count)

    output:
        tuple val(key), file("*.gz"), file("*_ga.txt"), file("merge_assignment.log"), emit: merge_assignment_out
        tuple val(key), file("*duplication_rate_stats.txt"), emit: duplication_rate_out
        file "*.bed", emit: temp_bed

    script:
        template "merge_assignment.sh"
}
