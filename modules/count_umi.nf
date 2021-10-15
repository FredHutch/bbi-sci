
/*************

Process: count_umis_by_sample

 Inputs:
    key - sample id
    gene_assign - text file with 3 columns: sample|cell, gene, gene type (exonic, intronic)
    logfile - running log

 Outputs:
    key - sample id
    logfile - running log
    umis_per_cell - count of umis per cell
    umis_per_cell_intronic - count of umis per cell only from intronic reads - stops here

 Pass through:
    cell_gene_count - gzipped text file with a count of cell, gene pairs

 Summary:
    calculate umis per sample - tabulate_per_cell_counts.py

 Downstream:
    make_matrix
    generate_qc_metrics

 Published:
    umis_per_cell - count of umis per cell
    umis_per_cell_intronic - count of umis per cell only from intronic reads

 Notes:

*************/

save_umi_per_cell = {params.output_dir + "/" + it - ~/.UMIs.per.cell.barcode.txt/ + "/umis_per_cell_barcode.txt"}
save_umi_per_int = {params.output_dir + "/" + it - ~/.UMIs.per.cell.barcode.intronic.txt/ + "/intronic_umis_per_cell_barcode.txt"}

process count_umis_by_sample {
    container "${params.container__python}"

    publishDir path: "${params.output_dir}/", saveAs: save_umi_per_int, pattern: "*intronic.txt", mode: 'copy'
    publishDir path: "${params.output_dir}/", saveAs: save_umi_per_cell, pattern: "*barcode.txt", mode: 'copy'

    input:
        tuple val(key), file(cell_gene_count), file(gene_assign), file(logfile)

    output:
        tuple val(key), file(cell_gene_count), file("count_umis_by_sample.log"), emit: ubss_out
        tuple val(key), file("*UMIs.per.cell.barcode.txt"), emit: umis_per_cell
        file "*UMIs.per.cell.barcode.intronic.txt", emit: umi_per_cell_intronic

    script:
        template "count_umi.sh"
}