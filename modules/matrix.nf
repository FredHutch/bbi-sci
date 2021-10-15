

/*************

Process: make_matrix

 Inputs:
    key - sample id
    gtf_path - path to gtf info folder
    cell_gene_count - gzipped text file with a count of cell, gene pairs
    logfile - running log

 Outputs:
    key - sample id
    logfile - running log
    umi_matrix - MatrixMarket format matrix of cells by genes
    cell_anno - Cell annotations for umi_matrix
    gene_anno - Gene annotations for umi_matrix
    gtf_path - path to gtf info folder

 Pass through:

 Summary:
    Generate a matrix of cells by genes - make_matrix.py

 Downstream:
    make_cds

 Published:
    umi_matrix - MatrixMarket format matrix of cell by umi
    cell_anno - Cell annotations for umi_matrix
    gene_anno - Gene annotations for umi_matrix

 Notes:

*************/


save_umi = {params.output_dir + "/" + it - ~/.umi_counts.mtx/ + "/umi_counts.mtx"}
save_cell_anno = {params.output_dir + "/" + it - ~/.cell_annotations.txt/ + "/cell_annotations.txt"}
save_gene_anno = {params.output_dir + "/" + it - ~/.gene_annotations.txt/ + "/gene_annotations.txt"}

process make_matrix {
   container "${params.container__tools}"

   publishDir path: "${params.output_dir}/", saveAs: save_umi, pattern: "*umi_counts.mtx", mode: 'copy'
   publishDir path: "${params.output_dir}/", saveAs: save_cell_anno, pattern: "*cell_annotations.txt", mode: 'copy'
   publishDir path: "${params.output_dir}/", saveAs: save_gene_anno, pattern: "*gene_annotations.txt", mode: 'copy'

   input:
      tuple val(key), file(cell_gene_count), file(logfile), val(gtf_path)

   output:
      tuple val(key), file("*cell_annotations.txt"), file("*umi_counts.mtx"), file("*gene_annotations.txt"), val(gtf_path), file("make_matrix.log")

   script:
      template "make_matrix.sh"

}
