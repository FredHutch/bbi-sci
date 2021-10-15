
/*************

Process: make_cds

 Inputs:
    key - sample id
    umi_matrix - MatrixMarket format matrix of cells by genes
    cell_anno - Cell annotations for umi_matrix
    gene_anno - Gene annotations for umi_matrix
    gtf_path - path to gtf info folder
    logfile - running log
    params.umi_cutoff

 Outputs:
    key - sample id
    scrub_matrix - matrix of counts output in proper format for scrublet
    cds_object - cds object in RDS format
    cell_qc - csv of cell quality control information
    logfile - running log

 Pass through:

 Summary:
    Generate a monocle3 cds object - make_cds.R

 Downstream:
    run_scrublet
    calc_cell_totals

 Published:

 Notes:

*************/

process make_cds {
    container "${params.container__Rscript}"

    input:
        tuple val(key), file(cell_data), file(umi_matrix), file(gene_data), val(gtf_path), file(logfile)

    output:
        tuple val(key), file("*for_scrub.mtx"), file("*.RDS"), file("*cell_qc.csv"), file("make_cds.log"), emit: cds_out
        path("*cell_qc.csv"), emit: cell_qcs

    """
    cat ${logfile} > make_cds.log
    printf "** Start process 'make_cds' at: \$(date)\n\n" >> make_cds.log
    printf "    Process versions:
        \$(R --version | grep 'R version')
            monocle3 version \$(Rscript -e 'packageVersion("monocle3")')\n\n" >> make_cds.log


    make_cds.R \
        "$umi_matrix"\
        "$cell_data"\
        "$gene_data"\
        "${gtf_path}/latest.genes.bed"\
        "$key"\
        "$params.umi_cutoff"


    printf "** End process 'make_cds' at: \$(date)\n\n" >> make_cds.log
    """
}
