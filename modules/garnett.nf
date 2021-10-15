
process apply_garnett {
    container "${params.container__Rscript}"

    input:
        tuple val(key), file(scrub_matrix), file(cds_object), file(cell_qc), file(logfile)

    output:
        tuple val(key), file(scrub_matrix), file("new_cds/*.RDS"), file(cell_qc), file("apply_garnett.log")
    
    script:
        template "apply_garnett.sh"

}
