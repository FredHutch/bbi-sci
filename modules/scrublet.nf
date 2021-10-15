
/*************

Process: run_scrublet

 Inputs:
    key - sample id
    scrub_matrix - matrix of counts output in proper format for scrublet
    logfile - running log

 Outputs:
    key - sample id
    scrub_csv - scrublet results csv
    scrublet_png - png histogram of scrublet scores
    logfile - running log

 Pass through:
    cds_object - cds object in RDS format
    cell_qc - csv of cell quality control information

 Summary:
    Run scrublet to generate doublet scores - run_scrublet.py

 Downstream:
    calc_duplication_rate
    generate_dashboard
    finish_log

 Published:
    scrublet_png - png histogram of scrublet scores

 Notes:

*************/

save_hist = {params.output_dir + "/" + it - ~/_scrublet_hist.png/ + "/" + it}

process run_scrublet {
   container "${params.container__python}"

   publishDir path: "${params.output_dir}/", saveAs: save_hist, pattern: "*png", mode: 'copy', overwrite: true

   input:
      tuple val(key), file(scrub_matrix), file(cds_object), file(cell_qc), file(logfile)

   output:
      tuple val(key), file("*scrublet_out.csv"), file(cds_object), file(cell_qc), emit: scrublet_out
      path ("*.png"), emit: scrub_pngs
      tuple val(key), file("run_scrublet.log"), emit: pipe_log

   script:
      template "scrublet.sh"

}
