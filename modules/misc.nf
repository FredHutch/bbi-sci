
/*************

Process: gather_info

 Inputs:
    good_sample_sheet - corrected csv sample sheet
    key - sample id
    params.star_file
    params.gene_file

 Outputs:
    key - sample id
    star_path - path to star index folder
    gtf_path - path to gtf info folder
    star_mem - GB needed for star alignment

 Pass through:
    name - file id (including lane and split info)
    trimmed_fastq - trimmed, gzipped fastq
    log_piece2 - piece of log to be concatenated for full log
    logfile - running log

 Summary:
    Gather the star and gtf paths and info for downstream

 Downstream:
    split_bam
    make_matrix
    align_reads

 Published:

 Notes:
   o  the 'spec' variable uses the awk split() function to remove
      '_fq_part', which is added when very large fastq files are
      split. The gsub() function removes unacceptable characters
      from the sample names.

*************/

process gather_info {
    container "${params.container__python}"

    input:
        file good_sample_sheet, stageAs: "sample_sheet.csv"
        tuple val(key), val(name), file(trimmed_fastq)

    output:
        tuple val(key), val(name), env(star_path), env(star_mem), file(trimmed_fastq), emit: align_prepped
        tuple val(key), env(gtf_path), emit: gtf_info

    script:
        template 'gather_info.sh'
}