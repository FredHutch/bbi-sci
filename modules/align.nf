
/*************

Process: align_reads

 Inputs:
    name - file id (including lane and split info)
    star_path - path to star index folder
    star_mem - GB needed for star alignment
    trimmed_fastq - trimmed, gzipped fastq
    logfile - running log
    cores_align - number of cores to use

 Outputs:
    align_out - folder of all alignment output - stops here
    name - file id (including lane and split info)
    aligned_bam - bam output from star alignment
    logfile - running log
    log_piece3 - piece of log to be concatenated for full log

 Pass through:
    key - sample id
    log_piece2 - piece of log to be concatenated for full log

 Summary:
    Align reads to genome using STAR

 Downstream:
    sort_and_filter

 Published:

 Notes:

*************/

// Cores for alignment set at 8 unless limit is lower
cores_align = params.max_cores < 8 ? params.max_cores : 8

process align_reads {
   container "${params.container__star}"

   publishDir "${params.output_dir}/logs", pattern: '*align.txt', mode: 'copy', overwrite: true

   memory { star_mem.toInteger()/cores_align + " GB" }
   cpus cores_align

   input:
      tuple val(key), val(name), val(star_path), val(star_mem), file(trimmed_fastq)

   output:
      path("align_out"), emit: align_output
      path("*align.txt"), emit: align_logs
      tuple val(key), val(name), file("align_out/*Aligned.out.bam"), emit: aligned_bams

   script:
      template 'align.sh'

}