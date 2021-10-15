
/*************

Process: trim_fastqs

 Inputs:
    input_fastq - all fastq files from params.demux_out folder
    logfile - running log

 Outputs:
    trim_out - output folder from trimming - stops here
    key - sample id
    name - file id (including lane and split info)
    trimmed_fastq - trimmed, gzipped fastq
    logfile - running log
    log_piece2 - piece of log to be concatenated for full log

 Pass through:
    input_fastq - all fastq files from params.demux_out folder

 Summary:
    Trim fastqs - trim_galore
    Continue log

 Downstream:
    gather_info
    process_hashes

 Published:

 Notes:
    Only moves forward if sample is in samp_list - where params.samples comes in

*************/

process trim_fastqs {
   container "${params.container__trim_galore}"

   publishDir "${params.output_dir}/logs", pattern: '*trim.txt', mode: 'copy', overwrite: true

   input:
      file input_fastq
      val sample_list

   output:
      path("trim_out"), emit: trim_output
      tuple val(key), val(name), path("trim_out/*.fq.gz"), emit: trimmed_fastq
      path('*trim.txt'), emit: trim_log
      tuple val(key), path(input_fastq), emit: fastqs_out

   when:
      !((input_fastq.name.split(/-L[0-9]{3}/)[0].split(/\.fq.part/)[0]) in "Undetermined") &&
       ((input_fastq.name.split(/-L[0-9]{3}/)[0].split(/\.fq.part/)[0]) in sample_list)

   script:
      name = input_fastq.name.replaceAll('.fastq.gz', '')
      key = input_fastq.baseName.split(/-L[0-9]{3}/)[0].split(/\.fq.part/)[0]
      template 'trim_fastqs.sh'
}
