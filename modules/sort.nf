
/*************

Process: sort_and_filter

 Inputs:
    name - file id (including lane and split info)
    aligned_bam - bam output from star alignment
    logfile - running log
    cores_sf - number of cores to use

 Outputs:
    name - file id (including lane and split info)
    sorted_bam - sorted and quality filtered bam
    log_piece4 - piece of log to be concatenated for full log

 Pass through:
    key - sample id
    log_piece2 - piece of log to be concatenated for full log
    log_piece3 - piece of log to be concatenated for full log

 Summary:
    Use samtools to filter for read quality 30 and sort bam

 Downstream:
    merge_bams
    combine_logs

 Published:

 Notes:

*************/

// Cores for sort and filter set at 10 unless limit is lower
cores_sf = params.max_cores < 10 ? params.max_cores : 10

process sort_and_filter {
   container "${params.container__samtools}"

   cores_sf = params.max_cores < 10 ? params.max_cores : 10
   cpus "${cores_sf}"

   input:
      tuple val(key), val(name), file(aligned_bam)

   output:
      tuple val(key), file("*.bam"), emit: sorted_bams
      tuple val(key), file("*_sf.txt"), emit: sf_logs

   script:
      template "sort_and_filter.sh"
}
