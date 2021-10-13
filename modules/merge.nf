
/*************

 Process: merge_bams

 Inputs:
    key - sample id
    logfile - running log
    sorted_bam - sorted and quality filtered bam

 Outputs:
    key - sample id
    merged_bam - sorted and quality filtered bam merged by sample
    read_count - file with the total reads listed
    logfile - running log

 Pass through:

 Summary:
    Use samtools to merge bams from the same sample
    Count the number of reads in the sample

 Downstream:
    split_bam
    calc_duplication_rate

 Published:
    merged_bam - sorted and quality filtered bam merged by sample

 Notes:

*************/

cores_merge = params.max_cores < 8 ? params.max_cores : 8

save_bam = {params.output_dir + "/" + it - ~/.bam/ + "/" + it}

process merge_bams {
    container "${params.container__samtools}"

    cpus cores_merge
    publishDir path: "${params.output_dir}/", saveAs: save_bam, pattern: "*.bam", mode: 'copy', overwrite: true
    publishDir path: "${params.output_dir}/logs", pattern: "*.read_count.txt", mode: 'copy', overwrite: true

    input:
        tuple val(key), file(sorted_bam)

    output:
        tuple val(key), file("*.bam"), emit: sample_bams
        tuple val(key), file("*.read_count.txt"), emit: read_count

    script:
        template "merge.sh"
}
