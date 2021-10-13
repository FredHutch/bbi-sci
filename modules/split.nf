

/*************

Process: split_bam

 Inputs:
    merged_bam - sorted and quality filtered bam merged by sample
    logfile - running log

 Outputs:
    merged_bam - sorted and quality filtered bam merged by sample - stops here
    split_bam - bams split by reference (chromosome)

 Pass through:
    key - sample id
    gtf_path - path to gtf info folder
    logfile - running log

 Summary:
    Use bamtools split to split bam by reference (chromosome) for speed
    Combine the small non-chromosomal references to keep file number reasonable

 Downstream:
    merge_assignment
    remove_dups_assign_genes

 Published:

 Notes:
    Potential improvement: find a way to split bams into more evenly sized chunks

*************/

process split_bam {
    container "${params.container_tools}"

    input:
        tuple val(key), file(merged_bam), val(gtf_path)

    output:
        tuple val(key), file("split_bams/*.bam"), val(gtf_path), emit: split_bams
        tuple val(key), file("remove_dups.log"), emit: split_bam_log
        tuple file merged_bam, emit: output

    script:
        template "split.sh"

}
