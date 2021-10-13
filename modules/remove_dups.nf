
/*************

Process: remove_dups_assign_genes

 Inputs:
    split_bam - bams split by reference (chromosome)
    gtf_path - path to gtf info folder
    split_umi_count - file with count of umis in split bam

 Outputs:
    split_gene_assign - text file with 3 columns: sample|cell, gene, gene type (exonic, intronic)
    split_bed - deduplicated sorted bed file

 Pass through:
    key - sample id

 Summary:
    1. Remove duplicate reads using rmdup.py
    2. Use bedtools map to map the dedupped bed file to all exons with options:
        -s forced strandedness
        -f 0.95 95% of read must overlap exon
        -c 7 map the name of the gene
        -o distinct concatenate list of gene names
        -delim "|" custom delimiter
        -nonamecheck Don't error if there are different naming conventions for the chromosomes
    3. Use bedtools map to map output to gene index
        -s forced strandedness
        -f 0.95 95% of read must overlap exon
        -c 4 map the name of the cell name
        -o distinct concatenate list of gene names
        -delim "|" custom delimiter
        -nonamecheck Don't error if there are different naming conventions for the chromosomes
    4. Sort and collapse
    5. Run assign-reads-to-genes.py to deal with exon v intron

 Downstream:
    merge_assignment

 Published:

 Notes:
    Potential speed up - remove non-genic reads before sort?

*************/

process remove_dups_assign_genes {
    container "${params.container__tools}"

    input:
        tuple val(key), file(split_bam), val(gtf_path)

    output:
        tuple val(key), file("*.bed"), file("*_ga.txt"), file("*_umi_count.txt")

    script:
        template "remove_dups.sh"

}