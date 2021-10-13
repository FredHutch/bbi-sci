#!/usr/bin/env nextflow

nextflow.enable.dsl=2

DEFAULT = "default"
default_star_file = "$baseDir/bin/star_file.txt"
default_gene_file = "$baseDir/bin/gene_file.txt"
default_rt2_barcode_file = "$baseDir/bin/barcode_files/rt2.txt"
default_rt3_barcode_file = "$baseDir/bin/barcode_files/rt.txt"

// Containers used
params.container__python = "python:3.7.7"
genome_tools_container = "ghcr.io/fredhutch/docker-genome-tools:latest"
params.container__tools = genome_tools_container
params.container__trim_galore = genome_tools_container
params.container__star = "quay.io/biocontainers/star:2.6.1d--0"
params.container__samtools = "quay.io/biocontainers/samtools:1.9"

// Parse input parameters
params.help = false
params.samples = false
params.star_file = DEFAULT
params.gene_file = DEFAULT
params.umi_cutoff = 100
params.rt_barcode_file=DEFAULT
params.max_cores = 16
params.hash_list = false
params.max_wells_per_sample = 20
params.garnett_file = false
params.skip_doublet_detect = false
params.level = 3
params.demux_out = false
params.sample_sheet = false
params.output_dir = false

// NF Modules
import {
    check_sample_sheet
} from 'modules/samplesheets'

import {
    trim_fastqs
} from 'modules/trim'

import {
    process_hashes
} from 'modules/hash'

def helpMessage() {
    log.info"""
        BBI sci-RNA-seq Pipeline
        --------------------------------

        For reproducibility, please specify all parameters to a config file
        by specifying -c CONFIG_FILE.config.

        Usage:
            nextflow run bbi-sci -c CONFIG_FILE

        Help:
            --help                                     Show this message and exit.

        Required parameters (specify in your config file):
            params.output_dir = OUTPUT DIRECTORY       Output directory.
            params.sample_sheet = SAMPLE_SHEET_PATH    Sample sheet of the format described in the README.
            params.demux_out = DEMUX OUTPUT DIR        Path to the demux_out folder from the bbi-dmux run.
            params.level = 3                           2 or 3 level sci?

        Optional parameters (specify in your config file):
            params.rt_barcode_file = "default"         The path to a custom RT barcode file. If "default", default BBI barcodes will be used.
            params.max_cores = 16                      The maximum number of cores to use - fewer will be used if appropriate.
            process.maxForks = 20                      The maximum number of processes to run at the same time on the cluster.
            process.queue = "trapnell-short.q"         The queue on the cluster where the jobs should be submitted.
            params.samples = [sample1, sample2]        Add to only run certain samples from trimming on. Default is to run all.
            params.star_file = PATH/TO/FILE            File with the genome to star maps, similar to the one included with the package.
            params.gene_file = PATH/TO/FILE            File with the genome to gene model maps, similar to the one included with the package.
            params.umi_cutoff = 100                    The umi cutoff to be called a cell in matrix output.
            params.hash_list = false                   Path to a tab-delimited file with at least two columns, first the hash name and second the hash barcode sequence. Default is false to indicate no hashing.
            params.max_wells_per_sample = 20           The maximum number of wells per sample - if a sample is in more wells, the fastqs will be split then reassembled for efficiency.
            params.garnett_file = false                Path to a csv with two columns, first is the sample name, and second is a path to the Garnett classifier to be applied to that sample. Default is false - no classification.
            params.skip_doublet_detect = false         Whether to skip doublet detection, i.e. scrublet - useful for very large datasets.

        Issues? Contact hpliner@uw.edu
    """.stripIndent()
}

// Generate a sample list with fixed naming
def generate_sample_list(samp_file) {
    def samp_list = []

    for (line in samp_file.readLines()) {
        samp_list.add(line.split(",")[1])
    }

    samp_list = samp_list.collect{"$it".replaceAll(/\s/, ".").replaceAll(/_/, ".").replaceAll(/-/, ".").replaceAll(/\\//, ".")}
    samp_list.removeElement("Sample.ID")
    // filter out samples defined in sample parameter (if defined)
    if (params.samples != false) {
        samp_list = samp_list.intersect(params.samples.collect{"$it".replaceAll(/\s/, ".").replaceAll(/_/, ".").replaceAll(/-/, ".").replaceAll(/\\//, ".")})
    }
}


// TODO: How to save versions of software used
workflow {
    // check required options
    if (params.help || !params.output_dir || !params.sample_sheet || !params.demux_out) {
        helpMessage()
        exit 1
    }

    sample_sheet_file = file(params.sample_sheet)

    check_sample_sheet(sample_sheet_file, star_file, rt_barcode_file)

    sample_list = generate_sample_list(sample_sheet_file)

    sample_sheet_file = check_sample_sheet.out

    input_fastq = Channel.fromPath("${params.demux_out}/*.fastq.gz")

    // Trim fastqs (using trim galore)
    trim_fastqs(input_fastq, sample_list)

    trimmed_fastq = trim_fastqs.out.trimmed_fastq
    fastqs_out = trim_fastqs.out.fastqs_out

    // Gather the star and gtf paths and info for downstream
    gather_info(sample_sheet_file, trimmed_fastq)
    align_prepped = gather_info.out.align_prepped
    gtf_info = gather_info.out.gtf_info

    if (params.hash_list) {
        // Group fastqs for finding hash barcodes
        grouped_fastqs = fastqs_out.groupTuple()
        process_hashes(grouped_fastqs)
    }

    align_reads(align_prepped)

    aligned_bams = align_reads.out.aligned_bams

    sort_and_filter(aligned_bams)

    sorted_bams = sort_and_filter.out.sorted_bams

    // TODO: combine logs

    to_merge = sorted_bams.groupTuple()
    merge_bams(to_merge)

    sample_bams = merge_bams.out.sample_bams
    assign_prepped = sample_bams.join(gtf_info)

    split_bam(align_prepped)

    split_bams = split_bam.out.split_bams.flatten()

    remove_dups_assign_genes(split_bams)

    for_cat_dups = remove_dups_assign_genes.out
        .groupTuple()
        .join(split_bam_log)
        .join(read_count)

}




/*************

Process: count_umis_by_sample

 Inputs:
    key - sample id
    gene_assign - text file with 3 columns: sample|cell, gene, gene type (exonic, intronic)
    logfile - running log

 Outputs:
    key - sample id
    logfile - running log
    umis_per_cell - count of umis per cell
    umis_per_cell_intronic - count of umis per cell only from intronic reads - stops here

 Pass through:
    cell_gene_count - gzipped text file with a count of cell, gene pairs

 Summary:
    calculate umis per sample - tabulate_per_cell_counts.py

 Downstream:
    make_matrix
    generate_qc_metrics

 Published:
    umis_per_cell - count of umis per cell
    umis_per_cell_intronic - count of umis per cell only from intronic reads

 Notes:

*************/

save_umi_per_cell = {params.output_dir + "/" + it - ~/.UMIs.per.cell.barcode.txt/ + "/umis_per_cell_barcode.txt"}
save_umi_per_int = {params.output_dir + "/" + it - ~/.UMIs.per.cell.barcode.intronic.txt/ + "/intronic_umis_per_cell_barcode.txt"}

process count_umis_by_sample {
    cache 'lenient'
    publishDir path: "${params.output_dir}/", saveAs: save_umi_per_int, pattern: "*intronic.txt", mode: 'copy'
    publishDir path: "${params.output_dir}/", saveAs: save_umi_per_cell, pattern: "*barcode.txt", mode: 'copy'

    input:
        set val(key), file(cell_gene_count), file(gene_assign), file(logfile) from merge_assignment_out

    output:
        set key, file(cell_gene_count), file("count_umis_by_sample.log") into ubss_out
        set key, file("*UMIs.per.cell.barcode.txt") into umis_per_cell
        file "*UMIs.per.cell.barcode.intronic.txt" into umi_per_cell_intronic

    """
    cat ${logfile} > count_umis_by_sample.log
    printf "** Start process 'count_umis_by_sample' at: \$(date)\n\n" >> count_umis_by_sample.log
    printf "    Process versions:
        \$(python --version)\n\n" >> count_umis_by_sample.log
    printf "    Process command:
        tabulate_per_cell_counts.py
            --gene_assignment_files "$gene_assign"
            --all_counts_file "${key}.UMIs.per.cell.barcode.txt"
            --intron_counts_file "${key}.UMIs.per.cell.barcode.intronic.txt"\n\n"      >> count_umis_by_sample.log


    tabulate_per_cell_counts.py \
        --gene_assignment_files "$gene_assign" \
        --all_counts_file "${key}.UMIs.per.cell.barcode.txt" \
        --intron_counts_file "${key}.UMIs.per.cell.barcode.intronic.txt"


    printf "    Process stats:
        Total cells                            : \$(wc -l ${key}.UMIs.per.cell.barcode.txt | awk '{print \$1;}')
        Total cells > 100 reads                : \$(awk '\$2>100{c++} END{print c+0}' ${key}.UMIs.per.cell.barcode.txt)
        Total cells > 1000 reads               : \$(awk '\$2>1000{c++} END{print c+0}' ${key}.UMIs.per.cell.barcode.txt)
        Total reads in cells with > 100 reads  : \$(awk '\$2>100{c=c+\$2} END{print c+0}' ${key}.UMIs.per.cell.barcode.txt)\n\n" >> count_umis_by_sample.log

    printf "** End process 'count_umis_by_sample' at: \$(date)\n\n" >> count_umis_by_sample.log
    """
}


/*************

Process: make_matrix

 Inputs:
    key - sample id
    gtf_path - path to gtf info folder
    cell_gene_count - gzipped text file with a count of cell, gene pairs
    logfile - running log

 Outputs:
    key - sample id
    logfile - running log
    umi_matrix - MatrixMarket format matrix of cells by genes
    cell_anno - Cell annotations for umi_matrix
    gene_anno - Gene annotations for umi_matrix
    gtf_path - path to gtf info folder

 Pass through:

 Summary:
    Generate a matrix of cells by genes - make_matrix.py

 Downstream:
    make_cds

 Published:
    umi_matrix - MatrixMarket format matrix of cell by umi
    cell_anno - Cell annotations for umi_matrix
    gene_anno - Gene annotations for umi_matrix

 Notes:

*************/

ubss_out.join(gtf_info2).set{make_matrix_prepped}
save_umi = {params.output_dir + "/" + it - ~/.umi_counts.mtx/ + "/umi_counts.mtx"}
save_cell_anno = {params.output_dir + "/" + it - ~/.cell_annotations.txt/ + "/cell_annotations.txt"}
save_gene_anno = {params.output_dir + "/" + it - ~/.gene_annotations.txt/ + "/gene_annotations.txt"}

process make_matrix {
    cache 'lenient'
    publishDir path: "${params.output_dir}/", saveAs: save_umi, pattern: "*umi_counts.mtx", mode: 'copy'
    publishDir path: "${params.output_dir}/", saveAs: save_cell_anno, pattern: "*cell_annotations.txt", mode: 'copy'
    publishDir path: "${params.output_dir}/", saveAs: save_gene_anno, pattern: "*gene_annotations.txt", mode: 'copy'

    input:
        set key, file(cell_gene_count), file(logfile), val(gtf_path) from make_matrix_prepped

    output:
        set key, file("*cell_annotations.txt"), file("*umi_counts.mtx"), file("*gene_annotations.txt"), val(gtf_path), file("make_matrix.log") into mat_output

    """
    cat ${logfile} > make_matrix.log
    printf "** Start process 'make_matrix' at: \$(date)\n\n" >> make_matrix.log

    echo '    Process command:
        make_matrix.py <(zcat $cell_gene_count)
            --gene_annotation "${gtf_path}/latest.gene.annotations"
            --key "$key"
        cat ${gtf_path}/latest.gene.annotations > "${key}.gene_annotations.txt"  ' >> make_matrix.log


    make_matrix.py <(zcat $cell_gene_count) --gene_annotation "${gtf_path}/latest.gene.annotations" --key "$key"
    cat "${gtf_path}/latest.gene.annotations" > "${key}.gene_annotations.txt"


    printf "\n** End process 'make_matrix' at: \$(date)\n\n" >> make_matrix.log
    """

}


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
    cache 'lenient'

    input:
        set key, file(cell_data), file(umi_matrix), file(gene_data), val(gtf_path), file(logfile) from mat_output

    output:
        set key, file("*for_scrub.mtx"), file("*.RDS"), file("*cell_qc.csv"), file("make_cds.log") into cds_out
        file("*cell_qc.csv") into cell_qcs

    """
    cat ${logfile} > make_cds.log
    printf "** Start process 'make_cds' at: \$(date)\n\n" >> make_cds.log
    printf "    Process versions:
        \$(R --version | grep 'R version')
            monocle3 version \$(Rscript -e 'packageVersion("monocle3")')\n\n" >> make_cds.log
    echo '    Process command:
        make_cds.R
            "$umi_matrix"
            "$cell_data"
            "$gene_data"
            "${gtf_path}/latest.genes.bed"
            "$key"
            "$params.umi_cutoff"\n' >> make_cds.log


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


process apply_garnett {
    cache 'lenient'

    input:
        set key, file(scrub_matrix), file(cds_object), file(cell_qc), file(logfile) from cds_out

    output:
        set key, file(scrub_matrix), file("new_cds/*.RDS"), file(cell_qc), file("apply_garnett.log") into for_scrub

"""
    cat ${logfile} > apply_garnett.log
    printf "** Start process 'apply_garnett' at: \$(date)\n\n" >> apply_garnett.log
    mkdir new_cds
    echo "No Garnett classifier provided for this sample" > garnett_error.txt
    if [ $params.garnett_file == 'false' ]
    then
        cp $cds_object new_cds/
    else
        apply_garnett.R $cds_object $params.garnett_file $key
    fi

    cat garnett_error.txt >> apply_garnett.log
    printf "\n** End process 'apply_garnett' at: \$(date)\n\n" >> apply_garnett.log

"""


}


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
    cache 'lenient'
    publishDir path: "${params.output_dir}/", saveAs: save_hist, pattern: "*png", mode: 'copy'

    input:
        set key, file(scrub_matrix), file(cds_object), file(cell_qc), file(logfile) from for_scrub

    output:
        set key, file("*scrublet_out.csv"), file(cds_object), file(cell_qc) into scrublet_out
        file ("*.png") into scrub_pngs
        set key, file("run_scrublet.log") into pipe_log

    """
    cat ${logfile} > run_scrublet.log
    printf "** Start process 'run_scrublet' at: \$(date)\n\n" >> run_scrublet.log
    printf "    Process versions:
        \$(python --version)
            \$(pip freeze | grep scrublet | tr '==' ' ')\n\n" >> run_scrublet.log

    if [ $params.skip_doublet_detect == 'false' ]
    then
        run_scrublet.py --key $key --mat $scrub_matrix
        echo '    Process command:
        run_scrublet.py --key $key --mat $scrub_matrix\n'  >> run_scrublet.log
    else
        run_scrublet.py --key $key --mat $scrub_matrix --skip
        echo '    Process command:
        run_scrublet.py --key $key --mat $scrub_matrix --skip\n'  >> run_scrublet.log
        printf "    Scrublet skipped by request\n\n" >> run_scrublet.log
    fi

    printf "** End process 'run_scrublet' at: \$(date)\n\n" >> run_scrublet.log

    printf "** Start processes to generate qc metrics and dashboard at: \$(date)\n\n" >> run_scrublet.log
    """

}


/*************

Process: reformat_qc

 Inputs:
    key - sample id
    cds_object - cds object in RDS format
    cell_qc - csv of cell quality control information
    scrub_csv - scrublet results csv
    dup_stats - file with duplication rate information for the sample

 Outputs:
    key - sample id
    cds_object - cds object in RDS format
    sample_stats - csv with sample-wise statistics
    cell_qc - csv of cell quality control information
    collision - file containing collision rate if barnyard sample

 Pass through:

 Summary:
    Add scrublet info to cell_qc and cds object
    Calculate collision rate for barnyard
    Calculate sample statistics

 Downstream:
    generate_qc_metrics
    zip_up_sample_stats
    collapse_collision

 Published:
    cds_object - cds object in RDS format
    sample_stats - csv with sample-wise statistics
    cell_qc - csv of cell quality control information

 Notes:

*************/

scrublet_out.join(duplication_rate_out).set{reformat_qc_in}

save_cds = {params.output_dir + "/" + it - ~/_cds.RDS/ - ~/temp_fold/ + "/" + it - ~/temp_fold/}
save_cell_qc = {params.output_dir + "/" + it - ~/_cell_qc.csv/ - ~/temp_fold/ + "/" + it - ~/temp_fold/}
save_samp_stats = {params.output_dir + "/" + it - ~/_sample_stats.csv/ + "/" + it}

process reformat_qc {
    cache 'lenient'
    publishDir path: "${params.output_dir}/", saveAs: save_cds, pattern: "temp_fold/*cds.RDS", mode: 'copy'
    publishDir path: "${params.output_dir}/", saveAs: save_cell_qc, pattern: "temp_fold/*cell_qc.csv", mode: 'copy'
    publishDir path: "${params.output_dir}/", saveAs: save_samp_stats, pattern: "*sample_stats.csv", mode: 'copy'

    input:
        set key, file(scrub_csv), file(cds_object), file(cell_qc), file(dup_stats) from reformat_qc_in

    output:
        set key, file("temp_fold/*.RDS"), file("temp_fold/*.csv") into rscrub_out
        file("*sample_stats.csv") into sample_stats
        file("*collision.txt") into collision


    """
    #!/usr/bin/env Rscript

    library(monocle3)

    dir.create("temp_fold")
    cds <- readRDS("$cds_object")
    cell_qc <- read.csv("$cell_qc")

    if(nrow(pData(cds)) > 0) {
        if("$params.skip_doublet_detect" == 'false') {
            scrublet_out <- read.csv("$scrub_csv", header=F)
            pData(cds)\$scrublet_score <- scrublet_out\$V1
            pData(cds)\$scrublet_call <- ifelse(scrublet_out\$V2 == 1, "Doublet", "Singlet")
            cell_qc\$scrublet_score <- scrublet_out\$V1
            cell_qc\$scrublet_call <- ifelse(scrublet_out\$V2 == 1, "Doublet", "Singlet")
        }
    }

    write.csv(cell_qc, quote=FALSE, file="temp_fold/$cell_qc")

    dup_stats <- read.table(paste0("$key", ".duplication_rate_stats.txt"))

    df <- data.frame(sample="$key", n.reads = dup_stats\$V2, n.umi = dup_stats\$V3,
                     duplication_rate = dup_stats\$V4,
                     doublet_count = sum(cell_qc\$scrublet_call == "Doublet", na.rm=TRUE),
                     doublet_perc = paste0(round(sum(cell_qc\$scrublet_call == "Doublet",
                                                     na.rm=TRUE)/nrow(cell_qc) * 100, 1), "%"),
                     doublet_NAs=sum(is.na(cell_qc\$scrublet_call)))

    write.csv(df, file=paste0("$key", "_sample_stats.csv"), quote=FALSE, row.names=FALSE)
    saveRDS(cds, file="temp_fold/$cds_object")

    if ("$key" == "Barnyard") {
        fData(cds)\$mouse <- grepl("ENSMUSG", fData(cds)\$id)
        fData(cds)\$human <- grepl("ENSG", fData(cds)\$id)

        pData(cds)\$mouse_reads <- Matrix::colSums(exprs(cds)[fData(cds)\$mouse,])
        pData(cds)\$human_reads <- Matrix::colSums(exprs(cds)[fData(cds)\$human,])
        pData(cds)\$total_reads <- pData(cds)\$mouse_reads + pData(cds)\$human_reads
        pData(cds)\$human_perc <- pData(cds)\$human_reads/pData(cds)\$total_reads
        pData(cds)\$mouse_perc <- pData(cds)\$mouse_reads/pData(cds)\$total_reads
        pData(cds)\$collision <- ifelse(pData(cds)\$human_perc >= .9 | pData(cds)\$mouse_perc >= .9, FALSE, TRUE)

        collision_rate <- round(sum(pData(cds)\$collision/nrow(pData(cds))) * 200, 1)
        fileConn<-file("Barn_collision.txt")
        writeLines(paste0("$key", "\t", collision_rate, "%"), fileConn)
        close(fileConn)
    } else {
        fileConn<-file("${key}_no_collision.txt")
        writeLines(paste0("$key", "\t", "NA"), fileConn)
        close(fileConn)
    }
    """

}


/*************

Process: generate_qc_metrics

 Inputs:
    key - sample id
    params.umi_cutoff
    cds_object - cds object in RDS format
    cell_qc - csv of cell quality control information
    umis_per_cell - count of umis per cell

 Outputs:
    cutoff - not currently used
    umap_png - png sample UMAP
    knee_png - png sample knee plot
    qc_png - png of cell qc stats

 Pass through:

 Summary:
    Generate a bunch of qc metrics and plots - generate_qc.R

 Downstream:
    generate_dashboard

 Published:
    umap_png - png sample UMAP
    knee_png - png sample knee plot
    qc_png - png of cell qc stats

 Notes:
    Need to test umi cutoff here and in cds function
    Need to either remove or use output cutoff

*************/

for_gen_qc = rscrub_out.join(umis_per_cell)
save_knee = {params.output_dir + "/" + it - ~/_knee_plot.png/ + "/" + it}
save_umap = {params.output_dir + "/" + it - ~/_UMAP.png/ + "/" + it}
save_cellqc = {params.output_dir + "/" + it - ~/_cell_qc.png/ + "/" + it}
save_garnett = {params.output_dir + "/" + it.split("_")[0] + "/" + it}

process generate_qc_metrics {
    cache 'lenient'
    publishDir path: "${params.output_dir}/", saveAs: save_umap, pattern: "*UMAP.png", mode: 'copy'
    publishDir path: "${params.output_dir}/", saveAs: save_knee, pattern: "*knee_plot.png", mode: 'copy'
    publishDir path: "${params.output_dir}/", saveAs: save_cellqc, pattern: "*cell_qc.png", mode: 'copy'
    publishDir path: "${params.output_dir}/", saveAs: save_garnett, pattern: "*Garnett.png", mode: 'copy'

    input:
        set key, file(cds_object), file(cell_qc), file(umis_per_cell) from for_gen_qc

    output:
        file("*.png") into qc_plots
        file("*.txt") into cutoff

    """
    generate_qc.R\
        $cds_object $umis_per_cell $key \
        --specify_cutoff $params.umi_cutoff\

    """
}

/*************

Process: zip_up_sample_stats

 Inputs:
    sample_stats - csv with sample-wise statistics - collected

 Outputs:
    all_sample_stats - concatenated table of sample stats from all samples

 Pass through:

 Summary:
    Concatenate duplication information into one table

 Downstream:
    generate_dashboard

 Published:
    all_sample_stats - concatenated table of sample stats from all samples

 Notes:

*************/

process zip_up_sample_stats {
    cache 'lenient'
    publishDir path: "${params.output_dir}/", pattern: "all_sample_stats.csv", mode: 'copy'

    input:
        file files from sample_stats.collect()

    output:
        file "*ll_sample_stats.csv" into all_sample_stats

    """

     sed -s 1d $files > all_sample_stats.csv

    """
}


/*************

Process: calc_cell_totals

 Inputs:
    cell_qc - csv of cell quality control information - collected

 Outputs:
    cell_counts - table cell totals above set UMI thresholds for all samples

 Pass through:

 Summary:
    Count cell totals above set UMI thresholds for all samples

 Downstream:
    generate_dashboard

 Published:

 Notes:

*************/

process calc_cell_totals {
    cache 'lenient'

    input:
        file cell_qc from cell_qcs.collect()

    output:
        file "*.txt" into cell_counts

    """

    for f in $cell_qc
    do
      awk 'BEGIN {FS=","}; \$2>100{c++} END{print FILENAME, "100", c-1}' \$f >> cell_counts.txt
      awk 'BEGIN {FS=","}; \$2>500{c++} END{print FILENAME, "500", c-1}' \$f >> cell_counts.txt
      awk 'BEGIN {FS=","}; \$2>1000{c++} END{print FILENAME, "1000", c-1}' \$f >> cell_counts.txt
    done

    """

}


/*************

Process: collapse_collision

 Inputs:
    collision - file containing collision rate if barnyard sample - collected

 Outputs:
    all_collision - concatenate collision values for all samples (all NA except Barnyard)

 Pass through:

 Summary:
    Concatenate collision values for all samples

 Downstream:
    generate_dashboard

 Published:

 Notes:

*************/

process collapse_collision {
    cache 'lenient'

    input:
        file col_file from collision.collect()

    output:
        file "*.txt" into all_collision

    """

    cat $col_file > all_collision.txt

    """
}


/*************

Process: generate_dashboard

 Inputs:
    cell_counts - table cell totals above set UMI thresholds for all samples
    all_collision - concatenate collision values for all samples (all NA except Barnyard)
    all_sample_stats - concatenated table of sample stats from all samples
    params.output_dir
    umap_png - png sample UMAP - combined as qc_plots
    knee_png - png sample knee plot - combined as qc_plots
    qc_png - png of cell qc stats - combined as qc_plots
    scrublet_png - png histogram of scrublet scores
    params.garnett_file

 Outputs:
    exp_dash - experimental dashboard

 Pass through:

 Summary:
    Collect plots and generate data file for experimental dashboard
    Assemble dashboard

 Downstream:
    generate_summary_log

 Published:
    exp_dash - experimental dashboard

 Notes:

*************/

process generate_dashboard {
    cache 'lenient'
    publishDir path: "${params.output_dir}/", pattern: "exp_dash", mode: 'copy'

    input:
        file all_sample_stats
        file cell_counts
        file all_collision
        file plots from qc_plots.collect()
        file scrublet_png from scrub_pngs.collect()

    output:
        file exp_dash into exp_dash_out

    """
    generate_dash_data.R $all_sample_stats $params.output_dir $cell_counts $all_collision $params.garnett_file

    mkdir exp_dash
    cp -R $baseDir/bin/skeleton_dash/* exp_dash/
    mv *.png exp_dash/img/

    mv *.js exp_dash/js/

    """
}


/*************

Process: finish_log

 Inputs:
    key - sample id
    logfile - running log
    exp_dash - experimental dashboard - input so runs last

 Outputs:
    full_log - Final full pipeline log
    summary_log - Summary log
    log_data - Logging info for dashboards

 Pass through:

 Summary:
    Add parameter info to front of pipeline - allows restart when changing minor parameters
    Generate summary log
    Generate log info for dashboards

 Downstream:
    zip_up_log_data

 Published:
    full_log - Final full pipeline log
    summary_log - Summary log
    log_data - Logging info for dashboards

 Notes:

*************/

save_logs = {params.output_dir + "/" + it - ~/_read_metrics.log/ - ~/_full.log/ + "/" + it}
save_txt_for_wrap = {params.output_dir + "/" + it - ~/_log_data.txt/ + "/" + it}

process finish_log {
    cache 'lenient'
    publishDir path: "${params.output_dir}/", saveAs: save_logs, pattern: "*.log", mode: 'copy'
    publishDir path: "${params.output_dir}/", saveAs: save_txt_for_wrap, pattern: "*.txt", mode: 'copy'

    input:
        set key, file(logfile) from pipe_log
        file exp_dash from exp_dash_out

    output:
        file("*_full.log") into full_log
        file("*_read_metrics.log") into summary_log
        file("*log_data.txt") into log_txt_for_wrap

    """
    head -n 2 ${logfile} > ${key}_full.log
    printf "Nextflow version: $nextflow.version\n" >> ${key}_full.log
    printf "Pipeline version: $workflow.manifest.version\n" >> ${key}_full.log
    printf "Git Repository, Version, Commit ID, Session ID: $workflow.repository, $workflow.revision, $workflow.commitId, $workflow.sessionId\n\n" >> ${key}_full.log
    printf "Command:\n$workflow.commandLine\n\n" >> ${key}_full.log
    printf "***** PARAMETERS *****: \n\n" >> ${key}_full.log
    printf "    params.run_dir:               $params.run_dir\n" >> ${key}_full.log
    printf "    params.output_dir:            $params.output_dir\n" >> ${key}_full.log
    printf "    params.sample_sheet:          $params.sample_sheet\n" >> ${key}_full.log
    printf "    params.demux_out:             $params.demux_out\n" >> ${key}_full.log
    printf "    params.level:                 $params.level\n" >> ${key}_full.log
    printf "    params.max_cores:             $params.max_cores\n" >> ${key}_full.log
    printf "    params.samples:               $params.samples\n" >> ${key}_full.log
    printf "    params.star_file:             $params.star_file\n" >> ${key}_full.log
    printf "    params.gene_file:             $params.gene_file\n" >> ${key}_full.log
    printf "    params.umi_cutoff:            $params.umi_cutoff\n" >> ${key}_full.log
    printf "    params.rt_barcode_file:       $params.rt_barcode_file\n" >> ${key}_full.log
    printf "    params.hash_list:             $params.hash_list\n" >> ${key}_full.log
    printf "    params.max_wells_per_sample:  $params.max_wells_per_sample\n\n" >> ${key}_full.log
    printf "    params.garnett_file:          $params.garnett_file\n\n" >> ${key}_full.log
    printf "    params.skip_doublet_detect:   $params.skip_doublet_detect\n\n" >> ${key}_full.log

    tail -n +2 ${logfile} >> ${key}_full.log
    printf "\n** End processes generate qc metrics and dashboard at: \$(date)\n\n" >> ${key}_full.log
    printf "***** END PIPELINE *****: \n\n" >> ${key}_full.log
    filename=${key}_full.log

    # Trimming:
    trim_start=`cat \$filename | grep 'sequences processed in total' | awk -F ' ' '{sum += \$1} END {print sum}'`
    trim_lost=`cat \$filename | grep 'Sequences removed because they became shorter' | awk -F ' ' '{sum += \$14} END {print sum}'`
    trim_end=\$((\$trim_start - \$trim_lost))

    # Alignment:
    align_start=`cat \$filename | grep 'Number of input reads' | awk -F '|' '{sum += \$2} END {print sum}'`
    align_mapped=`cat \$filename | grep 'Uniquely mapped reads number' | awk -F '|' '{sum += \$2} END {print sum}'`
    align_totals=(\$(cat \$filename | grep 'Number of input reads' | cut -d "|" -f 2 | awk '{print \$1}'))
    align_multimapped=`cat \$filename | grep 'Number of reads mapped to multiple loci' |  awk -F '|' '{sum += \$2} END {print sum}'`
    align_too_short_arr=(\$(cat \$filename | grep 'unmapped: too short' | cut -d "|" -f 2 | tr '%' ' ' | awk '{\$1=\$1/100;print}'))
    align_too_short=`a=0
    for i in \${align_too_short_arr[@]}
    do
        echo "\${align_too_short_arr[\$a]} * \${align_totals[\$a]}" | bc
        a=\$((a+1))
    done | awk '{sum += \$1} END {printf "%1.0f", sum}'`

    # Sort and Filter:
    sf_start=`cat \$filename | grep 'sort_and_filter starting reads' | awk -F ':' '{sum += \$2} END {print sum}'`
    sf_end=`cat \$filename | grep 'sort_and_filter ending reads' | awk -F ':' '{sum += \$2} END {print sum}'`

    # Dups:
    dup_start=`cat \$filename | grep 'remove_dups starting reads' | awk -F ':' '{sum += \$2} END {print sum}'`
    dup_end=`cat \$filename | grep 'remove_dups ending reads' | awk -F ':' '{sum += \$2} END {print sum}'`

    # Assignment:
    assigned_exonic=`cat \$filename | grep '    exonic     ' | awk -F ' ' '{sum += \$2} END {print sum}'`
    assigned_intronic=`cat \$filename | grep '    intronic     ' | awk -F ' ' '{sum += \$2} END {print sum}'`
    assigned_end=\$((\$assigned_exonic + \$assigned_intronic))

    # In real cells:
    reads_in_cells=`cat \$filename | grep 'Total reads in cells with > 100 reads' | awk -F ':' '{sum += \$2} END {print sum}'`

    printf "
            \\"${key}\\": {
            \\"sample\\": \\"${key}\\",
            \\"alignment_start\\" : \\"\$align_start\\",
            \\"alignment_mapped\\" : \\"\$align_mapped\\",
            \\"align_multimapped\\" : \\"\$align_multimapped\\",
            \\"align_too_short\\" : \\"\$align_too_short\\",
            \\"sf_start\\" : \\"\$sf_start\\",
            \\"sf_end\\" : \\"\$sf_end\\",
            \\"dup_start\\" : \\"\$dup_start\\",
            \\"dup_end\\" : \\"\$dup_end\\",
            \\"assigned_exonic\\" : \\"\$assigned_exonic\\",
            \\"assigned_intronic\\" : \\"\$assigned_intronic\\",
            \\"reads_in_cells\\" : \\"\$reads_in_cells\\" }
      " > ${key}_log_data.txt


    printf "***** PIPELINE READ STATS *****: \n\n" >> ${key}_read_metrics.log

    printf "%20s %20s %20s %20s %20s\n" "Process" "Starting reads" "Ending reads" "% lost" "% of total lost" >> ${key}_read_metrics.log
    printf "========================================================================================================\n" >> ${key}_read_metrics.log
    printf "%20s %20s %20s %20.2f %20.2f\n" "Trimming" \$trim_start \$trim_end \$(echo "(\$trim_start - \$trim_end)/\$trim_start * 100" | bc -l ) \$(echo "(\$trim_start - \$trim_end)/\$trim_start * 100" | bc -l ) >> ${key}_read_metrics.log
    printf "%20s %20s %20s %20.2f %20.2f\n" "Alignment" \$align_start \$sf_start \$(echo "(\$align_start - \$sf_start)/\$align_start * 100" | bc -l ) \$(echo "(\$align_start - \$sf_start)/\$trim_start * 100" | bc -l ) >> ${key}_read_metrics.log
    printf "%20s %20s %20s %20.2f %20.2f\n" "Filtering" \$sf_start \$sf_end \$(echo "(\$sf_start - \$sf_end)/\$sf_start * 100" | bc -l ) \$(echo "(\$sf_start - \$sf_end)/\$trim_start * 100" | bc -l ) >> ${key}_read_metrics.log
    printf "%20s %20s %20s %20.2f %20.2f\n" "Deduplication" \$dup_start \$dup_end \$(echo "(\$dup_start - \$dup_end)/\$dup_start * 100" | bc -l ) \$(echo "(\$dup_start - \$dup_end)/\$trim_start * 100" | bc -l ) >> ${key}_read_metrics.log
    printf "%20s %20s %20s %20.2f %20.2f\n" "Gene assignment" \$dup_end \$assigned_end \$(echo "(\$dup_end - \$assigned_end)/\$dup_end * 100" | bc -l ) \$(echo "(\$dup_end - \$assigned_end)/\$trim_start * 100" | bc -l ) >> ${key}_read_metrics.log

    printf "\nAlignment details: \n" >> ${key}_read_metrics.log
    printf "%25s %20s %20s\n" "" "Count" "Percent"  >> ${key}_read_metrics.log
    printf "========================================================================================================\n" >> ${key}_read_metrics.log
    printf "%25s %20s %20s\n" "Total reads processed:" \$align_start "" >> ${key}_read_metrics.log
    printf "%25s %20s %20.2f\n" "Reads uniquely mapped:" \$align_mapped \$(echo "(\$align_mapped)/\$align_start * 100" | bc -l ) >> ${key}_read_metrics.log
    printf "%25s %20s %20.2f\n" "Reads multi-mapped:" \$align_multimapped \$(echo "(\$align_multimapped)/\$align_start * 100" | bc -l )  >> ${key}_read_metrics.log
    printf "%25s %20s %20.2f\n" "Reads too short:" \$align_too_short \$(echo "(\$align_too_short)/\$align_start * 100" | bc -l ) >> ${key}_read_metrics.log


    cat ${key}_read_metrics.log >> ${key}_full.log

    """

}

/*************

Process: zip_up_log_data

 Inputs:
    summary_log - collected summary log files
    full_log - collected full log files

 Outputs:
    log_js - Logging info in a js format for dashboards

 Pass through:

 Summary:
    Generate log data js file for dashboard

 Downstream:
    End

 Published:
    log_data.js - Logging info for dashboards

 Notes:

*************/

process zip_up_log_data {
    cache 'lenient'
    publishDir path: "${params.output_dir}/exp_dash/js/", pattern: "*.js", mode: 'copy'

    input:
        file summary_log from summary_log.collect()
        file full_log from full_log.collect()

    output:
        file "*.js" into log_js

    """

    echo 'const log_data = {' > log_data.js
    for file in $summary_log
    do
        samp_name=\$(basename \$file | sed 's/_read_metrics.log//')
        echo "\\"\$samp_name\\" :  \\`" >> log_data.js
        cat \$file >> log_data.js
        echo "\\`," >> log_data.js
    done
    sed -i '\$ s/,\$//' log_data.js
    echo '}' >> log_data.js

    echo 'const full_log_data = {' >> log_data.js
    for file in $full_log
    do
        samp_name=\$(basename \$file | sed 's/_full.log//')
        echo "\\"\$samp_name\\" :  \\`" >> log_data.js
        cat \$file >> log_data.js
        echo "\\`," >> log_data.js
    done
    sed -i '\$ s/,\$//' log_data.js
    echo '}' >> log_data.js

    """
}

workflow.onComplete {
	println ( workflow.success ? "Done! Saving output" : "Oops .. something went wrong" )
}
