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
params.container__Rscript = genome_tools_container
params.container__basic = "ubuntu:21.04"

// Parse input parameters
params.help = false
params.samples = false
params.star_file = DEFAULT
params.gene_file = DEFAULT
params.umi_cutoff = 100
params.rt_barcode_file = DEFAULT
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
include {
    check_sample_sheet
} from './modules/samplesheets'

include {
    trim_fastqs
} from './modules/trim'

include {
    process_hashes
} from './modules/hash'

include {
    gather_info
} from './modules/misc'

include {
    align_reads
} from './modules/align'

include {
    sort_and_filter
} from './modules/sort'

include {
    merge_bams
} from './modules/merge'

include {
    split_bam
} from './modules/split'

include {
    remove_dups_assign_genes
} from './modules/remove_dups'

include {
    merge_assignment
} from './modules/merge-assignment'

include {
    count_umis_by_sample
} from './modules/count_umi'

include {
    make_matrix
} from './modules/matrix'

include {
    make_cds
} from './modules/cds'

include {
    apply_garnett
} from './modules/garnett'

include {
    run_scrublet
} from './modules/scrublet'

include {
    reformat_qc
    generate_qc_metrics
    zip_up_sample_stats
    calc_cell_totals
    collapse_collision
} from './modules/qc'

include {
    generate_dashboard
} from './modules/dashboard'

include {
    finish_log
    zip_up_log_data
} from './modules/log'

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
def generateSampleList(samp_file) {
    def samp_list = []

    for (line in samp_file.readLines()) {
        samp_list.add(line.split(",")[1])
    }

    samp_list = samp_list.collect{"$it".replaceAll(/\s/, ".").replaceAll(/_/, ".").replaceAll(/-/, ".").replaceAll(/\\//, ".")}
    samp_list.removeElement("Sample.ID")
    // filter out samples defined in sample parameter (if defined)
    if (params.samples) {
        samp_list = samp_list.intersect(params.samples.collect{"$it".replaceAll(/\s/, ".").replaceAll(/_/, ".").replaceAll(/-/, ".").replaceAll(/\\//, ".")})
        println "Using samples: ${samp_list}"
    }
    return samp_list
}

def getRtBarcodeFile() {
    if (params.rt_barcode_file == DEFAULT) {
        if (params.level == 2) {
            return file(default_rt2_barcode_file)
        }
        if (params.level == 3) {
            return file(default_rt3_barcode_file)
        }
    } else {
        return file(params.rt_barcode_file)
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
    star_file = params.star_file == DEFAULT ? file(default_star_file) : file(params.star_file)
    gene_file = params.gene_file == DEFAULT ? file(default_gene_file) : file(params.gene_file)
    rt_barcode_file = getRtBarcodeFile()

    input_fastq = Channel.fromPath("${params.demux_out}/*.fastq.gz", checkIfExists: true)
    input_fastq.ifEmpty("No fastq files found in ${params.demux_out}")

    sample_list = generateSampleList(sample_sheet_file)

    check_sample_sheet(sample_sheet_file, star_file, rt_barcode_file)
    sample_sheet_file = check_sample_sheet.out

    // Trim fastqs (using trim galore)
    trim_fastqs(input_fastq, sample_list)

    trimmed_fastq = trim_fastqs.out.trimmed_fastq
    fastqs_out = trim_fastqs.out.fastqs_out

    // Gather the star and gtf paths and info for downstream
    gather_info(sample_sheet_file, trimmed_fastq, star_file, gene_file)
  
    // Is there a better way to stage the genome file?
    align_prepped = gather_info.out.align_prepped.map{
        it[2] = file(it[2])
        return it
    }
    gtf_info = gather_info.out.gtf_info
    align_prepped.view()
    gtf_info.view()

    if (params.hash_list) {
        // Group fastqs for finding hash barcodes
        grouped_fastqs = fastqs_out.groupTuple()
        process_hashes(grouped_fastqs)
    }


    align_reads(align_prepped)

    // aligned_bams = align_reads.out.aligned_bams

    // sort_and_filter(aligned_bams)

    // sorted_bams = sort_and_filter.out.sorted_bams

    // // TODO: combine logs

    // to_merge = sorted_bams.groupTuple()
    // merge_bams(to_merge)

    // sample_bams = merge_bams.out.sample_bams
    // assign_prepped = sample_bams.join(gtf_info)

    // split_bam(align_prepped)

    // split_bams = split_bam.out.split_bams.flatten()

    // remove_dups_assign_genes(split_bams)

    // for_cat_dups = remove_dups_assign_genes.out
    //     .groupTuple()
    //     .join(split_bam.out.split_bam_log)
    //     .join(merge_bams.out.read_count)
    
    // merge_assignment(for_cat_dups)

    // merge_assignment_out = merge_assignment.out.merge_assignment_out

    // count_umis_by_sample(merge_assignment_out)

    // make_matrix_prepped = count_umis_by_sample.out.ubss_out.join(gtf_info)

    // make_matrix(make_matrix_prepped)

    // make_cds(make_matrix.out)
    
    // apply_garnett(make_cds.out.cds_out)

    // run_scrublet(apply_garnett.out)

    // reformat_qc_in = run_scrublet.out.scrublet_out.join(merge_assignment.out.duplication_rate)
    // reformat_qc(reformat_qc_in)

    // sample_stats = reformat_qc.out.sample_stats.collect()
    // collision = reformat_qc.out.collision.collect()
    // for_gen_qc = reformat_qc.out.rscrub.join(count_umis_by_sample.out.umis_per_cell)

    // generate_qc_metrics(for_gen_qc)

    // zip_up_sample_stats(sample_stats)

    // cell_qcs = make_cds.out.cell_qcs.collect()

    // calc_cell_totals(cell_qcs)

    // collapse_collision(collision)

    // skeleton_dash = file("$baseDir/bin/skeleton_dash")
    // generate_dashboard(
    //     zip_up_sample_stats.out,
    //     calc_cell_totals.out,
    //     collapse_collision.out,
    //     generate_qc_metrics.out.qc_plots.collect(),
    //     run_scrublet.out.scrub_pngs.collect(),
    //     skeleton_dash
    // )

    // // TODO: figure out what to do with the log stuff
    // // Could be saved into a json file instead of a log
    // finish_log(
    //     run_scrublet.out.pipe_log,
    //     generate_dashboard.out
    // )

    // zip_up_log_data(
    //     finish_log.out.summary_log.collect(),
    //     finish_log.out.full_log.collect()
    // )

}

workflow.onComplete {
	println ( workflow.success ? "Done! Saving output" : "Oops .. something went wrong" )
}
