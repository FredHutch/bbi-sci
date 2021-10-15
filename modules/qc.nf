

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


save_cds = {params.output_dir + "/" + it - ~/_cds.RDS/ - ~/temp_fold/ + "/" + it - ~/temp_fold/}
save_cell_qc = {params.output_dir + "/" + it - ~/_cell_qc.csv/ - ~/temp_fold/ + "/" + it - ~/temp_fold/}
save_samp_stats = {params.output_dir + "/" + it - ~/_sample_stats.csv/ + "/" + it}

process reformat_qc {
   container "${params.container__Rscript}"

   publishDir path: "${params.output_dir}/", saveAs: save_cds, pattern: "temp_fold/*cds.RDS", mode: 'copy'
   publishDir path: "${params.output_dir}/", saveAs: save_cell_qc, pattern: "temp_fold/*cell_qc.csv", mode: 'copy'
   publishDir path: "${params.output_dir}/", saveAs: save_samp_stats, pattern: "*sample_stats.csv", mode: 'copy'

   input:
      tuple val(key), file(scrub_csv), file(cds_object), file(cell_qc), file(dup_stats)

   output:
      set key, file("temp_fold/*.RDS"), file("temp_fold/*.csv"), emit: rscrub_out
      file("*sample_stats.csv"), emit: sample_stats
      file("*collision.txt"), emit: collision

   script:
      template "reformat_qc.R"

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

save_knee = {params.output_dir + "/" + it - ~/_knee_plot.png/ + "/" + it}
save_umap = {params.output_dir + "/" + it - ~/_UMAP.png/ + "/" + it}
save_cellqc = {params.output_dir + "/" + it - ~/_cell_qc.png/ + "/" + it}
save_garnett = {params.output_dir + "/" + it.split("_")[0] + "/" + it}

process generate_qc_metrics {
    container "${params.container__Rscript}"

    publishDir path: "${params.output_dir}/", saveAs: save_umap, pattern: "*UMAP.png", mode: 'copy'
    publishDir path: "${params.output_dir}/", saveAs: save_knee, pattern: "*knee_plot.png", mode: 'copy'
    publishDir path: "${params.output_dir}/", saveAs: save_cellqc, pattern: "*cell_qc.png", mode: 'copy'
    publishDir path: "${params.output_dir}/", saveAs: save_garnett, pattern: "*Garnett.png", mode: 'copy'

    input:
        tuple val(key), file(cds_object), file(cell_qc), file(umis_per_cell)

    output:
        file("*.png"), emit: qc_plots
        file("*.txt"), emit: cutoff

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
   container "${params.container__alpine}"

   publishDir path: "${params.output_dir}/", pattern: "all_sample_stats.csv", mode: 'copy'

   input:
      file(files)

   output:
      file "*ll_sample_stats.csv"

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
   container "${params.container__alpine}"

   input:
      file cell_qcs

   output:
      file "*.txt"

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
   container "${params.container__alpine}"

   input:
      file(col_files)

   output:
      file "*.txt"

   """

   cat $col_file > all_collision.txt

   """
}
