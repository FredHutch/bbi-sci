
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
   container "${params.container__Rscript}"

   publishDir path: "${params.output_dir}/", pattern: "exp_dash", mode: 'copy', overwrite: true

   input:
      file all_sample_stats
      file cell_counts
      file all_collision
      file plots
      file scrublet_png
      file skeleton_dash

   output:
      path exp_dash

   """
   generate_dash_data.R $all_sample_stats $params.output_dir $cell_counts $all_collision $params.garnett_file

   mkdir exp_dash
   cp -R $skeleton_dash/* exp_dash/
   mv *.png exp_dash/img/

   mv *.js exp_dash/js/

   """
}
