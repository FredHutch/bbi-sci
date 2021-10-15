
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