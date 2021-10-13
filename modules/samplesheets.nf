
/*************

Process: check_sample_sheet

 Inputs:
    params.sample_sheet
    params.star_file
    params.rt_barcode_file
    params.level
    params.max_wells_per_sample

 Outputs:
    good_sample_sheet - corrected csv sample sheet
    logfile - running log
    log_piece1 - piece of log to be concatenated for full log

 Pass through:

 Summary:
    Check and process sample sheet - check_sample_sheet.py
    Start log

 Downstream:
    gather_info
    trim_fastqs
    combine_logs

 Published:

 Notes:

*************/

process check_sample_sheet {
    container "${params.container__python}"

    input:
        file(sample_sheet)
        file(star_file)
        file(rt_barcode_file)

    output:
        path("*.csv")

    script:
"""
set -euo pipefail

printf "Process versions: \$(python --version)"

check_sample_sheet.py --sample_sheet $sample_sheet --star_file $star_file \
    --level $params.level --rt_barcode_file $rt_barcode_file \
    --max_wells_per_samp $params.max_wells_per_sample
"""
}
