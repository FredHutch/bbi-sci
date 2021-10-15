
/*************

Process: process_hashes

 Inputs:
    key - sample id
    input_fastq - all fastq files from params.demux_out folder
    params.hash_list

 Outputs:
    hash_log - log of the hash information
    hash_mtx - MatrixMarket matrix of hash information
    hash_cell - Cell info for matrix of hash information
    hash_hash - Hash info for matrix of hash information

 Pass through:

 Summary:
    Collect and process hash barcodes - process_hashes.py

 Downstream:

 Published:
    hash_mtx - MatrixMarket matrix of hash information
    hash_cell - Cell info for matrix of hash information
    hash_hash - Hash info for matrix of hash information

 Notes:
    Only when params.hash = true

*************/
save_hash_cell = {params.output_dir + "/" + it - ~/.hashumis_cells.txt/ + "/" + it}
save_hash_hash = {params.output_dir + "/" + it - ~/.hashumis_hashes.txt/ + "/" + it}
save_hash_mtx = {params.output_dir + "/" + it - ~/.hashumis.mtx/ + "/" + it}

process process_hashes {
    container "${params.container__tools}"

    publishDir path: "${params.output_dir}/", saveAs: save_hash_cell, pattern: "*hashumis_cells.txt", mode: 'copy', overwrite: true
    publishDir path: "${params.output_dir}/", saveAs: save_hash_hash, pattern: "*hashumis_hashes.txt", mode: 'copy', overwrite: true
    publishDir path: "${params.output_dir}/", saveAs: save_hash_mtx, pattern: "*.mtx", mode: 'copy', overwrite: true

    input:
        tuple val(key), file(input_fastq)

    output:
        file("*hash.log"), emit: hash_logs
        tuple file("*mtx"), file("*hashumis_cells.txt"), file("*hashumis_hashes.txt"), emit: hash_mats

    script:
"""
set -euo pipefail

process_hashes.py --hash_sheet $params.hash_list \
    --fastq <(zcat $input_fastq) --key $key

"""
}
