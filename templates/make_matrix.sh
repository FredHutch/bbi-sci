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
