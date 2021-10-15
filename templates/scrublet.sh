#!/bin/bash
set -euo pipefail

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