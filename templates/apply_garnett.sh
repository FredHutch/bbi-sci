#!/bin/bash
set -euo pipefail

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
