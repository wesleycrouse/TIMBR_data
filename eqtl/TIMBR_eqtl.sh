#!/bin/bash

mkdir -p results

NUM_LINES=$(wc -l < eqtl_list.txt)
PROBE_NUMBER=$(seq 1 $NUM_LINES)

for i in $PROBE_NUMBER
do
	PROBE_NAME=$(awk -v r=$i 'FNR==r {print $1;}' eqtl_list.txt)
	INTERVAL_NAME=$(awk -v r=$i 'FNR==r {print $2;}' eqtl_list.txt)
	
	job_call="Rscript TIMBR_eqtl.R "$PROBE_NAME" "$INTERVAL_NAME
	
	if [ ! -f "./results/results_crp_"$PROBE_NAME".RData" ]; then
		sbatch -t 2:00:00 --wrap="$job_call"
	fi
done
