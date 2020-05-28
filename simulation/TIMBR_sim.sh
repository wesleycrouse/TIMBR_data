#!/bin/bash

mkdir -p results

#number of observations per haplotype
for NJ in 50; do
	#number of alleles
	for K in `seq 1 8`; do
		#values of alpha
		for A in 1 10; do
			#number of jobs
			for i in `seq 1 100`; do
				#variance explained
				for v in 0.5; do
					job_id="TIMBR_sim_N.J_"$NJ"_K_"$K"_alpha_"$A"_v_"$v"_job_"$i
					job_call="Rscript TIMBR_sim.R "$job_id" "$NJ" "$K" "$A" "$v" 10"
					if [ ! -f "./results/"$job_id".RData" ]; then
						sbatch -t 24:00:00 --wrap="$job_call"
					fi
				done
			done
		done
	done
done
