#!/bin/sh
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --job-name=cw_array_stan_snmm
#SBATCH --export=ALL
#SBATCH --array=1-2
#SBATCH --mem-per-cpu=4000
/share/apps/R-3.3.1/builddir/bin/Rscript run_bot_data_analysis_cluster.R $SLURM_ARRAY_TASK_ID