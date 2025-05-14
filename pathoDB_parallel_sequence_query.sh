#!/bin/sh
#SBATCH --job-name=pathoDB
#SBATCH -t 3-00:00              # time limit: (D-HH:MM)
#SBATCH --mem=10G            # memory per node in MB
#SBATCH --cpus-per-task=1
#SBATCH --array=1-1000%100

query=${1}
type=${2} # options are "phage" or "plasmid" or "genes"

module load mash
module load blast
module load python/3.8.1

pathoDB_parallel_sequence_query='/lustre/projects/SethCommichaux/scripts/pathoDB_parallel_sequence_query.py'

python3.8 ${pathoDB_parallel_sequence_query} --percent_identity 90 --query_coverage 95 --query_genes_found 0.8 --bootstrap ${SLURM_ARRAY_TASK_ID} -f ${query} -t ${type} 


