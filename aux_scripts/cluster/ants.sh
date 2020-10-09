#PBS -l nodes=1:ncpus=1
#PBS -l walltime=12:00:00
cd /local/kidzinsk/gmf/experiments
source /local/kidzinsk/miniconda3/bin/activate /local/kidzinsk/miniconda3/envs/gmf
Rscript simulation.ants.R $1 $2 $3 $4
