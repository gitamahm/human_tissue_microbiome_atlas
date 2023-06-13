#!/bin/bash
#SBATCH --job-name=simba_viral.
#SBATCH --time=96:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
# send mail to this address
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=gita@stanford.edu
#SBATCH -p quake

MY_HOME=/oak/stanford/groups/quake/gita/raw
MY_HOME2=/oak/stanford/groups/quake/gita/raw/tab3-14_20210420/tab10
SNAKEFILE=$MY_HOME2/simba_viral.snakefile
CONFIGFILE=$MY_HOME/config.yaml
SLURM=$MY_HOME2/slurm/
DATE=$(date "+%Y_%m_%d_%H_%M_%S")
NJOBS=500
WAIT=120
RESTART=1
LOCALCORES=1

if [ "$1" = "unlock" ]; then
    snakemake all --snakefile $SNAKEFILE --use-conda --configfile $CONFIGFILE  --keep-target-files -j $NJOBS -w $WAIT -k --rerun-incomplete --unlock
elif [ "$1" = 'dry' ]; then
    snakemake all --snakefile $SNAKEFILE --use-conda --configfile $CONFIGFILE  --keep-target-files -j $NJOBS -w $WAIT -k --rerun-incomplete -n -p --cluster "sbatch --ntasks=1 --job-name={params.name} --cpus-per-task={threads} --partition={params.partition}  --mem={params.mem} -o $SLURM/$DATE.{params.name}.%j.log" 
elif [ "$1" = 'local' ]; then
    snakemake all --snakefile $SNAKEFILE --local-cores $LOCALCORES --use-conda --configfile $CONFIGFILE  --keep-target-files -j $NJOBS -w $WAIT -k --rerun-incomplete --restart-times $RESTART
elif [ "$1" = 'dag' ]; then
    snakemake all --snakefile $SNAKEFILE --use-conda --configfile $CONFIGFILE  ---dag | dot -Tpdf > dag.pdf
elif [ "$1" = 'delete' ]; then
    snakemake all --delete-all-output
else	
    snakemake all --snakefile $SNAKEFILE --local-cores $LOCALCORES --use-conda --configfile $CONFIGFILE --keep-target-files -j $NJOBS -w $WAIT -k --rerun-incomplete --restart-times $RESTART --cluster "sbatch --ntasks=1 --job-name={params.name} --time={params.time} --cpus-per-task={threads} --partition={params.partition}  --mem={params.mem} -o $SLURM/$DATE.{params.name}.%j.log" 
fi
