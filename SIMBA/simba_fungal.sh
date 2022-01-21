#!/bin/bash
#SBATCH --job-name=simba_fungal.
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
# send mail to this address
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=gitam@stanford.edu
#SBATCH -p quake

MY_HOME=/oak/stanford/groups/quake/gita/raw
MY_HOME2=/oak/stanford/groups/quake/gita/raw/tab2_20200508/tab2microbial/thirdAnalysis
SNAKEFILE=$MY_HOME2/simba_fungal.snakefile
CONFIGFILE=$MY_HOME/config.yaml
SLURM=$MY_HOME2/slurm/
DATE=$(date "+%Y_%m_%d_%H_%M_%S")
NJOBS=300
WAIT=120
RESTART=1
LOCALCORES=4

if [ "$1" = "unlock" ]; then
    snakemake all --snakefile $SNAKEFILE --use-conda --configfile $CONFIGFILE  --keep-target-files -j $NJOBS -w $WAIT -k --rerun-incomplete --unlock
elif [ "$1" = 'dry' ]; then
    snakemake all --snakefile $SNAKEFILE --use-conda --configfile $CONFIGFILE  --keep-target-files -j $NJOBS -w $WAIT -k --rerun-incomplete -n -p --cluster "sbatch --ntasks=1 --job-name={params.name} --cpus-per-task={threads} --partition={params.partition}  --mem={params.mem} -o $SLURM/$DATE.{params.name}.%j.log" 
elif [ "$1" = 'forceall' ]; then
    snakemake all --snakefile $SNAKEFILE --use-conda --configfile $CONFIGFILE  --keep-target-files -j $NJOBS -w $WAIT -k --rerun-incomplete -F  --cluster "sbatch --ntasks=1 --job-name={params.name} --cpus-per-task={threads} --partition={params.partition}  --mem={params.mem} -o $SLURM/$DATE.{params.name}.%j.log" 
elif [ "$1" = 'local' ]; then
    snakemake all --snakefile $SNAKEFILE --use-conda --configfile $CONFIGFILE  --keep-target-files -j $NJOBS -w $WAIT -k --rerun-incomplete --restart-times $RESTART
elif [ "$1" = 'dag' ]; then
    snakemake all --snakefile $SNAKEFILE --use-conda --configfile $CONFIGFILE  ---dag | dot -Tpdf > dag.pdf
elif [ "$1" = 'delete' ]; then
    snakemake all --delete-all-output
else	
    snakemake all --snakefile $SNAKEFILE  --local-cores $LOCALCORES --use-conda --configfile $CONFIGFILE --keep-target-files -j $NJOBS -w $WAIT -k --rerun-incomplete --restart-times $RESTART --cluster "sbatch --ntasks=1 --job-name={params.name} --time={params.time} --cpus-per-task={threads} --partition={params.partition}  --mem={params.mem} -o $SLURM/$DATE.{params.name}.%j.log" 
fi
