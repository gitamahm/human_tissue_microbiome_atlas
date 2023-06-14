The shell scripts execute the viral, bacterial and fungal branches of the pipeline. We recommend the following steps: Running the simba_viral.sh followed by the .htmsimba_mico.sh and then the simba_fungal.sh. There are several parameters that are user-specific, which I will describe below in brackets. We have left examples from our own runs in these shell scripts to help with formatting your own.  
Note that in addition to Snakemake and Conda installations, there is a need for a config.yaml file that contains parameters used by Snakemake. 

For more details on SLURM commands please see documentations such as this: https://login.scg.stanford.edu/tutorials/job_scripts/ 
For more details on Snakemake commands please refer to the Snakemake tutorials: https://snakemake.readthedocs.io/en/stable/tutorial/tutorial.html 

#SBATCH --job-name=[provide a job title]
#SBATCH --time=[provide duration for the job]
#SBATCH --ntasks=[number of tasks, specifies how many instances of your command are executed]
#SBATCH --cpus-per-task=[number of CPUs per task]
#SBATCH --mail-type=[types of job statuses to have SLURM notify you about]
#SBATCH --mail-user=[job submitter email address where job statuses can be sent to]
#SBATCH -p [compute partition to request resources from]

MY_HOME=[local file path for hosting the config.yaml file]
MY_HOME2=[local file path for hosting the snakefiles]
SNAKEFILE=$MY_HOME2/simba_viral.snakefile
CONFIGFILE=$MY_HOME/config.yaml
SLURM=$MY_HOME2/slurm/
DATE=[date format]
NJOBS=[max number of jobs to send out]
WAIT=[how many seconds to wait]
RESTART=[how many times to retry if a job fails]
LOCALCORES=[number of cores to use for local jobs]

The following commands will let Snakemake either "unlock" a directory, perform a "dry" run, run "local" as opposed to using compute clusters, build a "dag" file, "delete" or run "all". 

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
