The shell scripts execute the viral, bacterial and fungal branches of the pipeline. We recommend the following steps: Running the **`simba_viral.sh`** followed by **`simba_mico.sh`** and then the **`simba_fungal.sh`**. There are several parameters that are user-specific, which I will describe below in brackets. We have left examples from our own runs in these shell scripts.  
Note that in addition to Snakemake and Conda installations, there is a need for a **`config.yaml`** file that contains parameters used by Snakemake. 

For more information on SLURM commands please see documentations such as [this](https://login.scg.stanford.edu/tutorials/job_scripts/).

```
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
```

The following commands will let Snakemake either "unlock" a directory, perform a "dry" run, run "local" as opposed to using compute clusters, build a "dag" file, "delete" or run "all". 
For more information on Snakemake commands please refer to the [Snakemake tutorials](https://snakemake.readthedocs.io/en/stable/tutorial/tutorial.html). 
```
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
```
To run the **`simba_viral.sh`** script, for example, use the following to get a dry run. 
```
bash simba_viral.sh dry
```
If there are no errors, proceed with the following which will execute the viral branch of the pipeline. Note that there are several instances such as **--partition={params.partition}** where the parameter is being set within the accompanying snakefile, in this case **`simba_viral.snakefile`**. Each rule within each of these snakefiles contains our set parameters for partition, number of threads (or cpus), job name, duration, and memory. Parameters such as **partition** are user-specific and need to be modified within these snakefiles if these files are to be executed on a computer cluster. 
```
sbatch simba_viral.sh 
```
If the directory is locked, simply run:
```
bash simba_viral.sh unlock
```
The various rules of snakemake will call on three in-house python scripts which are included in the **pyScripts** directory. A few of these rules will require the **mainEnv.yaml** to create the right virtual environment, which is also shared in this folder. Additionally, the rules rely on  [**UMI-tools** (v1.0.1)](https://umi-tools.readthedocs.io/en/latest/QUICK_START.html) to filter and extract cell barcodes and UMIs from raw fastqs using the extract and whitelist commands. [**STAR** (v2.7)](https://github.com/alexdobin/STAR) is also needed to align resulting sequences against a reference library containing the human genome (GRCh38.p13) and ERCCs. Finally, [**BLAST** (v2.10)](https://www.ncbi.nlm.nih.gov/books/NBK569861/) should also be installed. 
