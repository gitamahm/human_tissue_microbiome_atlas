This directory will contain all the necessary files to execute the shell scripts that will enable Snakemake to parallelize various jobs on compute clusters based on the provided Snakefiles. The **`config.yaml`** included in this SIMBA directory is set up assuming that SIMBA directory is the working directory within your home directory, containing all the necessary files including databases, raw fastq reads, in-house python scripts, and other files. Thus, if you need to change the file structure in any way, you can modify the paths in **`config.yaml`** file which snakefiles (those with **.snakefile** extensions) rely on. 

The shell scripts (those with **.sh** extensions) execute the viral, bacterial and fungal branches of the pipeline. We recommend the following steps: Running the **`simba_viral.sh`** followed by **`simba_mico.sh`** and then the **`simba_fungal.sh`**. There are several parameters that are user-specific, which I will describe below in brackets. We have left examples from our own runs in these shell scripts.  Note that you need Snakemake and Conda installations. The versions we tested with were Snakemake 5.10.0 and Conda 4.14.0. Please see [conda installation guide](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) as well as [Snakemake installation guide](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html). 

For more information please see tutorials on [SLURM commands](https://login.scg.stanford.edu/tutorials/job_scripts/).
```
#SBATCH --job-name=[provide a job title]
#SBATCH --time=[provide duration for the job]
#SBATCH --ntasks=[number of tasks, specifies how many instances of your command are executed]
#SBATCH --cpus-per-task=[number of CPUs per task]
#SBATCH --mail-type=[types of job statuses to have SLURM notify you about]
#SBATCH --mail-user=[job submitter email address where job statuses can be sent to]
#SBATCH -p [compute partition to request resources from]
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
If the directory is locked, run:
```
bash simba_viral.sh unlock
```
The various rules of Snakemake will call on three in-house python scripts which are included in the **pyScripts** directory. A few of these rules that rely on in-house scripts will require the **mainEnv.yaml** to create the right virtual environment, which is also shared in this folder. Additionally, the rules rely on  **UMI-tools** (v1.0.1) ([see UMI-tools installation guide](https://umi-tools.readthedocs.io/en/latest/QUICK_START.html)) to filter and extract cell barcodes and UMIs from raw fastqs using the extract and whitelist commands. **STAR** (v2.7) is also needed to align resulting sequences against a reference library containing the human genome (GRCh38.p13) and ERCCs ([see STAR installation guide](https://github.com/alexdobin/STAR)). Finally, **BLAST** (v2.10) should also be installed ([see BLAST installation guide](https://www.ncbi.nlm.nih.gov/books/NBK569861/)). 

Finally, the rules will use the config file to find paths to various databases to BLASTn sequences against. Those databases are the following:
- **`nucVirTot`** (the viral RefSeq database used in this study)
- **`micoDB`** (what is referred to as microbeDB in this study and contains fungal and bacterial representatives). 
- **`bactRefseqSlim_nuc`** (the reduced version of the bacterial RefSeq database used in this study as an intermediate database)
- **`fungalRefseqSlim_nuc`** (the reduced version of the fungal RefSeq database used in this study as an intermediate database)
- **`ntDatabaseDir`** (the nt database used in this study)
- **`human_release34`** (needed for STAR alignments to the human genome and ERCCs)

Because these databases are too large to upload to Github, we are hosting them in google drive folder [databases](https://drive.google.com/drive/u/1/folders/1s4lG2Yq7qXH-iJhCHkvh5BvoADvCeubn). After downloading them, move them to the SIMBA directory. In addition to these databases, the unprocessed fastq reads from 10X sequencing (format: {sample}_R1_001.fastq.gz"), placed in the **`raw_reads`** directory. 

In addition to these databases, we have provided example input **`raw_reads`** and output files in google drive folder [example_files](https://drive.google.com/drive/u/1/folders/1s4lG2Yq7qXH-iJhCHkvh5BvoADvCeubn) for both viral and bacterial branches of the pipeline. The pipeline can run on a local computer, however, it is really designed for being run on computer clusters due to the computationally intensive nature of BLAST. Depending on the branch of pipeline and the input dataset, the completion of each Snakefile can take from a few hours to several weeks. The viral branch is the fastest, due to having the smallest reference databases. 
