import glob
from collections import defaultdict
shell.prefix("set +euo pipefail;")
REFDIR = "/oak/stanford/groups/quake/gita/raw/tab3-14_20210420/tab10"
READSDIR = REFDIR + "/raw_reads"
SCRIPTSDIR = "/oak/stanford/groups/quake/gita/raw/pyScripts"
STARREFDIR="/oak/stanford/groups/quake/gita/raw/database/STARgenome/human_release34/starIndexWithERCC/"
STAR = "STAR"
CHRNAME = STARREFDIR + "chrName.txt"
GTFFILE = "/oak/stanford/groups/quake/gita/raw/database/STARgenome/human_release34/hg38.ERCC.gtf"
#this format is 10X 
SAMPLE, = glob_wildcards(READSDIR + "/{sample}_R1_001.fastq.gz")

rule all:
	input:
              	expand(REFDIR + "/umiFiltered/{sample}_extracted_R1.fastq.gz",sample=SAMPLE),
                expand(REFDIR + "/umiFiltered/{sample}_extracted_R2.fastq.gz",sample=SAMPLE),
                expand(REFDIR + "/umiFiltered/{sample}_whitelist.txt", sample=SAMPLE),

		expand(REFDIR + "/mapped/{sample}.Aligned.sortedByCoord.out.bam",sample=SAMPLE),
                expand(REFDIR + "/mapped/{sample}.Aligned.out.bam",sample=SAMPLE),
                expand(REFDIR + "/mapped/{sample}.Unmapped.out.mate1", sample=SAMPLE),

		expand(REFDIR + '/virRefseqBlastn/{sample}.csv',sample=SAMPLE),
		expand(REFDIR + '/virRefseqBlastn/{sample}_deduplicated_filtered.csv',sample=SAMPLE),
		expand(REFDIR + '/virRefseqBlastn/{sample}_secondBlast_filtered.fasta',sample=SAMPLE),

		expand(REFDIR + '/virNTblastn/{sample}.csv',sample=SAMPLE),
		expand(REFDIR + '/virNTblastn/{sample}_deduplicated.csv',sample=SAMPLE),
		expand(REFDIR + '/virNTblastn/{sample}_secondBlast.fasta',sample=SAMPLE)

	params: 
		name="all", 
		partition="normal,quake,owners", 
		mem="4000"
	threads: 1

rule umi_whitelist:
    input: r1= READSDIR + "/{sample}_R1_001.fastq.gz"
    output: o1 = REFDIR + "/umiFiltered/{sample}_whitelist.txt"
    params: name="umi_whitelist",partition='quake,normal,owners',mem='16000',time='24:00:00'
    threads: 4
    benchmark:REFDIR + "/benchmarks/umiFiltered/{sample}.benchmark.txt"
    shell:"umi_tools whitelist --stdin {input.r1} --bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNNNN --log2stderr > {output.o1}"

#barcode pattern has 16 bps for barcode and 12 bps for UMI
rule umi_filter:
    input: r1= READSDIR + "/{sample}_R1_001.fastq.gz", r2= READSDIR + "/{sample}_R2_001.fastq.gz",r3 = REFDIR + "/umiFiltered/{sample}_whitelist.txt"
    output: o1 = REFDIR + "/umiFiltered/{sample}_extracted_R1.fastq.gz",  o2 = REFDIR + "/umiFiltered/{sample}_extracted_R2.fastq.gz"
    params: name="umi_filter",partition='quake,normal,owners',mem='16000',time='24:00:00'
    threads: 4
    benchmark:REFDIR + "/benchmarks/umiFiltered/{sample}.benchmark.txt"
    shell:"umi_tools extract --bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNNNN --stdin {input.r1} --stdout {output.o1} --read2-in {input.r2} --read2-out={output.o2} --whitelist={input.r3} --error-correct-cell"

rule align:
     #this line is different from smartseq version of the pipeline
     input: CHRNAME, R1=REFDIR + "/umiFiltered/{sample}_extracted_R2.fastq.gz",starref=STARREFDIR,gtf=GTFFILE
     output: REFDIR + "/mapped/{sample}.Aligned.sortedByCoord.out.bam", REFDIR +"/mapped/{sample}.Aligned.out.bam",REFDIR + "/mapped/{sample}.Unmapped.out.mate1"
     threads: 8
     params: name="align", mem="32000", dir="mapped/{sample}", partition="normal,quake,owners",time='03:00:00'
     benchmark:REFDIR + "/benchmarks/align/{sample}.benchmark.txt"
     shell:
       """
       {STAR}  --genomeDir {input.starref} \
               --outFileNamePrefix {params.dir}. \
               --readFilesIn {input.R1}\
               --runThreadN 8 \
               --readFilesCommand zcat \
               --outReadsUnmapped Fastx \
               --outSAMtype BAM Unsorted SortedByCoordinate \
               --outSAMattributes All \
               --outMultimapperOrder Random 
       """

rule fasta:
    input: r1= REFDIR + "/mapped/{sample}.Unmapped.out.mate1"
    output: o1 = REFDIR + "/humanFiltered/{sample}.fasta"
    params: name="fasta",partition='quake,normal,owners',mem='4000',time='00:05:00'
    threads: 1
    benchmark:REFDIR + "/benchmarks/fasta/{sample}.benchmark.txt"
    shell:"perl -ne 'y/@/>/;print($_.<>)&&<>&&<>' {input.r1} > {output.o1}"

rule poly:
    input:REFDIR + "/humanFiltered/{sample}.fasta"
    output:REFDIR + "/polyFiltered/{sample}.fasta"
    params: name="poly",partition='quake,normal,owners',mem='100000',time='12:00:00'
    threads: 16
    conda:config['mainEnv_yaml']
    benchmark:REFDIR + "/benchmarks/polyFiltered/{sample}.benchmark.txt"
    script:SCRIPTSDIR + "/poly.py"

rule virBlastn_refseq:
    input: r1 = REFDIR + "/polyFiltered/{sample}.fasta"
    output: o2 = REFDIR + '/virRefseqBlastn/{sample}.tab'
    params: name='virRefseq_blastn',partition='quake,normal,owners',mem='32000',time='24:00:00'
    threads: 8
    benchmark:REFDIR + "/benchmarks/virRefseq_blastn/{sample}.benchmark.txt"
    shell: """blastn -query {input.r1} -db {config[nucVirTotal]}  -num_threads {threads} -outfmt '6 qseqid sseqid stitle bitscore pident evalue gapopen qstart qend sstart send length mismatch staxids' -evalue 1e-5 -max_target_seqs 1 -max_hsps 1 -out {output.o2}"""

rule streamline_virRefseq_blastn:
    input:REFDIR + "/virRefseqBlastn/{sample}.tab", REFDIR + "/polyFiltered/{sample}.fasta" 
    output:REFDIR + '/virRefseqBlastn/{sample}.csv',REFDIR + '/virRefseqBlastn/{sample}_deduplicated_filtered.csv',REFDIR + '/virRefseqBlastn/{sample}_secondBlast_filtered.fasta'
    params:name='streamline_virRefseq',partition='quake,normal,owners',mem='16000', percentile='90',num_rows='5000',time='02:00:00'
    threads:4
    conda:config['mainEnv_yaml']
    benchmark:REFDIR + "/benchmarks/streamline_virRefseq_blastn/{sample}.benchmark.txt"
    script:SCRIPTSDIR + "/dedup_v5.py"

rule virBlastn_NT:
    input: r1 = REFDIR + '/virRefseqBlastn/{sample}_secondBlast_filtered.fasta'
    output: o1 = REFDIR + '/virNTblastn/{sample}.tab'
    params: name='virBlastn_NT',partition='quake,normal,owners',mem='80000',time='48:00:00'
    threads: 22
    benchmark:REFDIR + "/benchmarks/virBlastnNT/{sample}.benchmark.txt"
    shell: """blastn -query {input.r1} -db {config[ntDatabaseDir]}  -num_threads {threads} -outfmt '6 qseqid sseqid stitle bitscore pident evalue gapopen qstart qend sstart send length mismatch staxids' -evalue 1e-5 -max_target_seqs 1 -max_hsps 1 -out {output.o1}"""

#want the csv of all blast results as opposed to individual .tab file per sample. The first input is inner merged with second input just to get sequeneces and csvs. 
rule streamline_virNT_blastn:
    input: REFDIR + '/virNTblastn/{sample}.tab', REFDIR + '/virRefseqBlastn/{sample}_secondBlast_filtered.fasta'
    output:REFDIR + '/virNTblastn/{sample}.csv',REFDIR + '/virNTblastn/{sample}_deduplicated.csv',REFDIR + '/virNTblastn/{sample}_secondBlast.fasta'
    params:name='streamline_virNTblastn',partition='quake,normal,owners',mem='16000',time='02:00:00'
    threads:4
    benchmark:REFDIR + "/benchmarks/streamline_virNTblastn/{sample}.benchmark.txt"
    conda:config['mainEnv_yaml']
    script: SCRIPTSDIR + "/dedup_v5.py"

