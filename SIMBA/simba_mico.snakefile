import glob
from collections import defaultdict
shell.prefix("set +euo pipefail;")
REFDIR = "/oak/stanford/groups/quake/gita/raw/tab3-14_20210420/tab10"
SCRIPTSDIR = "/oak/stanford/groups/quake/gita/raw/pyScripts"
#this format is 10x
SAMPLE, = glob_wildcards(REFDIR + "/polyFiltered/{sample}.fasta")

localrules: all,streamline_micoDB,streamline_blastn_micotoBactRefseq,streamline_mico_NT_blastn

rule all:
	input:
		expand(REFDIR + '/micoDBBlastn/{sample}.csv',sample=SAMPLE),
		expand(REFDIR + '/micoDBBlastn/{sample}_deduplicated.csv',sample=SAMPLE),
		expand(REFDIR + '/micoDBBlastn/{sample}_secondBlast.fasta',sample=SAMPLE),

                expand(REFDIR + '/micoAgainstRefseq_blastn/{sample}.csv',sample=SAMPLE), 
                expand(REFDIR + '/micoAgainstRefseq_blastn/{sample}_deduplicated_filtered.csv', sample=SAMPLE),
                expand(REFDIR + '/micoAgainstRefseq_blastn/{sample}_secondBlast_filtered.fasta', sample=SAMPLE),
                expand(REFDIR + '/micoAgainstRefseq_blastn/{sample}_nonBacterial.fasta',sample=SAMPLE),

		expand(REFDIR + '/micoNT_blastn/{sample}.csv',sample=SAMPLE),
		expand(REFDIR + '/micoNT_blastn/{sample}_deduplicated.csv', sample=SAMPLE),
		expand(REFDIR + '/micoNT_blastn/{sample}_secondBlast.fasta',sample=SAMPLE)

	params:name="all",partition="normal,quake,owners",mem="4000"
	threads: 1


#taking reads and blasting them against the micoDB (combination of microbeDB and mycoDB)
rule micoDB_blastn:
    input: r1 = REFDIR + "/polyFiltered/{sample}.fasta"
    output: o1 = REFDIR + '/micoDBBlastn/{sample}.tab'
    params: name='micoDB_blastn',partition='quake,normal,owners',mem="16000", time='24:00:00'
    threads: 4
    benchmark:REFDIR + "/benchmarks/micoDB_blastn/{sample}.benchmark.txt"
    shell: """blastn -query {input.r1} -db {config[micoDB]}  -num_threads {threads} -outfmt '6 qseqid sseqid stitle bitscore pident evalue gapopen qstart qend sstart send length mismatch staxids' -evalue 1e-5 -max_target_seqs 1 -out {output.o1} """


#making blastable reads ready for a second round of blast against the NT databases
rule streamline_micoDB:
    input:REFDIR + '/micoDBBlastn/{sample}.tab', REFDIR + "/polyFiltered/{sample}.fasta
    output:REFDIR + '/micoDBBlastn/{sample}.csv',REFDIR + '/micoDBBlastn/{sample}_deduplicated.csv',REFDIR + '/micoDBBlastn/{sample}_secondBlast.fasta'
    params:name='streamline_micoDB',partition='quake,normal,owners',mem='100000',time='00:30:00'
    threads: 4
    conda:config['mainEnv_yaml']
    benchmark:REFDIR + "/benchmarks/streamline_micoDB/{sample}.benchmark.txt"
    script: SCRIPTSDIR + "/dedup_v5.py"

#blasting reads that mapped to microbeDB against total NT database from NCBI
rule blastn_micotoBactRefseq:
     input: REFDIR +'/micoDBBlastn/{sample}_secondBlast.fasta'
     output:REFDIR + '/micoAgainstRefseq_blastn/{sample}_deduplicated_secondBlast_result_blastn.tab'
     params:name='bact_refseq_blastn',partition='quake,normal,owners',mem='100000',time='24:00:00'
     #group:"g1"
     threads:8
     benchmark:REFDIR + "/benchmarks/blastn_micotoBactRefseq/{sample}.benchmark.txt"
     shell:"""blastn -query {input} -db {config[bactRefseqSlim_nuc]} -num_threads {threads} -outfmt '6 qseqid sseqid stitle bitscore pident evalue gapopen qstart qend sstart send length mismatch staxids' -evalue 1e-5 -max_target_seqs 1 -out {output}"""
  
rule streamline_blastn_micotoBactRefseq:
     input:REFDIR + '/micoAgainstRefseq_blastn/{sample}_deduplicated_secondBlast_result_blastn.tab', REFDIR + '/micoDBBlastn/{sample}_secondBlast.fasta'
     output:REFDIR + '/micoAgainstRefseq_blastn/{sample}.csv', REFDIR + '/micoAgainstRefseq_blastn/{sample}_deduplicated_filtered.csv',REFDIR + '/micoAgainstRefseq_blastn/{sample}_secondBlast_filtered.fasta',REFDIR + '/micoAgainstRefseq_blastn/{sample}_nonBacterial.fasta'
     params:name='streamline_blastn_micotoBactRefseq',partition='quake,normal,owners',mem='100000',time='02:00:00'
     threads: 4
     conda:config['mainEnv_yaml']
     benchmark:REFDIR + "/benchmarks/streamline_blastn_micotoBactRefseq/{sample}.benchmark.txt"
     script: SCRIPTSDIR + "/dedup_v5_noBact.py"


#blastn of the refseq blastn results against the NT database
rule mico_NT_blastn:
    input: r1 = REFDIR + '/micoAgainstRefseq_blastn/{sample}_secondBlast_filtered.fasta'
    output: o1 = REFDIR + '/micoNT_blastn/{sample}.tab'
    params: name='mico_NT',partition='quake,normal,owners',mem='80000',time='48:00:00'
    threads: 16
    benchmark:REFDIR + "/benchmarks/micoNT_blastn/{sample}.benchmark.txt"
    shell: """blastn -query {input.r1} -db {config[ntDatabaseDir]}  -num_threads {threads} -outfmt '6 qseqid sseqid stitle bitscore pident evalue gapopen qstart qend sstart send length mismatch staxids' -evalue 1e-5 -max_target_seqs 1 -out {output.o1}"""


#no second blast is required here, just want the csv of all blast results as opposed to individual .tab file per sample
rule streamline_mico_NT_blastn:
    input: REFDIR + '/micoNT_blastn/{sample}.tab', REFDIR + "/polyFiltered/{sample}.fasta
    output:REFDIR + '/micoNT_blastn/{sample}.csv',REFDIR + '/micoNT_blastn/{sample}_deduplicated.csv',REFDIR + '/micoNT_blastn/{sample}_secondBlast.fasta'
    params:name='streamline_virNTblastn',partition='quake,normal,owners',mem='100000',time='01:00:00'
    threads:4
    conda:config['mainEnv_yaml']
    script: SCRIPTSDIR + "/dedup_v5.py"


