import glob
from collections import defaultdict
shell.prefix("set +euo pipefail;")
REFDIR = "/oak/stanford/groups/quake/gita/raw/tab2_20200508/tab2microbial/thirdAnalysis"
SCRIPTSDIR = "/oak/stanford/groups/quake/gita/raw/pyScripts"
SAMPLE, = glob_wildcards(REFDIR + '/micoAgainstRefseq_blastn/{sample}_nonBacterial.fasta')

rule all:
	input:
		expand(REFDIR + '/nonBactMicoAgainstFungRefseq_blastn/{sample}.csv',sample=SAMPLE),
		expand(REFDIR + '/nonBactMicoAgainstFungRefseq_blastn/{sample}_deduplicated.csv', sample=SAMPLE),
		expand(REFDIR + '/nonBactMicoAgainstFungRefseq_blastn/{sample}_secondBlast.fasta', sample=SAMPLE),

		expand(REFDIR + '/fungi_NT_blastn/{sample}.csv', sample=SAMPLE),
		expand(REFDIR + '/fungi_NT_blastn/{sample}_deduplicated.csv',sample=SAMPLE),
		expand(REFDIR + '/fungi_NT_blastn/{sample}_secondBlast.fasta',sample=SAMPLE)
	params: 
		name="all", 
		partition="normal,quake,owners", 
		mem="4000"
	threads: 1
 
rule nonBactMicoAgainstFungRefseq_blastn:
     input:REFDIR + '/micoAgainstRefseq_blastn/{sample}_nonBacterial.fasta'
     output:REFDIR + '/nonBactMicoAgainstFungRefseq_blastn/{sample}_deduplicated_secondBlast_result_blastn.tab'
     params:name='nonBactMicoAgainstFungRefseq_blastn',partition='quake,normal,owners',mem='8000',time='04:00:00'
     threads:4
     benchmark:REFDIR + "/benchmarks/nonBactMicoAgainstFungRefseq_blastn/{sample}.benchmark.txt"
     shell:"""blastn -query {input} -db {config[fungalRefseqSlim_nuc]} -num_threads {threads} -outfmt '6 qseqid sseqid stitle bitscore pident evalue gapopen qstart qend sstart send length mismatch staxids' -evalue 1e-5 -max_target_seqs 1 -out {output}"""

rule streamline_micotoFungRefseq:
     input:REFDIR + '/nonBactMicoAgainstFungRefseq_blastn/{sample}_deduplicated_secondBlast_result_blastn.tab', REFDIR + "/polyFiltered/{sample}.fasta
     output:REFDIR + '/nonBactMicoAgainstFungRefseq_blastn/{sample}.csv', REFDIR + '/nonBactMicoAgainstFungRefseq_blastn/{sample}_deduplicated.csv',REFDIR + '/nonBactMicoAgainstFungRefseq_blastn/{sample}_secondBlast.fasta'
     params:name='streamline_micotoFungRefseq',partition='quake,normal,owners',mem='200000',time='00:30:00'
     threads: 4
     conda:config['mainEnv_yaml']
     benchmark:REFDIR + "/benchmarks/streamline_micotoFungRefseq/{sample}.benchmark.txt"
     script: SCRIPTSDIR + "/dedup_v5.py"


#blastn of the refseq blastn results against the NT database
rule fungi_NT_blastn:
    input: r1 = REFDIR + '/nonBactMicoAgainstFungRefseq_blastn/{sample}_secondBlast.fasta'
    output: o1 = REFDIR + '/fungi_NT_blastn/{sample}.tab'
    params: name='fungi_NT_blastn',partition='quake,normal,owners',mem='80000',time='48:00:00'
    threads: 16
    benchmark:REFDIR + "/benchmarks/fungi_NT_blastn/{sample}.benchmark.txt"
    shell: """blastn -query {input.r1} -db {config[ntDatabaseDir]}  -num_threads {threads} -outfmt '6 qseqid sseqid stitle bitscore pident evalue gapopen qstart qend sstart send length mismatch staxids' -evalue 1e-5 -max_target_seqs 1 -out {output.o1}"""

#no second blast is required here, just want the csv of all blast results as opposed to individual .tab file per sample
rule streamline_fungi_NT_blastn:
    input: REFDIR + '/fungi_NT_blastn/{sample}.tab', REFDIR + "/polyFiltered/{sample}.fasta
    output:REFDIR + '/fungi_NT_blastn/{sample}.csv',REFDIR + '/fungi_NT_blastn/{sample}_deduplicated.csv',REFDIR + '/fungi_NT_blastn/{sample}_secondBlast.fasta'
    params:name='streamline_fungi_NT_blastn',partition='quake,normal,owners',mem='200000',time='00:30:00'
    threads:4
    conda:config['mainEnv_yaml']
    script: SCRIPTSDIR + "/dedup_v5.py"
