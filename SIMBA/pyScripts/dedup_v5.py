import pandas as pd
import os

from pandas import ExcelWriter
from pandas import ExcelFile
import numpy as np
import itertools
from collections import Counter

import glob
import re

import scipy as sc
from scipy import spatial
import itertools

from Bio import SeqIO
#from ete3 import NCBITaxa
#ncbi = NCBITaxa()
import phylopandas as ph

file=snakemake.input[0]
#converting .tab fasta results to .csv 
#will only select files that are not empty (i.e have blast output)
if os.stat(file).st_size != 0:
	df0=pd.read_table(file, delimiter='\t', names = ['seqName', 'refName', 'pathogen', 'bitscore', 'pident', 'evalue', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'length','mismatch', 'staxids'])
	df=df0
	#adding df to the filename, rather work with a dataframe with column headings from above
	df0.to_csv(snakemake.output[0], index=False)  
	
    	#getting the number of duplicated sequence names in the blast outputs (duplicated because they have the same cell barcode and umi)
	b = df.pivot_table(index=['seqName'], aggfunc='size').to_frame()
	c = b.rename(columns={0:'duplicates'})

    	#dropping the duplicates in the original dataframe
	d = df.drop_duplicates(subset='seqName', keep='first')
    	#merging the de-duplicated dataframe with another that contains information how many times a sequence was repeated.
	e = d.merge(c, on='seqName', how= 'inner')

    	#Reading in the original fasta file, quickly converting it to dataframe, merging that dataframe with blast output dataframe to identify seqs that need to be blasted against the larger databases
	df1 =ph.read_fasta(snakemake.input[1])

	df1 =df1.rename(columns={'id':'seqName', 'sequence': 'seq'})
	df1 =df1[['seqName', 'seq']]
    	#merging with the dataframe containing the deduplicated fasta output in addition to the query sequence which we need for other blast jobs
	final=df1.merge(e, how='inner', on='seqName')
	finaldf = final
    	#in addition to saving this csv file, we want to also write a fasta file containing the filtered query sequences (with blast output against the pathogeDB)
	final.to_csv(snakemake.output[1], index=False)

    	#creating the fasta file
	with open(snakemake.output[2], "w+") as output_file:
		for row in finaldf.index:
			seqName=finaldf.iloc[row,:]['seqName']
			seq = finaldf.iloc[row,:]['seq']
			output_file.write(">" + seqName + "\n" + seq + "\n")


#what to do with empty .tab blast results: we're going to create empty csv's so that snakemake doesn't fail going forward. 
else:
	pd.DataFrame({}).to_csv(snakemake.output[0], index=False)
	pd.DataFrame({}).to_csv(snakemake.output[1], index=False)
        #creating the fasta file
	with open(snakemake.output[2], "w+") as output_file:
		# this is a dummy fasta file so blast doesn't produce an error. It should result in an empty tab file. 
		output_file.write(">dummy_sequence" + "\n" + 'ACT' + "\n")



