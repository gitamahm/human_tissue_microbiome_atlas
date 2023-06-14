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
import phylopandas as ph


"""
filters a dataframe to eliminate sequences with polyA, polyC, polyG, or polT (as long as there are at least 16 of then in a row)
df_fil --> filtered dataset
"""
def filterPoly(df):
	df_fil = df[~df["seq"].str.contains('AAAAAAAAAAAAAAAA|GGGGGGGGGGGGGGGG|CCCCCCCCCCCCCCCC|TTTTTTTTTTTTTTTT')]
	return(df_fil)

"""
This is a function for writing fasta files from dataframe. 
It takes the dataframe, for each sequence, writes the cell barcode, and the sequence. 
"""
def makeFasta(df, filepath):
	with open(filepath,'w+') as output:
        	for index,row in enumerate(df.index):
            		output.write('>' + df.iloc[index]['id'] + '\n' + df.iloc[index]['seq'] + '\n')
	return()

#the fasta file that we want to filter
file=snakemake.input[0]
#the filtered fasta file we want to create, without the polys
outputfile=snakemake.output[0]

#will only select files that are not empty 
if os.stat(file).st_size != 0:
	df = ph.read_fasta(file)
	df.rename(columns={'sequence':'seq'}, inplace=True)
	df_fil =filterPoly(df)
	#what fraction of reads survived the poly filter
	surviving_fraction=df_fil.shape[0]/df.shape[0]
	print(surviving_fraction)	
    	#POST FILTERING: if no reads survive, we're writing a dummy fasta file so snakemake doesn't fail
	if df_fil.shape[0]==0:
        	with open(outputfile,'w+') as output_file:
                	output_file.write(">dummy_sequence" + "\n" + 'ACT' + "\n")
	else:
		makeFasta(df_fil, outputfile)
#what to do with empty files: we're going to create empty "filtered" files so that snakemake doesn't fail going forward. 
else:
	with open(outputfile, "w+") as output_file:
		output_file.write(">dummy_sequence" + "\n" + 'ACT' + "\n")


