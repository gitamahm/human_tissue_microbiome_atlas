# The hidden microbiome of human tissues and cell types


In `post_SIMBA_processing_alldonors_and_validation.ipynb` , I further processed SIMBA outputs. Some of the key steps include:
- adding taxonomic/lineage information to each microbial hit
- performing decontamination and QC based on alignment length and percent identity of each hit
- merging cell barcodes with those found in the Ts objects to get cell type annotations 

In `single_cell_microbial_decontamination.ipynb`, I performed additional microbial decontamination of single cell data. 

In `tabula_sapiens_microbiome_analysis.ipynb`, I take resulting high-confidence microbial hits within human tissues from curated in previous steps to create 
- chordmaps connecting microbes to human tissues (see below diagram)
- upsetplots 
- microbe human cell-type network (see below diagram)
- EBV-human cell type sankey plot

In `comparison_with_kraken2.ipynb`, I calculated the precision, recall and F1 of taxonomic assignments made by using Kraken2. 

Finally, in `ts_intersection_with_external_datasets.ipynb`, I integrated data from the 
- Human Microbiome Project (HMP)
- the tumor microbiome dataset by Nejman et al. (2020)
- our validation dataset consisting of bulk DNA sequencing of virus-like and bacterial-like particles derived from human tissues
- PATRIC database to get pathogenicity and antimicrobial resistance profiles for various bacterial species 

<img src="https://github.com/gitamahm/human_tissue_microbiome_atlas/blob/master/images/tissue_microbiome.png" width="100%" height="100%">

## microbial flow routes
<img src="https://github.com/gitamahm/human_tissue_microbiome_atlas/blob/master/images/sankey.png" width="100%" height="100%">

## microbes shared between tissues
<img src="https://github.com/gitamahm/human_tissue_microbiome_atlas/blob/master/images/chordmap.png" width="50%" height="50%">

## microbes in association with human cell types
<img src="https://github.com/gitamahm/human_tissue_microbiome_atlas/blob/master/images/network.png" width="100%" height="100%">

