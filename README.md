# The hidden microbiome of human tissues and cell types
#### The human tissue microbiome atlas reveals numerous species of microbes spanning three domains of life
<img src="https://github.com/gitamahm/human_tissue_microbiome_atlas/blob/master/images/tissue_microbiome.png" width="100%" height="100%">

SIMBA
A computational pipeline written in Snakemake for the identification of bacteria, viruses and fungi in single human cells. We used the raw fastqs from the [Tabula Sapiens](https://github.com/czbiohub/tabula-sapiens) single-cell transcriptomic dataset. 

Notebooks:
- **`post_SIMBA_processing_alldonors_and_validation.ipynb`**
    - adding taxonomic/lineage information
    - performing decontamination and QC based on alignment length and percent identity
    - merging cell barcodes with those found in the Ts objects to get cell type annotations
- **`single_cell_microbial_decontamination.ipynb`**
    - single-cell decontamination pipeline 
- **`tabula_sapiens_microbiome_analysis.ipynb`**
    - used for visualizing chordmaps 
    - upsetplots 
    - microbe & human cell-type network 
    - EBV-human cell type sankey plot
- **`comparison_with_kraken2.ipynb`**
    - calculated the precision, recall and F1 of Kraken2 taxonomic assignments  
- **`ts_intersection_with_external_datasets.ipynb`**
    - integrating data from the Human Microbiome Project (HMP)
    - the tumor microbiome dataset by Nejman et al. (2020)
    - our validation dataset consisting of bulk DNA sequencing of virus-like and bacterial-like particles derived from human tissues
    - the PATRIC database for pathogenicity and antimicrobial resistance profiles 
    - used for creating the following sankey plot of microbial flow routes


## microbial flow routes
#### Mapping the possible microbial flow routes from external-facing microbiomes to internal tissues and tumors

<img src="https://github.com/gitamahm/human_tissue_microbiome_atlas/blob/master/images/sankey.png" width="100%" height="100%">

## microbes shared between tissues
#### Human tissue microbiome tend to be unique to each donor and tissue
<img src="https://github.com/gitamahm/human_tissue_microbiome_atlas/blob/master/images/chordmap.png" width="50%" height="50%">

## microbes in association with human cell types
#### Network of human cell types in association with microbes reveals cell types with greatest microbial diversity 
<img src="https://github.com/gitamahm/human_tissue_microbiome_atlas/blob/master/images/network.png" width="100%" height="100%">

