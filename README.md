# The hidden microbiome of human tissues and cell types

<img src="https://github.com/gitamahm/human_tissue_microbiome_atlas/blob/master/images/tissue_microbiome.png" width="100%" height="100%">

The notebooks 
- **`post_SIMBA_processing_alldonors_and_validation.ipynb`**
    - adding taxonomic/lineage information
    - performing decontamination and QC based on alignment length and percent identity
    - merging cell barcodes with those found in the Ts objects to get cell type annotations
- **`single_cell_microbial_decontamination.ipynb`**
    - single-cell decontamination pipeline 
- **`tabula_sapiens_microbiome_analysis.ipynb`**
    - chordmaps (see below diagram)
    - upsetplots 
    - microbe & human cell-type network (see below diagram)
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
<img src="https://github.com/gitamahm/human_tissue_microbiome_atlas/blob/master/images/sankey.png" width="100%" height="100%">

## microbes shared between tissues
<img src="https://github.com/gitamahm/human_tissue_microbiome_atlas/blob/master/images/chordmap.png" width="50%" height="50%">

## microbes in association with human cell types
<img src="https://github.com/gitamahm/human_tissue_microbiome_atlas/blob/master/images/network.png" width="100%" height="100%">

