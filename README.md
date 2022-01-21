# The hidden microbiome of human tissues and cell types
#### The human tissue microbiome atlas reveals numerous species of microbes spanning three domains of life
<img src="https://github.com/gitamahm/human_tissue_microbiome_atlas/blob/master/images/tissue_microbiome.png" width="100%" height="100%">

In the accompanying manuscript, we created an atlas of the human tissue microbiome with cell type resolution, and in this repository, we provide our pipeline called SIMBA (written in Snakemake) for the identification of bacteria, viruses and fungi in human cellular transcriptomes. We used the raw fastqs from the [Tabula Sapiens](https://github.com/czbiohub/tabula-sapiens) single-cell transcriptomic dataset. Additionally, we have included several Jupyter notebooks used for further processing SIMBA's outputs and generating figures. Theses notebooks include:

- **`post_SIMBA_processing_alldonors_and_validation.ipynb`**
    - adding taxonomic/lineage information
    - performing decontamination and QC steps
- **`single_cell_microbial_decontamination.ipynb`**
    - single-cell decontamination pipeline 
- **`tabula_sapiens_microbiome_analysis.ipynb`**
    - used for visualizing chordmaps 
    - upset plots 
    - microbe & human cell-type network 
    - EBV-human cell type sankey plot
- **`comparison_with_kraken2.ipynb`**
    - calculated the precision, recall and F1 of Kraken2 taxonomic assignments  
- **`ts_intersection_with_external_datasets.ipynb`**
    - integrating data from the Human Microbiome Project (HMP)
    - the tumor microbiome dataset by Nejman et al. (2020)
    - our validation dataset 
    - the PATRIC database  

