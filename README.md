# The hidden microbiome of human tissues and cell types
#### The human tissue microbiome atlas reveals numerous species of microbes spanning three domains of life
<img src="https://github.com/gitamahm/human_tissue_microbiome_atlas/blob/master/images/tissue_microbiome_atlas.png" width="100%" height="100%">
<iframe src="https://github.com/gitamahm/human_tissue_microbiome_atlas/tree/master/images/f7_ann_host_microbe_network.html" width="100%" height="500px">

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
- **`ts_intersection_with_external_datasets.ipynb`**
    - integrating data from various sources (see `external_dataset` for input files)
    - the Human Microbiome Project (HMP)
    - the tumor microbiome dataset by Nejman et al. (2020)
    - our validation dataset 
    - the PATRIC database  

Below is the list of package requirements for analysis and visualization: 

    - scipy (v1.5.2)
    - pandas (v1.0.3)
    - numpy (v1.19.2)
    - matplotlib (v3.3.0)
    - seaborn (v0.11.1)
    - scanpy (v1.5.1)
    - anndata (v0.7.3)
    - phylopandas (v0.8.0)
    - holoviews (v1.14.0)
    - bokeh (v2.2.3)
    - venn (v0.1.3)
    - plotly (v4.14.1)
    - pyvis (v0.1.8.2)
    - networkx (v2.5)
    - upsetplot (v0.4.1) 
