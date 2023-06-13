# The hidden microbiome of human tissues and cell types
#### The human tissue microbiome atlas reveals numerous species of microbes spanning three domains of life
<img src="https://github.com/gitamahm/human_tissue_microbiome_atlas/blob/master/TSM.png" width="100%" height="100%">

In the accompanying manuscript, we created an atlas of the human tissue microbiome with cell type resolution, and in this repository, we provide our pipeline called SIMBA for the identification of bacteria, viruses and fungi in human cellular transcriptomes. We used the raw fastqs from the [Tabula Sapiens](https://github.com/czbiohub/tabula-sapiens) single-cell transcriptomic dataset. Additionally, we have included several Jupyter notebooks used for further processing SIMBA's outputs and generating figures. Theses notebooks include:

- **`p01_databaseCreationSC_microbeDB_fungalDB.ipynb`**
    - the purpose of this notebook is to create microbeDB and other databases used in the SIMBA pipeline
- **`p02_creating_truth_datasets_for_testing_simba.ipynb`**
	- the purpose of this notebook is to create truth datasets for testing the precision and recall of SIMBA and Kraken2
- **`p03_testing_truth_datasets.ipynb`**
	- In this notebook we estimate the precision and recall of SIMBA (as well as Kraken2) for detecting bacteria, fungi and viruses using various truth datasets created in the previous notebook. 
- **`p04_experimental_negative_controls.ipynb'**
	- In this notebook, we will assemble a list of contaminats from ~200 negative control samples from instruments and reagent used during tissue transport, digestion, and sequencing. These are from 5 10X lanes devoted to deeply sequencing negative controls
- **`p05_post_SIMBA_processing_alldonors_and_EHTM.ipynb`**
	- In this notebook, I perform additional processing of the SIMBA outputs for 10X and EHTM datasets. Some of the key steps include 
	- adding taxonomic/lineage information 
	- merging cell barcodes with those found in the Tabula Sapiens (TS) objects to get annotation data
	- formatting the TS objects such that cell type names are shortened and more neatly presented in subsequent figures
	- removed species found in negative controls from bulk tissue dataset. 
- **`p06_cell_ranger_counts_gh.ipynb`**
	- In this notebook we analyze the results of Cell Ranger that was run on TSP1-14 fastq files. Various processing steps are performed. 
- **`p07_machine_learning_modeling_contamination_and_other_filters.ipynb`**
	- This notebook delves into creating several columns (Filters 1-5). Filter 1-3 are contamination removal based on experimental negative controls, Salter et al genera, and fungi found in Narunsky-Haziza et al. Filter 4 is based on machine learning models of contamination using various inputs. Filter 5 involves filtering hits based on quality of alignments.
- **`p08_habitat_assignment_using_GEM_HMP_UHGG.ipynb`**
	- In this notebook we will use the HUE datasets comprising the Human Microbiome Project (HMP), Unified Human Gastrointestinal Genomes collection (UHGG), and the Genomes from Earth Microbiomes (GEM) datasets in order to identify putative habitats of origin for the taxa that we find at the intersection of these datasets and the Tabula Sapiens Microbiome (TSM) dataset. This notebook will be creating filter 6.  
- **`p09_davinci_language_model_habitat_of_origin.ipynb`**
	- This notebook is used to determine habitats of origin for 1500+ species not found in several large datasets (HUE datasets) using OpenAI's Davinci model. We will also benchmark this model against the same prompt for its self-consistency and its precision and recall against a truth dataset (100 random species selected from the GEM dataset). 
- **`p10_impact_of_filters.ipynb`**
	- This notebook is used to understand the impact of the 7 filters described in previous notebooks on the total number of hits, species, and genera.  
- **`p11_single_cell_microbiome_stats.ipynb`**
	- This notebook is used to explore stats surrounding the number of cells, empty droplets positive for hits from different filters. 
- **`p12_tabula_sapiens_multibiome_analysis.ipynb`**
	- This notebook is used to generate many of the plots found in the manuscript (e.g. chordmaps, upsetplots, networks, sankey) using already fully processed and QC'd SIMBA outputs. 
- **`p13_ts_intersection_with_external_datasets.ipynb`**
	- The purpose of this notebook is to compare the Tabula Sapiens Microbiome (TSM) dataset to the Human Microbiome Project (HMP), the tumor microbiome dataset by Nejman et al., and the EHTM dataset. Additionally, I obtain data from the PATRIC database on pathogenicity.


For a full list of libraries and various dependencies used in these notebooks, please see mainEnv2.yaml file, which can be used to create the right environment. See https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html for more information. 
