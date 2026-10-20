# OME Analysis Toolkit
A collection of scripts internal to the Ocean Molecular Ecology group at NOAA PMEL/UW CICOES for repeated data processing tasks.

## scripts\_sequencingPrep

### convert\_96wellPlate\_2\_longFormat.Rmd
This R markdown script imports the OME sequencing plate map file and converts to a long format for use downstream. For each OME sequencing run, one tab of the google sheet contains a clean copy of the sequencing plate maps. Though this script does seek out and identify each 96-well plate that it can detect in the file, the expected format is a plate identifier in the upper left hand corner (i.e. P#/P30 or MU) with 1-12 labeled columns proceeding to the right and A-H labeled rows proceeding down. Each well within the plate map should include the sample name. And immediately below the plate map prior to the next detected 96-well plate should be a cell indicating the marker in the format "OME P30 Sequencing MiFish". There is a user variable section in the Rmd where a marker\_map is created, which is where the user can align the words matching "MiFish" with the desired marker controlled vocab, i.e. "MiFish 12S". After changing the import/export information, the Rmd should be straightforward to run.


## scripts\_dataManagement

## scripts\_dataAnalysis

### 01\_MakePhyloseq\_ExploreData.Rmd
This R markdown script imports the OME sequencing data from the decontamination pipeline output for some data exploration and visualization using [phyloseq](https://joey711.github.io/phyloseq/index.html) and [microViz](https://david-barnett.github.io/microViz/). For each metabarcoding region of interest the script will read in the taxa, asv, metadata tables and fasta file, create a phyloseq object that is saved as an R data file (.rds), and produces a bargraph and ordination plot based on a metadata column of interest (e.g. sample collection device). To run the script you will need to ensure you've also downloaded the 'OME_analysis_utilities.R' file to load all libraries and functions. The path to this file will need to be edited in the first code chunk (`setup`). In the second (`load-seq-data`) and third (`make-phyloseq-obj`) code chunks you will need to add change to the filenames you are using (likely just the metarcoding region). Within the same code chunk ( ) and function (`create_phyloseq_object()`) you can select particular metadata columns of interest (`col_select`), specify and output directory (`output_dir`), and change the filename for the output (`filename`). In the next code chunk (`filer-phyloseq-metadata`), a column for pcr rep is added based on the sample naming schema, and all unknown taxa at the domain level are removed. Microviz will be used to create a bargraph and ordination plot which will be saved in a figures folder.

### OME\_analysis\_utilities.R
Includes functions to download and load packages, load data from decontam output, and to create a phyloseq object. More details coming soon...


### Legal Disclaimer

_This repository is a software product and is not official communication of the National Oceanic and Atmospheric Administration (NOAA), or the United States Department of Commerce (DOC). All NOAA GitHub project code is provided on an 'as is' basis and the user assumes responsibility for its use. Any claims against the DOC or DOC bureaus stemming from the use of this GitHub project will be governed by all applicable Federal law. Any reference to specific commercial products, processes, or services by service mark, trademark, manufacturer, or otherwise, does not constitute or imply their endorsement, recommendation, or favoring by the DOC. The DOC seal and logo, or the seal and logo of a DOC bureau, shall not be used in any manner to imply endorsement of any commercial product or activity by the DOC or the United States Government._
