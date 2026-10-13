# OME Analysis Toolkit
A collection of scripts internal to the Ocean Molecular Ecology group at NOAA PMEL/UW CICOES for repeated data processing tasks.

## scripts\_sequencingPrep

### convert\_96wellPlate\_2\_longFormat.Rmd
This R markdown script imports the OME sequencing plate map file and converts to a long format for use downstream. For each OME sequencing run, one tab of the google sheet contains a clean copy of the sequencing plate maps. Though this script does seek out and identify each 96-well plate that it can detect in the file, the expected format is a plate identifier in the upper left hand corner (i.e. P#/P30 or MU) with 1-12 labeled columns proceeding to the right and A-H labeled rows proceeding down. Each well within the plate map should include the sample name. And immediately below the plate map prior to the next detected 96-well plate should be a cell indicating the marker in the format "OME P30 Sequencing MiFish". There is a user variable section in the Rmd where a marker\_map is created, which is where the user can align the words matching "MiFish" with the desired marker controlled vocab, i.e. "MiFish 12S". After changing the import/export information, the Rmd should be straightforward to run.


## scripts\_dataManagement

## scripts\_dataAnalysis



### Legal Disclaimer

_This repository is a software product and is not official communication of the National Oceanic and Atmospheric Administration (NOAA), or the United States Department of Commerce (DOC). All NOAA GitHub project code is provided on an 'as is' basis and the user assumes responsibility for its use. Any claims against the DOC or DOC bureaus stemming from the use of this GitHub project will be governed by all applicable Federal law. Any reference to specific commercial products, processes, or services by service mark, trademark, manufacturer, or otherwise, does not constitute or imply their endorsement, recommendation, or favoring by the DOC. The DOC seal and logo, or the seal and logo of a DOC bureau, shall not be used in any manner to imply endorsement of any commercial product or activity by the DOC or the United States Government._
