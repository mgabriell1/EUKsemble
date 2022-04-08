# EUK_classification_pipeline
Scripts for identification of eukaryotic contigs in metagenomic assemblies

## Requirements
The following dependencies environments are required (the user can specify the name of a conda environment name for each tool):
- Seqkit (https://bioinf.shenwei.me/seqkit/)
- Tiara (https://github.com/ibe-uw/tiara/blob/main/docs/detailed-installation.md#detailed-installation) 
- Whokaryote (https://github.com/LottePronk/whokaryote#installation)
- DeepMicrobeFinder (https://github.com/chengsly/DeepMicrobeFinder#installation)
- R

The R scripts used will automatically install (if not present) the packages used:
  - Tidyverse
  - Taxonomizr

Beware that Taxonomizr will require to download a large database on its first usage

## Usage
1. Specify all variables including assembly name and location, minimum size required for Kaiju's and k-mer eukaryotic classification and the details regarding the location of the dependencies in `EUKs_classification_pipeline-CONFIG.conf`. \
Caveats:
    - Folder names have to end with `/`
    - Minimum contigs sizes are intended as kbp (e.g. `3 == 3000`)
3. Then execute `./EUKs_classification_pipeline.sh EUKs_classification_pipeline-CONFIG.conf`
