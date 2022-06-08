# EUKsemble
Pipeline for identification of eukaryotic contigs in metagenomic assemblies. Best performances are obtained by using an ensemble approach of different k-mer classifiers and a reference-based one.  

This pipeline can run on a desktop computer if Kaiju's predictions are performed on the available webserver. Otherwise around 175GB of RAM are required for this step.

## Requirements
The following dependencies environments are required (the user can specify the name of a conda environment name for each tool):
- Seqkit (https://bioinf.shenwei.me/seqkit/)
- Tiara (https://github.com/ibe-uw/tiara/blob/main/docs/detailed-installation.md#detailed-installation) 
- Whokaryote (https://github.com/LottePronk/whokaryote#installation)
- DeepMicrobeFinder (https://github.com/chengsly/DeepMicrobeFinder#installation)
- Kaiju (https://github.com/bioinformatics-centre/kaiju, only if Kaiju predictions are performed locally)
- R

The R scripts used will automatically install (if not present) the packages used:
  - readr
  - dplyr
  - tidyr
  - Taxonomizr (only if Kaiju predictions are performed on the webserver)

Beware that Taxonomizr will require to download a large database on its first usage

## Usage
1. Specify all variables including assembly name and location, minimum size required for Kaiju's and k-mer eukaryotic classification and the details regarding the location of the dependencies in `EUKs_classification_pipeline-CONFIG.conf`. \
Caveats:
    - Folder names have to end with `/`
    - Minimum contigs sizes are intended as bp
    - Minimum contig length for k-mer-based classification must be equal or greater than Kaiju's one
2. Then execute `./EUKs_classification_pipeline.sh EUKs_classification_pipeline-CONFIG.conf`

By default the results are saved in a subfolder next to the classified assembly called `{Assembly name}_EUK_classification/`, where `{Assembly name}` is the name of the assembly to be classified. This behaviour can be changed in the `OTHER PARAMETERS` section of the config file.
N.B.: In case the prediction is carried out on the webserver the pipeline has to be run twice: once, to prepare the files, and then after obtaining the contigs predictions provided by Kaiju's webserver. 

## Citation
Paper in progress. 
For now: "EUKsemble, ensemble strategy for eukaryotic retrieval in metagenomes, Gabrielli and Pinto, 2022".
Also: do not forget to cite the software that made this pipeline possible (look at the requirements)

If you have found this pipeline useful let us know!

