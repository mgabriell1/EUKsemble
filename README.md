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

## Important parameters and usage 
# Most imprtant parameters
To use `EUKsemble` these are the parameters to be changed: 
- Details on assembly to be classified:
	- `FOLDER`: folder contaning the contigs file
	- `DATA`: name of the contigs file
	- `FILE_EXT`: extention of the contigs file (e.g., `.fasta`)
- Minimum contigs length to be classified by k-mer-based majority (`MINSIZE_CONTIGS_KMERCLASS`) and Kaiju (`MINSIZE_CONTIGS_KAIJUCLASS`). The length is intended in bp (e.g., 1000 = 1 kbp). 
The minimum length can be different between the strategy to exploit the higher confidence of Kaiju's classification with short contigs. For this reason `MINSIZE_CONTIGS_KMERCLASS` needs to be equal or greater than `MINSIZE_CONTIGS_KAIJUCLASS`.
- Whether to include non classified contigs in the output fasta (`INCLUDE_NA`). This choice depends is left to the user as it depends on the specific workflow and the following amount of refinement. 
Our benchmark shows that while the vast majority of prokaryotic contigs at 1 kbp are classified by Kaiju, only about half of the eukaryotic ones are, meaning that the non-classified contigs might include a large fractino of eukaryotic contigs. 
- Paths to conda environments and binaries

# Usage
1. Specify all variables including assembly name and location, minimum size required for Kaiju's and k-mer eukaryotic classification and the details regarding the location of the dependencies in `EUKs_classification_pipeline-CONFIG.conf`. \
2. Then execute `./EUKs_classification_pipeline.sh EUKs_classification_pipeline-CONFIG.conf`

By default the results are saved in a subfolder next to the classified assembly called `{Assembly name}_EUK_classification/`, where `{Assembly name}` is the name of the assembly to be classified. This behaviour can be changed in the `OTHER PARAMETERS` section of the config file.
N.B.: In case the prediction is carried out on the webserver the pipeline has to be run twice: once, to prepare the files, and then after obtaining the contigs predictions provided by Kaiju's webserver. 

## Output
By default `EUKsemble` will create a folder alongside the file containing the classified contigs (can be changed editing `RESULTS_FOLDER`). This folder will contain:
- Four folders containing the results of the single tools (e.g., `deepmicrobefinder-results_min3000bp`, `kaiju-results_min1000bp`,`tiara-results_min3000bp`,`whokaryote-results_min3000bp`). The folder of k-mer-based tools might large as they include the contig files classified by each tools which might be removed if not necessary.
- Two or four Files containing the contigs longer than specified and the list of contigs IDs. The name of these files will start as specified in `DATA` followed by the minimum length of the included contigs (e.g., `DATA_min3000bp`)
- One tabular file with the classification results of all the tools on each contig named `DATA.Classification_details_Kmer_minXXXXbp_Kaiju_minYYYYbp.tsv` (`DATA` depends on the input provided as well as the minimum contig lengths). This file will contain 9 columns containing the contig ID, the classification by all the tools, the number of k-mer-based classification and the fraction of EUK votes and the majority voting result (either with or without including Kaiju).
For a normal usage the final result is included in the column `MajorityKmer_Kaiju_class` which indicates the classification based on Kaiju and k-mer-based majority voting.
- A text file containing the list of contigs selected as eukaryotic (including `NA` if specified). The file name is `DATA.EUK_contigsIDs_Kmer_minXXXXbp_Kaiju_minYYYYbp.txt` (in case `NA` are included the filename will include `EUK_NA`).
- A fasta file containg the contigs classified as eukaryotic (including `NA` if specified). The file name is `DATA.EUK_Kmer_minXXXXbp_Kaiju_minYYYYbp.txt` (in case `NA` are included the filename will include `EUK_NA`).


## Citation
Paper in progress. 
For now: "EUKsemble, an ensemble strategy for eukaryotic retrieval in metagenomes, Gabrielli and Pinto, 2022".
Also: do not forget to cite the software that made this pipeline possible (look at the requirements)

If you have found this pipeline useful let us know!

