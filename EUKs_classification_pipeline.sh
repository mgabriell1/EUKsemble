#!/bin/bash

###### READ CONFIGURATION #####
source $1

let MINSIZE_CONTIGS_KMERCLASS_bp=MINSIZE_CONTIGS_KMERCLASS*1000
let MINSIZE_CONTIGS_KAIJUCLASS_bp=MINSIZE_CONTIGS_KAIJUCLASS*1000

######################################################################################
eval "$(conda shell.bash hook)"

###### PREPARE FILES #####
printf "\nClassfying assembly: $FOLDER$DATA$FILE_EXT \n\n"
mkdir -p $FOLDER$RESULTS_FOLDER

if [ ! -f $FOLDER$RESULTS_FOLDER$DATA"_min"$MINSIZE_CONTIGS_KMERCLASS"kbp"$FILE_EXT ] && [ ! -f $FOLDER$RESULTS_FOLDER$DATA"_min"$MINSIZE_CONTIGS_KAIJUCLASS"kbp"$FILE_EXT ]; then
    printf "Prepare files for classification \n"
    conda activate $SEQKIT_ENV
    seqkit seq $FOLDER$DATA$FILE_EXT -m $MINSIZE_CONTIGS_KMERCLASS_bp -o $FOLDER$RESULTS_FOLDER$DATA"_min"$MINSIZE_CONTIGS_KMERCLASS"kbp"$FILE_EXT 
    grep ">" $FOLDER$RESULTS_FOLDER$DATA"_min"$MINSIZE_CONTIGS_KMERCLASS"kbp"$FILE_EXT | awk 'sub(/^>/, "")'  > $FOLDER$RESULTS_FOLDER$DATA"_min"$MINSIZE_CONTIGS_KMERCLASS"kbp_contigsIDs.txt"

    if [ $MINSIZE_CONTIGS_KMERCLASS != $MINSIZE_CONTIGS_KAIJUCLASS ]; then
        seqkit seq $FOLDER$DATA$FILE_EXT -m $MINSIZE_CONTIGS_KAIJUCLASS_bp -o $FOLDER$RESULTS_FOLDER$DATA"_min"$MINSIZE_CONTIGS_KAIJUCLASS"kbp"$FILE_EXT
        grep ">" $FOLDER$RESULTS_FOLDER$DATA"_min"$MINSIZE_CONTIGS_KAIJUCLASS"kbp"$FILE_EXT | awk 'sub(/^>/, "")'  > $FOLDER$RESULTS_FOLDER$DATA"_min"$MINSIZE_CONTIGS_KAIJUCLASS"kbp_contigsIDs.txt"
    fi
    
	mkdir -p $KAIJURESULTS_FOLDER
	
    printf "Next steps for classification:\n 1. Perform Kaiju classification at https://kaiju.binf.ku.dk/server with the generated ${DATA}_min${MINSIZE_CONTIGS_KAIJUCLASS}kbp${FILE_EXT} file (after compression) and default settings\n 2. Save Kaiju's output file in $KAIJURESULTS_FOLDER \n 3. Re-run the pipeline script for k-mer classification and final output\n"
    exit
fi

printf '%*s\n' "${COLUMNS:-$(tput cols)}" '' | tr ' ' =

###### K-MER CLASSIFICATIONS #####

conda activate $WHOKARYOTE_ENV #Definition from whokaryote github

printf "Whokaryote classification \n"
if [ ! -d $FOLDER$RESULTS_FOLDER"whokaryote-results_min"$MINSIZE_CONTIGS_KMERCLASS"kbp" ]; then
    whokaryote.py --contigs $FOLDER$RESULTS_FOLDER$DATA"_min"$MINSIZE_CONTIGS_KMERCLASS"kbp"$FILE_EXT  \
        --outdir $FOLDER$RESULTS_FOLDER"whokaryote-results_min"$MINSIZE_CONTIGS_KMERCLASS"kbp" \
        --minsize $MINSIZE_CONTIGS_KMERCLASS_bp --f 
else
    printf "Whokaryote classification already present. Skipped \n"
fi
printf '%*s\n' "${COLUMNS:-$(tput cols)}" '' | tr ' ' =

conda activate $TIARA_ENV 

printf "Tiara classification \n"
if [ ! -f $FOLDER$RESULTS_FOLDER"tiara-results_min"$MINSIZE_CONTIGS_KMERCLASS"kbp/tiara-out_classification.txt" ]; then
    tiara -i $FOLDER$RESULTS_FOLDER$DATA"_min"$MINSIZE_CONTIGS_KMERCLASS"kbp"$FILE_EXT  \
        --to_fasta class all --threads 4 --probabilities --verbose \
        -p 0.65 0.65 \
        -m $MINSIZE_CONTIGS_KMERCLASS_bp \
        --output $FOLDER$RESULTS_FOLDER"tiara-results_min"$MINSIZE_CONTIGS_KMERCLASS"kbp/tiara-out_classification.txt"
else
    printf "Tiara classification already present. Skipped \n"
fi
    
printf '%*s\n' "${COLUMNS:-$(tput cols)}" '' | tr ' ' =

conda activate $DEEPMICROBEFINDER_ENV #Definition from DeepMicrobeFinder github

printf "DeepMicrobeFinder classification \n"
if [ ! -d $FOLDER$RESULTS_FOLDER"deepmicrobefinder-results_min"$MINSIZE_CONTIGS_KMERCLASS"kbp" ]; then
    python $DEEPMICROBEFINDER_FOLDER"predict.py" \
        -i $FOLDER$RESULTS_FOLDER$DATA"_min"$MINSIZE_CONTIGS_KMERCLASS"kbp"$FILE_EXT  \
        -e one-hot \
        -d $DEEPMICROBEFINDER_FOLDER"models/one-hot-models/" \
        -m hybrid \
        -o $FOLDER$RESULTS_FOLDER"deepmicrobefinder-results_min"$MINSIZE_CONTIGS_KMERCLASS"kbp"
        
    conda activate r-env

    $RSCRIPT_PATH Scripts/DeepMicrobeFinder_predictionAssignment_extInput.R $FOLDER$RESULTS_FOLDER"deepmicrobefinder-results_min"$MINSIZE_CONTIGS_KMERCLASS"kbp" $DATA"_min"$MINSIZE_CONTIGS_KMERCLASS"kbp"$FILE_EXT"_pred_one-hot_hybrid.txt"
else
    printf "DeepMicrobeFinder classification already present. Skipped \n"
fi

printf '%*s\n' "${COLUMNS:-$(tput cols)}" '' | tr ' ' =


##### ASSIGN TAXONOMY TO KAIJU WEBSERVER RESULTS #####

printf "Assign taxonomy from Kaiju webserver results \n"
$RSCRIPT_PATH Scripts/KaijuTaxaonomicAssignment_taxonomizr_extInput.R \
    $KAIJURESULTS_FOLDER \
    $KAIJURESULTS_FILENAME \
    $KAIJURESULTS_EXT \
    $MINSIZE_CONTIGS_KAIJUCLASS \
    $ACCESSIONTAXADB_TAXONOMIZR

printf '%*s\n' "${COLUMNS:-$(tput cols)}" '' | tr ' ' =

##### PERFORM CLASSIFICATION #####

printf "Perform majority classification \n"
$RSCRIPT_PATH Scripts/EUKs_majority_classification_extInput.R \
    $FOLDER$RESULTS_FOLDER$DATA"_min"$MINSIZE_CONTIGS_KAIJUCLASS"kbp_contigsIDs.txt" \
    $FOLDER$RESULTS_FOLDER"whokaryote-results_min"$MINSIZE_CONTIGS_KMERCLASS"kbp/featuretable_predictions_T.tsv" \
    $FOLDER$RESULTS_FOLDER"tiara-results_min"$MINSIZE_CONTIGS_KMERCLASS"kbp/tiara-out_classification.txt" \
    $FOLDER$RESULTS_FOLDER"deepmicrobefinder-results_min"$MINSIZE_CONTIGS_KMERCLASS"kbp/"$DATA"_min"$MINSIZE_CONTIGS_KMERCLASS"kbp"$FILE_EXT"_pred_one-hot_hybrid.txt" \
    $KAIJURESULTS_FOLDER$KAIJURESULTS_FILENAME"_min"$MINSIZE_CONTIGS_KAIJUCLASS"kbp_taxa"$KAIJURESULTS_EXT \
    $FOLDER$RESULTS_FOLDER \
    $MINSIZE_CONTIGS_KMERCLASS \
    $MINSIZE_CONTIGS_KAIJUCLASS

printf '%*s\n' "${COLUMNS:-$(tput cols)}" '' | tr ' ' =

##### EXTRACT EUKARYOTIC CONTIGS #####

printf "Extract eukaryotic contigs in new fasta file \n"
conda activate $SEQKIT_ENV
seqkit grep -f $FOLDER$RESULTS_FOLDER$DATA"_EUKs_contigsIDs_Kmer_min"$MINSIZE_CONTIGS_KMERCLASS"kbp_Kaiju_min"$MINSIZE_CONTIGS_KAIJUCLASS"kbp.txt" $FOLDER$RESULTS_FOLDER$DATA"_min"$MINSIZE_CONTIGS_KAIJUCLASS"kbp"$FILE_EXT -o $FOLDER$RESULTS_FOLDER$DATA"_Kmer_min"$MINSIZE_CONTIGS_KMERCLASS"kbp_Kaiju_min"$MINSIZE_CONTIGS_KAIJUCLASS"kbp_EUK"$FILE_EXT
