#!/bin/bash

###### READ CONFIGURATION #####
source $1


######################################################################################
eval "$(conda shell.bash hook)"

###### PREPARE FILES #####
printf "\nClassfying assembly: $FOLDER$DATA$FILE_EXT \n\n"
mkdir -p $FOLDER$RESULTS_FOLDER

if [ ! -f $FOLDER$RESULTS_FOLDER$DATA"_min"$MINSIZE_CONTIGS_KMERCLASS"bp"$FILE_EXT ] && [ ! -f $FOLDER$RESULTS_FOLDER$DATA"_min"$MINSIZE_CONTIGS_KAIJUCLASS"bp"$FILE_EXT ]; then
    printf "Prepare files for classification \n"
    conda activate $SEQKIT_ENV
    seqkit seq $FOLDER$DATA$FILE_EXT -m $MINSIZE_CONTIGS_KMERCLASS -o $FOLDER$RESULTS_FOLDER$DATA"_min"$MINSIZE_CONTIGS_KMERCLASS"bp"$FILE_EXT 
    grep ">" $FOLDER$RESULTS_FOLDER$DATA"_min"$MINSIZE_CONTIGS_KMERCLASS"bp"$FILE_EXT | awk 'sub(/^>/, "")'  > $FOLDER$RESULTS_FOLDER$DATA"_min"$MINSIZE_CONTIGS_KMERCLASS"bp_contigsIDs.txt"

    if [ $MINSIZE_CONTIGS_KMERCLASS != $MINSIZE_CONTIGS_KAIJUCLASS ]; then
        seqkit seq $FOLDER$DATA$FILE_EXT -m $MINSIZE_CONTIGS_KAIJUCLASS -o $FOLDER$RESULTS_FOLDER$DATA"_min"$MINSIZE_CONTIGS_KAIJUCLASS"bp"$FILE_EXT
        grep ">" $FOLDER$RESULTS_FOLDER$DATA"_min"$MINSIZE_CONTIGS_KAIJUCLASS"bp"$FILE_EXT | awk 'sub(/^>/, "")'  > $FOLDER$RESULTS_FOLDER$DATA"_min"$MINSIZE_CONTIGS_KAIJUCLASS"bp_contigsIDs.txt"
    fi
    
	mkdir -p $KAIJURESULTS_FOLDER
	
	if [ $KAIJU_LOCAL != "TRUE" ]; then
	    printf "Next steps for classification:\n 1. Perform Kaiju classification at https://kaiju.binf.ku.dk/server with the generated ${DATA}_min${MINSIZE_CONTIGS_KAIJUCLASS}bp${FILE_EXT} file (after compression) and default settings\n 2. Save Kaiju's output file in $KAIJURESULTS_FOLDER \n 3. Re-run the pipeline script for k-mer classification and final output\n"
	    exit
    fi
fi

%printf '%*s\n' "${COLUMNS:-$(tput cols)}" '' | tr ' ' =
printf "==============================================="

###### KAIJU CLASSIFICATIONS #####

if [ $KAIJU_LOCAL == "TRUE" ]; then
	printf "Kaiju classification \n"
	
	if [ ! -f $KAIJURESULTS_FOLDER$KAIJURESULTS_FILENAME".taxa"$KAIJURESULTS_EXT ]; then
		conda activate $KAIJU_ENV
		kaiju -t $KAIJU_DB_NODES \
			-f $KAIJU_DB_FMI \
			-i $FOLDER$RESULTS_FOLDER$DATA"_min"$MINSIZE_CONTIGS_KAIJUCLASS"bp"$FILE_EXT \
			-z $THREADS \
			-o $KAIJURESULTS_FOLDER$KAIJURESULTS_FILENAME$KAIJURESULTS_EXT \
			-e $KAIJU_PARAMS_e -E $KAIJU_PARAMS_E -s $KAIJU_PARAMS_s -v
			
		kaiju-addTaxonNames -i $KAIJURESULTS_FOLDER$KAIJURESULTS_FILENAME$KAIJURESULTS_EXT -o $KAIJURESULTS_FOLDER$KAIJURESULTS_FILENAME".taxa"$KAIJURESULTS_EXT \
			-t $KAIJU_DB_NODES \
			-n $KAIJU_DB_NAMES -p 
		
	else
		printf "Kaiju classification already present. Skipped \n"
	fi
fi

%printf '%*s\n' "${COLUMNS:-$(tput cols)}" '' | tr ' ' =
printf "==============================================="

###### K-MER CLASSIFICATIONS #####

source activate $WHOKARYOTE_ENV #Definition from whokaryote github

printf "Whokaryote classification \n"
if [ ! -d $FOLDER$RESULTS_FOLDER"whokaryote-results_min"$MINSIZE_CONTIGS_KMERCLASS"bp" ]; then
    whokaryote.py --contigs $FOLDER$RESULTS_FOLDER$DATA"_min"$MINSIZE_CONTIGS_KMERCLASS"bp"$FILE_EXT  \
        --outdir $FOLDER$RESULTS_FOLDER"whokaryote-results_min"$MINSIZE_CONTIGS_KMERCLASS"bp" \
        --minsize $MINSIZE_CONTIGS_KMERCLASS --f 
else
    printf "Whokaryote classification already present. Skipped \n"
fi
%printf '%*s\n' "${COLUMNS:-$(tput cols)}" '' | tr ' ' =
printf "==============================================="

conda activate $TIARA_ENV 

printf "Tiara classification \n"
if [ ! -f $FOLDER$RESULTS_FOLDER"tiara-results_min"$MINSIZE_CONTIGS_KMERCLASS"bp/tiara-out_classification.txt" ]; then
    tiara -i $FOLDER$RESULTS_FOLDER$DATA"_min"$MINSIZE_CONTIGS_KMERCLASS"bp"$FILE_EXT  \
        --to_fasta class all --threads $THREADS --probabilities --verbose \
        -p 0.65 0.65 \
        -m $MINSIZE_CONTIGS_KMERCLASS \
        --output $FOLDER$RESULTS_FOLDER"tiara-results_min"$MINSIZE_CONTIGS_KMERCLASS"bp/tiara-out_classification.txt"
else
    printf "Tiara classification already present. Skipped \n"
fi
    
%printf '%*s\n' "${COLUMNS:-$(tput cols)}" '' | tr ' ' =
printf "==============================================="

conda activate $DEEPMICROBEFINDER_ENV #Definition from DeepMicrobeFinder github

printf "DeepMicrobeFinder classification \n"
if [ ! -d $FOLDER$RESULTS_FOLDER"deepmicrobefinder-results_min"$MINSIZE_CONTIGS_KMERCLASS"bp" ]; then
    python $DEEPMICROBEFINDER_FOLDER"predict.py" \
        -i $FOLDER$RESULTS_FOLDER$DATA"_min"$MINSIZE_CONTIGS_KMERCLASS"bp"$FILE_EXT  \
        -e one-hot \
        -d $DEEPMICROBEFINDER_FOLDER"models/one-hot-models/" \
        -m hybrid \
        -o $FOLDER$RESULTS_FOLDER"deepmicrobefinder-results_min"$MINSIZE_CONTIGS_KMERCLASS"bp"
        
    $RSCRIPT_PATH Scripts/DeepMicrobeFinder_predictionAssignment_extInput.R $FOLDER$RESULTS_FOLDER"deepmicrobefinder-results_min"$MINSIZE_CONTIGS_KMERCLASS"bp" $DATA"_min"$MINSIZE_CONTIGS_KMERCLASS"bp"$FILE_EXT"_pred_one-hot_hybrid.txt"
else
    printf "DeepMicrobeFinder classification already present. Skipped \n"
fi

%printf '%*s\n' "${COLUMNS:-$(tput cols)}" '' | tr ' ' =
printf "==============================================="

##### ASSIGN TAXONOMY TO KAIJU WEBSERVER RESULTS #####
if [ $KAIJU_LOCAL != "TRUE" ]; then
	printf "Assign taxonomy from Kaiju webserver results \n"
	$RSCRIPT_PATH Scripts/KaijuTaxaonomicAssignment_taxonomizr_extInput.R \
	    $KAIJURESULTS_FOLDER \
	    $KAIJURESULTS_FILENAME \
	    $KAIJURESULTS_EXT \
	    $MINSIZE_CONTIGS_KAIJUCLASS \
	    $ACCESSIONTAXADB_TAXONOMIZR
	
	%printf '%*s\n' "${COLUMNS:-$(tput cols)}" '' | tr ' ' =
	printf "==============================================="
fi 
##### PERFORM CLASSIFICATION #####

printf "Perform majority classification \n"
$RSCRIPT_PATH Scripts/EUKs_majority_classification_extInput.R \
    $FOLDER$RESULTS_FOLDER$DATA"_min"$MINSIZE_CONTIGS_KAIJUCLASS"bp_contigsIDs.txt" \
    $FOLDER$RESULTS_FOLDER"whokaryote-results_min"$MINSIZE_CONTIGS_KMERCLASS"bp/featuretable_predictions_T.tsv" \
    $FOLDER$RESULTS_FOLDER"tiara-results_min"$MINSIZE_CONTIGS_KMERCLASS"bp/tiara-out_classification.txt" \
    $FOLDER$RESULTS_FOLDER"deepmicrobefinder-results_min"$MINSIZE_CONTIGS_KMERCLASS"bp/"$DATA"_min"$MINSIZE_CONTIGS_KMERCLASS"bp"$FILE_EXT"_pred_one-hot_hybrid.txt" \
    $KAIJURESULTS_FOLDER$KAIJURESULTS_FILENAME".taxa"$KAIJURESULTS_EXT \
    $FOLDER$RESULTS_FOLDER \
    $MINSIZE_CONTIGS_KMERCLASS \
    $MINSIZE_CONTIGS_KAIJUCLASS \
    $INCLUDE_NA \
    $KAIJU_LOCAL

%printf '%*s\n' "${COLUMNS:-$(tput cols)}" '' | tr ' ' =
printf "==============================================="

##### EXTRACT EUKARYOTIC CONTIGS #####

printf "Extract eukaryotic contigs in new fasta file \n"
conda activate $SEQKIT_ENV
if [ $INCLUDE_NA == "TRUE" ]; then
	seqkit grep -f $FOLDER$RESULTS_FOLDER$DATA"_EUKs_NAs_contigsIDs_Kmer_min"$MINSIZE_CONTIGS_KMERCLASS"bp_Kaiju_min"$MINSIZE_CONTIGS_KAIJUCLASS"bp.txt" $FOLDER$RESULTS_FOLDER$DATA"_min"$MINSIZE_CONTIGS_KAIJUCLASS"bp"$FILE_EXT -o $FOLDER$RESULTS_FOLDER$DATA"_Kmer_min"$MINSIZE_CONTIGS_KMERCLASS"bp_Kaiju_min"$MINSIZE_CONTIGS_KAIJUCLASS"bp_EUK_NA"$FILE_EXT
else
	seqkit grep -f $FOLDER$RESULTS_FOLDER$DATA"_EUKs_contigsIDs_Kmer_min"$MINSIZE_CONTIGS_KMERCLASS"bp_Kaiju_min"$MINSIZE_CONTIGS_KAIJUCLASS"bp.txt" $FOLDER$RESULTS_FOLDER$DATA"_min"$MINSIZE_CONTIGS_KAIJUCLASS"bp"$FILE_EXT -o $FOLDER$RESULTS_FOLDER$DATA"_Kmer_min"$MINSIZE_CONTIGS_KMERCLASS"bp_Kaiju_min"$MINSIZE_CONTIGS_KAIJUCLASS"bp_EUK"$FILE_EXT
fi
