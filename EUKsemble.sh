#!/bin/bash

###### READ CONFIGURATION #####
source $1


######################################################################################
eval "$(conda shell.bash hook)"

###### PREPARE FILES #####
printf "\nClassifying assembly: $FOLDER$DATA$FILE_EXT \n\n"
mkdir -p $FOLDER"/"$RESULTS_FOLDER


if [ ${DB_CLASSIFIER,,} == "kaiju" ]; then
	mkdir -p $KAIJURESULTS_FOLDER
elif [ ${DB_CLASSIFIER,,} == "cat" ]; then
	mkdir -p $CATRESULTS_FOLDER
else
	printf "Reference-based classification method not recognized. \n"
	exit
fi
	

if [ ! -f $FOLDER"/"$RESULTS_FOLDER"/"$DATA"_min"$MINSIZE_CONTIGS_KMERCLASS"bp"$FILE_EXT ] && [ ! -f $FOLDER"/"$RESULTS_FOLDER"/"$DATA"_min"$MINSIZE_CONTIGS_DBCLASS"bp"$FILE_EXT ]; then
    printf "Prepare files for classification \n"
    conda activate $SEQKIT_ENV
    seqkit seq $FOLDER"/"$DATA$FILE_EXT -m $MINSIZE_CONTIGS_KMERCLASS -o $FOLDER"/"$RESULTS_FOLDER"/"$DATA"_min"$MINSIZE_CONTIGS_KMERCLASS"bp"$FILE_EXT 
    grep ">" $FOLDER"/"$RESULTS_FOLDER"/"$DATA"_min"$MINSIZE_CONTIGS_KMERCLASS"bp"$FILE_EXT | awk 'sub(/^>/, "")'  > $FOLDER"/"$RESULTS_FOLDER"/"$DATA"_min"$MINSIZE_CONTIGS_KMERCLASS"bp_contigsIDs.txt"

    if [ $MINSIZE_CONTIGS_KMERCLASS != $MINSIZE_CONTIGS_DBCLASS ]; then
        seqkit seq $FOLDER"/"$DATA$FILE_EXT -m $MINSIZE_CONTIGS_DBCLASS -o $FOLDER"/"$RESULTS_FOLDER"/"$DATA"_min"$MINSIZE_CONTIGS_DBCLASS"bp"$FILE_EXT
        grep ">" $FOLDER"/"$RESULTS_FOLDER"/"$DATA"_min"$MINSIZE_CONTIGS_DBCLASS"bp"$FILE_EXT | awk 'sub(/^>/, "")'  > $FOLDER"/"$RESULTS_FOLDER"/"$DATA"_min"$MINSIZE_CONTIGS_DBCLASS"bp_contigsIDs.txt"
    fi
	
	if [ $KAIJU_LOCAL != "TRUE" ]; then
	    printf "Next steps for classification:\n 1. Perform Kaiju classification at https://kaiju.binf.ku.dk/server with the generated ${DATA}_min${MINSIZE_CONTIGS_DBCLASS}bp${FILE_EXT} file (after compression) and default settings\n 2. Save Kaiju's output file in $KAIJURESULTS_FOLDER \n 3. Re-run the pipeline script for k-mer classification and final output\n"
	    exit
    fi
fi

printf "===============================================\n"

###### KAIJU CLASSIFICATIONS #####

if [ $KAIJU_LOCAL == "TRUE" ] && [ ${DB_CLASSIFIER,,} == "kaiju" ]; then
	printf "Kaiju classification \n"
	
	if [ ! -f $KAIJURESULTS_FOLDER"/"$KAIJURESULTS_FILENAME".taxa"$KAIJURESULTS_EXT ]; then
		conda activate $KAIJU_ENV
		kaiju -t $KAIJU_DB_NODES \
			-f $KAIJU_DB_FMI \
			-i $FOLDER"/"$RESULTS_FOLDER"/"$DATA"_min"$MINSIZE_CONTIGS_DBCLASS"bp"$FILE_EXT \
			-z $THREADS \
			-o $KAIJURESULTS_FOLDER"/"$KAIJURESULTS_FILENAME$KAIJURESULTS_EXT \
			-e $KAIJU_PARAMS_e -E $KAIJU_PARAMS_E -s $KAIJU_PARAMS_s -v
			
		kaiju-addTaxonNames -i $KAIJURESULTS_FOLDER"/"$KAIJURESULTS_FILENAME$KAIJURESULTS_EXT -o $KAIJURESULTS_FOLDER"/"$KAIJURESULTS_FILENAME".taxa"$KAIJURESULTS_EXT \
			-t $KAIJU_DB_NODES \
			-n $KAIJU_DB_NAMES -p 
		
	else
		printf "Kaiju classification already present. Skipped \n"
	fi
fi

###### CAT CLASSIFICATIONS #####

if [ ${DB_CLASSIFIER,,} == "cat" ]; then
	printf "CAT classification \n"
	
	conda activate $CAT_ENV
	if [ ! -f $CATRESULTS_FOLDER"/"$CATRESULTS_FILENAME".contig2classification.official_names.txt" ]; then
		$CAT_EXEC contigs -c $FOLDER"/"$RESULTS_FOLDER"/"$DATA"_min"$MINSIZE_CONTIGS_DBCLASS"bp"$FILE_EXT \
			-d $CAT_DB \
			-t $CAT_TAXONOMY \
			-n $THREADS \
			--index_chunks 1 --verbose $CAT_EXTRAPARAMS \
			--out_prefix $CATRESULTS_FOLDER"/"$CATRESULTS_FILENAME
			
		$CAT_EXEC add_names -i $CATRESULTS_FOLDER"/"$CATRESULTS_FILENAME".contig2classification.txt" \
			-t $CAT_TAXONOMY \
			--only_official \
			-o $CATRESULTS_FOLDER"/"$CATRESULTS_FILENAME".contig2classification.official_names.txt"
				
	else
		printf "CAT classification already present. Skipped \n"
	fi
fi

printf "===============================================\n"

###### K-MER CLASSIFICATIONS #####

conda activate $WHOKARYOTE_ENV #Definition from whokaryote github

printf "Whokaryote classification \n"
if [ ! -d $FOLDER"/"$RESULTS_FOLDER"/whokaryote-results_min"$MINSIZE_CONTIGS_KMERCLASS"bp" ]; then
    whokaryote.py --contigs $FOLDER"/"$RESULTS_FOLDER"/"$DATA"_min"$MINSIZE_CONTIGS_KMERCLASS"bp"$FILE_EXT  \
        --outdir $FOLDER"/"$RESULTS_FOLDER"/whokaryote-results_min"$MINSIZE_CONTIGS_KMERCLASS"bp" \
        --minsize $MINSIZE_CONTIGS_KMERCLASS --f 
else
    printf "Whokaryote classification already present. Skipped \n"
fi

printf "===============================================\n"

conda activate $TIARA_ENV 

printf "Tiara classification \n"
if [ ! -f $FOLDER"/"$RESULTS_FOLDER"/tiara-results_min"$MINSIZE_CONTIGS_KMERCLASS"bp/tiara-out_classification.txt" ]; then
    tiara -i $FOLDER"/"$RESULTS_FOLDER"/"$DATA"_min"$MINSIZE_CONTIGS_KMERCLASS"bp"$FILE_EXT  \
        --to_fasta class all --threads $THREADS --probabilities --verbose \
        -p $TIARA_FIRST_CLASS $TIARA_SECOND_CLASS \
        -m $MINSIZE_CONTIGS_KMERCLASS \
        --output $FOLDER"/"$RESULTS_FOLDER"/tiara-results_min"$MINSIZE_CONTIGS_KMERCLASS"bp/tiara-out_classification.txt"
else
    printf "Tiara classification already present. Skipped \n"
fi
    
printf "===============================================\n"

conda activate $DEEPMICROCLASS_ENV #Definition from DeepMicroClass github

printf "DeepMicroClass classification \n"
if [ ! -f $FOLDER"/"$RESULTS_FOLDER"/deepmicroclass-results_min"$MINSIZE_CONTIGS_KMERCLASS"bp/"$DATA"_min"$MINSIZE_CONTIGS_KMERCLASS"bp"$FILE_EXT"_pred_onehot_hybrid.class.tsv" ]; then
	if [ ${REMOTE_MODEL} == "TRUE" ]; then
	    DeepMicroClass predict \
	        -i $FOLDER"/"$RESULTS_FOLDER"/"$DATA"_min"$MINSIZE_CONTIGS_KMERCLASS"bp"$FILE_EXT  \
	        -e onehot \
	        -md hybrid \
	        -o $FOLDER"/"$RESULTS_FOLDER"/deepmicroclass-results_min"$MINSIZE_CONTIGS_KMERCLASS"bp"
	else 
		DeepMicroClass predict \
			-i $FOLDER"/"$RESULTS_FOLDER"/"$DATA"_min"$MINSIZE_CONTIGS_KMERCLASS"bp"$FILE_EXT  \
	        -e onehot \
	        -md hybrid \
	        -m $LOCAL_MODEL_LOC \
	        -o $FOLDER"/"$RESULTS_FOLDER"/deepmicroclass-results_min"$MINSIZE_CONTIGS_KMERCLASS"bp"
	fi
        
    $RSCRIPT_PATH $EUKSEMBLE_INSTALLATION_FOLDER"/Scripts/DeepMicroClass_predictionAssignment_extInput.R" $FOLDER"/"$RESULTS_FOLDER"/deepmicroclass-results_min"$MINSIZE_CONTIGS_KMERCLASS"bp" $DATA"_min"$MINSIZE_CONTIGS_KMERCLASS"bp"$FILE_EXT"_pred_onehot_hybrid.tsv"
else
    printf "DeepMicroClass classification already present. Skipped \n"
fi

printf "===============================================\n"

##### ASSIGN TAXONOMY TO KAIJU WEBSERVER RESULTS #####
if [ $KAIJU_LOCAL != "TRUE" ]; then
	printf "Assign taxonomy from Kaiju webserver results \n"
	$RSCRIPT_PATH $EUKSEMBLE_INSTALLATION_FOLDER"/Scripts/KaijuTaxaonomicAssignment_taxonomizr_extInput.R" \
	    $KAIJURESULTS_FOLDER"/" \
	    $KAIJURESULTS_FILENAME \
	    $KAIJURESULTS_EXT \
	    $MINSIZE_CONTIGS_DBCLASS \
	    $ACCESSIONTAXADB_TAXONOMIZR"/"
	
	printf "===============================================\n"
fi 
##### PERFORM CLASSIFICATION #####

printf "Perform majority classification \n"
$RSCRIPT_PATH $EUKSEMBLE_INSTALLATION_FOLDER"/Scripts/Majority_classification_extInput.R" \
    $FOLDER"/"$RESULTS_FOLDER"/"$DATA"_min"$MINSIZE_CONTIGS_DBCLASS"bp_contigsIDs.txt" \
    $FOLDER"/"$RESULTS_FOLDER"/whokaryote-results_min"$MINSIZE_CONTIGS_KMERCLASS"bp/featuretable_predictions_T.tsv" \
    $FOLDER"/"$RESULTS_FOLDER"/tiara-results_min"$MINSIZE_CONTIGS_KMERCLASS"bp/tiara-out_classification.txt" \
    $FOLDER"/"$RESULTS_FOLDER"/deepmicroclass-results_min"$MINSIZE_CONTIGS_KMERCLASS"bp/"$DATA"_min"$MINSIZE_CONTIGS_KMERCLASS"bp"$FILE_EXT"_pred_onehot_hybrid.class.tsv" \
    $KAIJURESULTS_FOLDER"/"$KAIJURESULTS_FILENAME".taxa"$KAIJURESULTS_EXT \
    $CATRESULTS_FOLDER"/"$CATRESULTS_FILENAME".contig2classification.official_names.txt"\
    $FOLDER"/"$RESULTS_FOLDER"/" \
    $MINSIZE_CONTIGS_KMERCLASS \
    $MINSIZE_CONTIGS_DBCLASS \
    $INCLUDE_NA \
    $DB_CLASSIFIER \
    $KAIJU_LOCAL

printf "===============================================\n"

##### EXTRACT EUKARYOTIC CONTIGS #####

printf "Extract eukaryotic contigs in new fasta file \n"
conda activate $SEQKIT_ENV
if [ $INCLUDE_NA == "TRUE" ]; then
	seqkit grep -f $FOLDER"/"$RESULTS_FOLDER"/"$DATA".EUK_NA_contigsIDs_KmerMajority_min"$MINSIZE_CONTIGS_KMERCLASS"bp_"${DB_CLASSIFIER^^}"_min"$MINSIZE_CONTIGS_DBCLASS"bp.txt" $FOLDER"/"$RESULTS_FOLDER"/"$DATA"_min"$MINSIZE_CONTIGS_DBCLASS"bp"$FILE_EXT -o $FOLDER"/"$RESULTS_FOLDER"/"$DATA".EUK_NA_KmerMajority_min"$MINSIZE_CONTIGS_KMERCLASS"bp_"${DB_CLASSIFIER^^}"_min"$MINSIZE_CONTIGS_DBCLASS"bp"$FILE_EXT
else
	seqkit grep -f $FOLDER"/"$RESULTS_FOLDER"/"$DATA".EUK_contigsIDs_KmerMajority_min"$MINSIZE_CONTIGS_KMERCLASS"bp_"${DB_CLASSIFIER^^}"_min"$MINSIZE_CONTIGS_DBCLASS"bp.txt" $FOLDER"/"$RESULTS_FOLDER"/"$DATA"_min"$MINSIZE_CONTIGS_DBCLASS"bp"$FILE_EXT -o $FOLDER"/"$RESULTS_FOLDER"/"$DATA".EUK_KmerMajority_min"$MINSIZE_CONTIGS_KMERCLASS"bp_"${DB_CLASSIFIER^^}"_min"$MINSIZE_CONTIGS_DBCLASS"bp"$FILE_EXT
fi
