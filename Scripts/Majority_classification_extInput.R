if (!suppressPackageStartupMessages(require("dplyr", quietly = TRUE)))
    install.packages("dplyr")
if (!suppressPackageStartupMessages(require("readr", quietly = TRUE)))
    install.packages("readr")
if (!suppressPackageStartupMessages(require("stringr", quietly = TRUE)))
    install.packages("stringr")
if (!suppressPackageStartupMessages(require("tidyr", quietly = TRUE)))
    install.packages("tidyr")

args = commandArgs(trailingOnly = TRUE)

contigs_classFile <- args[1]
whokaryote_classFile <- args[2]
tiara_classFile <- args[3]
deepmicrobefinder_classFile <- args[4]
kaiju_classFile <- args[5]
cat_classFile <- args[6]
out_folder <- args[7]
kmer_minLength<- args[8]
refdb_minLength <- args[9]
include_na <- args[10]
refdb_classifier <- args[11]
kaiju_local <- args[12]

results_filename <- sub("*.*/", "", contigs_classFile)
results_filename <- sub("_min.*", "", results_filename)

##### LOAD CONTIGS IDS #####

contigs_class <- read_tsv(contigs_classFile, col_names = FALSE, show_col_types = FALSE) %>% rename(ContigID = X1)
cat("Contigs IDs loaded\n")

##### LOAD CLASSIFICATION RESULTS #####

whokaryote_class <- read_delim(whokaryote_classFile, "\t", escape_double = FALSE, trim_ws = TRUE, show_col_types = FALSE) %>%
  rename(ContigID = 'contig', Whokaryote_class = predicted) %>%
  mutate(Whokaryote_class  = ifelse(Whokaryote_class  == "eukaryote", "EUK", "OTHER"))
cat("Whokaryote classification loaded\n")
tiara_class <- read_delim(tiara_classFile, "\t", escape_double = FALSE, trim_ws = TRUE, show_col_types = FALSE) %>%
  rename(ContigID = 'sequence_id', Tiara_class = class_fst_stage) %>%
  mutate(Tiara_class = ifelse(Tiara_class %in% c("eukarya","organelle"), "EUK", "OTHER"))
cat("Tiara classification loaded\n")
deepmicrobefinder_class <- read_delim(deepmicrobefinder_classFile, "\t", escape_double = FALSE, trim_ws = TRUE, show_col_types = FALSE) %>%
  rename(ContigID = 'Sequence Name', DeepMicrobeFinder_class = Prediction) %>%
  mutate(DeepMicrobeFinder_class = ifelse(DeepMicrobeFinder_class == "Eukaryote", "EUK", "OTHER"))
cat("Deepmicrobefinder classification loaded\n")

if (tolower(refdb_classifier) == "kaiju"){
	if (kaiju_local == "TRUE"){
	kaiju_class <- read_delim(kaiju_classFile, delim = "\t", col_names = FALSE, show_col_types = FALSE) %>%
		rename(ContigID = X2) %>%
		mutate(Kaiju_class = ifelse(str_detect(X8, "Eukaryota"), "EUK", "OTHER"))
    cat("Kaiju local classification loaded\n")
	} else{
		kaiju_class <- read_delim(kaiju_classFile, delim = "\t", escape_double = FALSE, col_names = TRUE, trim_ws = TRUE, show_col_types = FALSE) %>%
	  		rename(Kaiju_class = 'superkingdom') %>%
	  		mutate(Kaiju_class  = ifelse(Kaiju_class  == "Eukaryota", "EUK", "OTHER"))
        cat("Kaiju webserver classification loaded\n")
	}
} else if (tolower(refdb_classifier) == "cat"){
	cat_class <- read_delim(cat_classFile, delim = "\t", col_names = TRUE, show_col_types = FALSE) %>%
		rename(ContigID = `# contig`) %>%
		mutate(CAT_class = ifelse(str_detect(superkingdom, "Eukaryota"), "EUK", "OTHER")) 
    cat("CAT classification loaded\n")
}



##### KMER MAJORITY CLASSIFICATION #####

contigs_class <- deepmicrobefinder_class %>% select(ContigID, DeepMicrobeFinder_class) %>% right_join(contigs_class, by = "ContigID")
contigs_class <- tiara_class %>% select(ContigID, Tiara_class) %>% right_join(contigs_class, by = "ContigID")
contigs_class <- whokaryote_class %>% select(ContigID, Whokaryote_class) %>% right_join(contigs_class, by = "ContigID")

contigs_class <- contigs_class %>% mutate( 
  MajorityKmer_Nclassifications = (as.numeric(!is.na(DeepMicrobeFinder_class)) + as.numeric(!is.na(Tiara_class)) + as.numeric(!is.na(Whokaryote_class))),
  MajorityKmer_votesFrac = (str_detect(replace_na(DeepMicrobeFinder_class, ''),"EUK") + 
                                  str_detect(replace_na(Tiara_class, ''),"EUK") +  
                                  str_detect(replace_na(Whokaryote_class, ''),"EUK"))/MajorityKmer_Nclassifications,
  MajorityKmer_class = ifelse(MajorityKmer_votesFrac > 0.5, "EUK", "OTHER")) 


##### REFERENCE-BASED CLASSIFICATION #####

if (tolower(refdb_classifier) == "kaiju"){
	contigs_class <- kaiju_class %>% select(ContigID, Kaiju_class) %>% right_join(contigs_class, by = "ContigID")
	contigs_class <- contigs_class %>% 
	  mutate(MajorityKmer_Kaiju_class = ifelse(is.na(Kaiju_class), MajorityKmer_class, Kaiju_class))
	write_delim(contigs_class, paste0(out_folder, results_filename, ".Classification_details_KmerMajority_min",kmer_minLength,"bp_KAIJU_min",refdb_minLength,"bp.tsv"), delim = "\t")
	
	if (tolower(include_na) == "true"){
		EUKs_ContigID <- contigs_class %>% filter(MajorityKmer_Kaiju_class != "OTHER") %>% select(ContigID)
		write_delim(EUKs_ContigID, paste0(out_folder, results_filename, ".EUK_NA_contigsIDs_KmerMajority_min",kmer_minLength,"bp_KAIJU_min",refdb_minLength,"bp.txt"), col_names = FALSE)
	} else if (tolower(include_na) == "false") {
		EUKs_ContigID <- contigs_class %>% filter(MajorityKmer_Kaiju_class == "EUK") %>% select(ContigID)
		write_delim(EUKs_ContigID, paste0(out_folder, results_filename, ".EUK_contigsIDs_KmerMajority_min",kmer_minLength,"bp_KAIJU_min",refdb_minLength,"bp.txt"), col_names = FALSE)
	} else {
		cat("Specify if the output should include only the detected eukaryotes or also all the not classified contigs\n")
	}
} else if (tolower(refdb_classifier) == "cat"){
	contigs_class <- cat_class %>% select(ContigID, CAT_class) %>% right_join(contigs_class, by = "ContigID")
	contigs_class <- contigs_class %>% 
	  mutate(MajorityKmer_CAT_class = ifelse(is.na(CAT_class), MajorityKmer_class, CAT_class))
	write_delim(contigs_class, paste0(out_folder, results_filename, ".Classification_details_KmerMajority_min",kmer_minLength,"bp_CAT_min",refdb_minLength,"bp.tsv"), delim = "\t")
	
	if (tolower(include_na) == "true"){
		EUKs_ContigID <- contigs_class %>% filter(MajorityKmer_CAT_class != "OTHER") %>% select(ContigID)
		write_delim(EUKs_ContigID, paste0(out_folder, results_filename, ".EUK_NA_contigsIDs_KmerMajority_min",kmer_minLength,"bp_CAT_min",refdb_minLength,"bp.txt"), col_names = FALSE)
	} else if (tolower(include_na) == "false") {
		EUKs_ContigID <- contigs_class %>% filter(MajorityKmer_CAT_class == "EUK") %>% select(ContigID)
		write_delim(EUKs_ContigID, paste0(out_folder, results_filename, ".EUK_contigsIDs_KmerMajority_min",kmer_minLength,"bp_CAT_min",refdb_minLength,"bp.txt"), col_names = FALSE)
	} else {
		cat("Specify if the output should include only the detected eukaryotes or also all the not classified contigs\n")
	}
}
