if (!require("tidyverse", quietly = TRUE))
    install.packages("tidyverse")

args = commandArgs(trailingOnly = TRUE)

contigs_classFile <- args[1]
whokaryote_classFile <- args[2]
tiara_classFile <- args[3]
deepmicrobefinder_classFile <- args[4]
kaiju_classFile <- args[5]
out_folder <- args[6]
kmer_minLength<- args[7]
kaiju_minLength <- args[8]

results_filename <- sub("*.*/", "", contigs_classFile)
results_filename <- sub("_min.*", "", results_filename)

##### LOAD CONTIGS IDS #####

contigs_class <- read_csv(contigs_classFile, col_names = FALSE) %>% rename(ContigID = X1)


##### LOAD CLASSIFICATION RESULTS #####

whokaryote_class <- read_delim(whokaryote_classFile, "\t", escape_double = FALSE, trim_ws = TRUE) %>%
  rename(ContigID = 'contig', Whokaryote_class = predicted) %>%
  mutate(Whokaryote_class  = ifelse(Whokaryote_class  == "eukaryote", "EUK", "OTHER"))
tiara_class <- read_delim(tiara_classFile, "\t", escape_double = FALSE, trim_ws = TRUE) %>%
  rename(ContigID = 'sequence_id', Tiara_class = class_fst_stage) %>%
  mutate(Tiara_class = ifelse(Tiara_class == "eukarya", "EUK", "OTHER"))
deepmicrobefinder_class <- read_delim(deepmicrobefinder_classFile, "\t", escape_double = FALSE, trim_ws = TRUE) %>%
  rename(ContigID = 'Sequence Name', DeepMicrobeFinder_class = Prediction) %>%
  mutate(DeepMicrobeFinder_class = ifelse(DeepMicrobeFinder_class == "Eukaryote", "EUK", "OTHER"))

kaiju_class <- read_delim(kaiju_classFile, delim = "\t", escape_double = FALSE, col_names = TRUE, trim_ws = TRUE) %>%
  rename(Kaiju_class = 'superkingdom') %>%
  mutate(Kaiju_class  = ifelse(Kaiju_class  == "Eukaryota", "EUK", "OTHER"))


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


##### KAIJU CLASSIFICATION #####

contigs_class <- kaiju_class %>% select(ContigID, Kaiju_class) %>% right_join(contigs_class, by = "ContigID")
contigs_class <- contigs_class %>% 
  mutate(MajorityKmer_Kaiju_class = ifelse(is.na(Kaiju_class), MajorityKmer_class, Kaiju_class))
write_delim(contigs_class, paste0(out_folder, results_filename, "_EUKs_Classification_details_Kmer_min",kmer_minLength,"bp_Kaiju_min",kaiju_minLength,"bp.tsv"), delim = "\t")

EUKs_ContigID <- contigs_class %>% filter(MajorityKmer_Kaiju_class == "EUK") %>% select(ContigID)
write_delim(EUKs_ContigID, paste0(out_folder, results_filename, "_EUKs_contigsIDs_Kmer_min",kmer_minLength,"bp_Kaiju_min",kaiju_minLength,"bp.txt"), col_names = FALSE)
