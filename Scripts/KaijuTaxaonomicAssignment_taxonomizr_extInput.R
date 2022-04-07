if (!require("tidyverse", quietly = TRUE))
    install.packages("tidyverse")
if (!require("taxonomizr", quietly = TRUE))
    install.packages("taxonomizr")

args = commandArgs(trailingOnly = TRUE)

folder <- args[1]
filenames <- args[2]
ext <- args[3]
minContig_size <- args[4]
accessionTaxaDB_folder <- args[5]

prepareDatabase(paste0(accessionTaxaDB_folder,'accessionTaxa.sql')) 

for (filename in filenames){
  kaiju_out<-read.table(paste(folder,filename,ext,sep=""), sep = "\t",fill = TRUE, 
                        colClasses = c("character","character","numeric","numeric","character","character","character"),
                        col.names = c("Classified_Unclassified",	"ContigID",	"NCBI_TaxID",	"Score",	"TaxID_bestMatching" ,"Accession_bestMatching", "Seq_bestMatching"))
  cat(paste(filename, "loaded"))
  
  kaiju_out <- kaiju_out %>% mutate(row = row_number())
  
  kaiju_outC <- kaiju_out %>% filter(Classified_Unclassified == "C")
  kaiju_outC_TaxID <- kaiju_outC %>% pull(NCBI_TaxID)
  
  Kaiju_outC_taxa <-getTaxonomy(kaiju_outC_TaxID, paste0(accessionTaxaDB_folder,'accessionTaxa.sql'))
  
  kaiju_outC <- kaiju_outC %>% select(ContigID, row) %>% bind_cols(as_tibble(Kaiju_outC_taxa))
  
  kaiju_out_taxa <- left_join(kaiju_out, kaiju_outC, by = c("ContigID" = "ContigID", "row" = "row")) %>% select(-row)
  
  filename_out <- paste(folder,filename,"_min",minContig_size,"kbp_taxa",ext,sep="")
  write_tsv(kaiju_out_taxa,filename_out)
  cat(" and written with taxa \n")
  
}

