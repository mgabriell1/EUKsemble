if (!require("dplyr", quietly = TRUE))
    install.packages("dplyr")
if (!require("readr", quietly = TRUE))
    install.packages("readr")

args = commandArgs(trailingOnly = TRUE)

folder   <- args[1]
filename <- args[2]

setwd(folder)

contigs <- read_delim(filename, 
                      delim = "\t", escape_double = FALSE, 
                      trim_ws = TRUE)

contigs <- contigs %>% 
  mutate(Prediction = case_when(
    Eukaryote > EukaryoteVirus & Eukaryote > Plasmid & Eukaryote > Prokaryote & Eukaryote > ProkaryoteVirus ~ "Eukaryote",
    EukaryoteVirus > Eukaryote & EukaryoteVirus > Plasmid & EukaryoteVirus > Prokaryote & EukaryoteVirus > ProkaryoteVirus ~ "EukaryoteVirus",
    Plasmid > Eukaryote & Plasmid > EukaryoteVirus & Plasmid > Prokaryote & Plasmid > ProkaryoteVirus ~ "Plasmid",
    Prokaryote > Eukaryote & Prokaryote > EukaryoteVirus & Prokaryote > Plasmid & Prokaryote > ProkaryoteVirus ~ "Prokaryote",
    ProkaryoteVirus > Eukaryote & ProkaryoteVirus > EukaryoteVirus & ProkaryoteVirus > Plasmid & ProkaryoteVirus > Prokaryote ~ "ProkaryoteVirus"
  ))
#contigs %>% group_by(Prediction) %>% 
#  summarise(nContigs = n(), fracContigs = nContigs/nrow(contigs))

filename_class <- substr(filename,1,nchar(filename)-4)

write_delim(contigs, paste0(filename_class,".class.txt"), delim = "\t")
