library(xlsx)
library(dplyr)
library(knitr)
library(kableExtra)
library(rstudioapi)


# Set workdir to source file location
setwd(dirname(getActiveDocumentContext()$path))

# Load genera from classification spreadsheet
classification_df <- read.xlsx(file = "Table1.1_Classification_of_Rutaceae_genera.xlsx", sheetIndex = 1)

classification_df <- classification_df[with(classification_df, order(Appelhans_Subfamily, Kubitzki_Alliance, Kubitzki_Group,Genus)), ]

classification_df <- classification_df %>%
  select(Kubitzki_Alliance, Kubitzki_Group, Genus, Authority, species_number, species_num_reference, species_num_notes, Appelhans_Subfamily)

classification_df[is.na(classification_df)] <- " "

classification_df %>%
  kbl(row.names = F) %>%
  kable_styling() %>%
  pack_rows(index = table(classification_df$Appelhans_Subfamily))

# At this point table was copied from plot window and formatted in excel in the 'Classification_of_Rutaceae_genera.xlsx' document.






