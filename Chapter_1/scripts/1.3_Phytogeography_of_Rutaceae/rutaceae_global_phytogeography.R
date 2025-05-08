library(rstudioapi)
library(dplyr)
library(tidyverse)
library(progress)
library(kewr)
library(bazr)
library(sf)
library(ggplot2)
library(ape)
library(ggtree)
library(phytools)
library(geiger)
library(reshape2)
library(MonoPhy)
library(cowplot)
library(gridExtra)
library(scales)
library(TreeTools)
library(viridis)
library(igraph)
library(ggnetwork)
library(ips)
library(terra)
library(xlsx)
library(ggplotify)


# Set workdir to source file location
setwd(dirname(getActiveDocumentContext()$path))

# Load World Geographical Scheme for Recording Plant Distributions
data(tdwg_level1)
data(tdwg_level3)
st_crs(tdwg_level3) <- st_crs(tdwg_level1)
st_crs(tdwg_level3) <- st_crs(tdwg_level3)

# get rid of antarctic from lvl 1
tdwg_level1 <- rbind(tdwg_level1[1,], tdwg_level1[2:8,])

# Load GBIF collection records 
raw_gbif_data <- read.csv("01_GBIF_raw/occurrence.csv", header = T)


########################
### Data preparation ###
########################

# Can skip to visualization after this has been run once...

ifelse(!dir.exists(file.path("02_data_preparation")), dir.create(file.path("02_data_preparation")), FALSE)


# Trim dataset
slim_gbif_data <- raw_gbif_data[which(is.na(raw_gbif_data[1,]) == FALSE)]
slim_gbif_data <- slim_gbif_data[which(slim_gbif_data[1,] != "")]
slim_gbif_data <- slim_gbif_data %>%
  filter(taxonRank == "SPECIES", taxonomicStatus == "ACCEPTED", hasCoordinate == "TRUE")
# Get the species names in the retained dataset
taxon_names <- unique(paste(slim_gbif_data$genericName, slim_gbif_data$specificEpithet, sep = " "))

# Search POWO for Rutaceae species
Rutaceae_species <- search_powo("Rutaceae", filters=c("species", "accepted"), limit=3000)
Rutaceae_species_names <- purrr::map(Rutaceae_species$results, ~str_extract(.x$name, ".*"))
Rutaceae_species_ids <- purrr::map(Rutaceae_species$results, ~str_extract(.x$fqId, "[\\d\\-]+$"))

# Remove names from records that aren't in POWO
taxon_names[match(setdiff(taxon_names, Rutaceae_species_names), taxon_names)] <- NA
taxon_names <- na.omit(taxon_names)

# # Remove names from POWO lookup ids for optimisation
# Rutaceae_species_ids_df <- unlist(Rutaceae_species_ids)[match(intersect(unlist(Rutaceae_species_names), taxon_names), unlist(Rutaceae_species_names))]
# Rutaceae_species_names_df <- unlist(Rutaceae_species_names)[match(intersect(unlist(Rutaceae_species_names), taxon_names), unlist(Rutaceae_species_names))]
# 
# filtered_data <-data.frame("Species" = Rutaceae_species_names_df,
#                            "ID" = Rutaceae_species_ids_df)

# Function to do the lookups (to get species distributions)
fcn <- function(id) {
  lookup_powo(id, distribution=TRUE)
}

# Run the function
Rutaceae_lookups <- purrr::map(Rutaceae_species_ids, fcn)

## CODE SAVE CHECKPOINT ##
saveRDS(Rutaceae_lookups, file = "02_data_preparation/Rutaceae_lookups.rds")
saveRDS(slim_gbif_data, file = "02_data_preparation/slim_gbif_data.rds")


## CODE LOAD CHECKPOINT ##
Rutaceae_lookups <- readRDS(file = "02_data_preparation/Rutaceae_lookups.rds")
slim_gbif_data <- readRDS(file = "02_data_preparation/slim_gbif_data.rds")

# Convert the lookups to a dataframe
Rutaceae_checklist <- map_dfr(Rutaceae_lookups, tidy)

# Create a dictionary dataframe with the locations for all species
distribution_dictionary = data.frame(row.names = c("species", "native_location"))
for (i in 1:length(Rutaceae_checklist$modified)) {
  name <- Rutaceae_checklist[[25]][[i]]
  native_dist <- Rutaceae_checklist[[33]][[i]][[1]][[1]]$name
  distribution_dictionary <- rbind.data.frame(distribution_dictionary, data.frame("species"=name,"native_location"=native_dist))
}

# For loop to find the WGSRPD region from lat/long coords for every record and write it into a new column in the slim_gbif_data
point_in_polygon <- function(pt, ply){
  st_crs(pt) <- st_crs(ply)
  st_intersects(pt, st_sfc(ply))
}

for (i in 1:length(slim_gbif_data$gbifID)){
  # First insert the species name in same format as the distribution dictionary
  record_species_name <- paste(slim_gbif_data$genericName[i], slim_gbif_data$specificEpithet[i], sep = " ")
  slim_gbif_data$record_species_name[i] <- record_species_name
  
  # Extract the collection point
  point <- st_as_sf(data.frame("lat"=as.numeric(slim_gbif_data$decimalLatitude[i]),
                               "lon"=as.numeric(slim_gbif_data$decimalLongitude[i])),
                    coords = c("lon", "lat"), crs = 4326)
  
  # Get the row index for the polygon that contains the point
  row_index <- which(unlist(purrr::map(lapply(tdwg_level3$geometry, point_in_polygon, pt = point), ~str_extract(.x, ".*"))) == 1)
  
  # Write the region containing the point to "tdwg_lvl3_locality", if no match, then write NA
  if (length(row_index) == 1){
    slim_gbif_data$tdwg_lvl3_locality[i] <- as.character(tdwg_level3$LEVEL3_NAM[row_index])
  } else {
    slim_gbif_data$tdwg_lvl3_locality[i] <- NA
  }
  print(paste0("Progress: Record ", i, "/", length(slim_gbif_data$gbifID)))
}

# Save/load the raw output
saveRDS(distribution_dictionary, file = "02_data_preparation/distribution_dictionary.rds")
saveRDS(slim_gbif_data, file = "02_data_preparation/slim_gbif_data_distributions_raw.rds")

distribution_dictionary <- readRDS(file = "02_data_preparation/distribution_dictionary.rds")
slim_gbif_data_distributions_raw <- readRDS(file = "02_data_preparation/slim_gbif_data_distributions_raw.rds")

# Filter out records whose coords did not fall in a region
slim_gbif_data_regions_filtered <- slim_gbif_data_distributions_raw %>%
  filter(!is.na(tdwg_lvl3_locality))

# Filter out records whose 'tdwg_lvl3_locality' doesn't match the native range for the species in the distribution dictionary
idxs_2_keep <- c()
for (i in 1:length(slim_gbif_data_regions_filtered$record_species_name)){
  dist_i <- distribution_dictionary$native_location[which(distribution_dictionary$species %in% slim_gbif_data_regions_filtered$record_species_name[i])]
  if (slim_gbif_data_regions_filtered$tdwg_lvl3_locality[i] %in% dist_i){
    idxs_2_keep <- c(idxs_2_keep, i)
  }
  print(i)
}
slim_gbif_data_regions_filtered <- slim_gbif_data_regions_filtered[idxs_2_keep,]

# Fetch the tdgw lvl 1 locality for each record
slim_gbif_data_regions_filtered$tdwg_lvl1_locality <- tdwg_level1$LEVEL1_NAM[tdwg_level3$LEVEL1_COD[match(slim_gbif_data_regions_filtered$tdwg_lvl3_locality, tdwg_level3$LEVEL3_NAM)]]
distribution_dictionary$tdwg_lvl1_location <- tdwg_level1$LEVEL1_NAM[tdwg_level3$LEVEL1_COD[match(distribution_dictionary$native_location, tdwg_level3$LEVEL3_NAM)]]

## Some final data cleaning
# Convert Geijera paniculata to Coatesia paniculata
a <- slim_gbif_data_regions_filtered[grep("Geijera paniculata", slim_gbif_data_regions_filtered$scientificName), ]
a[] <- lapply(a, gsub, pattern = "Geijera", replacement = "Coatesia")
slim_gbif_data_regions_filtered[grep("Geijera paniculata", slim_gbif_data_regions_filtered$scientificName), ] <- a
# Remove Feroniella
slim_gbif_data_regions_filtered <- slim_gbif_data_regions_filtered %>% filter(genus != "Feroniella")
# Add artificial rows for Desmotes, Limnocitrus and Merrillia
a <- as.data.frame(matrix(data = c(c(rep(NA,72), "Desmotes incomparabilis", "Panamá", "SOUTHERN AMERICA"),
                                   c(rep(NA,72), "Desmotes incomparabilis", "Colombia", "SOUTHERN AMERICA"),
                                   c(rep(NA,72), "Limnocitrus littoralis", "Jawa", "ASIA-TROPICAL"),
                                   c(rep(NA,72), "Limnocitrus littoralis", "Vietnam", "ASIA-TROPICAL"),
                                   c(rep(NA,72), "Merrillia caloxylon", "Malaya", "ASIA-TROPICAL"),
                                   c(rep(NA,72), "Merrillia caloxylon", "Thailand", "ASIA-TROPICAL")), nrow = 6, ncol = 75, byrow = T))
colnames(a) <- colnames(slim_gbif_data_regions_filtered)
slim_gbif_data_regions_filtered <- rbind(slim_gbif_data_regions_filtered, a)

# Save filtered data
saveRDS(slim_gbif_data_regions_filtered, file = "02_data_preparation/slim_gbif_data_regions_filtered.rds")
write.table(slim_gbif_data_regions_filtered, file = "02_data_preparation/slim_gbif_data_regions_filtered.txt", quote = FALSE, row.names = FALSE, sep = "\t")
#saveRDS(distribution_dictionary, file = "02_data_preparation/distribution_dictionary.rds")

###################################
### Data prep for visualisation ###
###################################

### Load input data ###

# Create a directory for the plots and output datafiles
ifelse(!dir.exists(file.path("04_raw_plots")), dir.create(file.path("04_raw_plots")), FALSE)
ifelse(!dir.exists(file.path("05_data_files")), dir.create(file.path("05_data_files")), FALSE)

# Load GBIF data with distributions
slim_gbif_data_regions_filtered <- readRDS(file = "02_data_preparation/slim_gbif_data_regions_filtered.rds")

# Load the POWO distribution dictionary
distribution_dictionary <- readRDS(file = "02_data_preparation/distribution_dictionary.rds")
colnames(distribution_dictionary) <- c('record_species_name', 'tdwg_lvl3_locality', 'tdwg_lvl1_locality')

# Load World Geographical Scheme for Recording Plant Distributions
data(tdwg_level1)
data(tdwg_level3)
st_crs(tdwg_level3) <- st_crs(tdwg_level1)
st_crs(tdwg_level3) <- st_crs(tdwg_level3)
# ... get rid of antarctic from lvl 1
tdwg_level1 <- rbind(tdwg_level1[1,], tdwg_level1[2:8,])



### Create the tree to map to ###

# Load Joyce 2023 IQTREE phylogeny Rutaceae subtree...
rut_tree <- read.tree(file = "03_Joyce2023_trees/Supp7_IQTREE_230119_Rutaceae_subtree.tre")

# Collapse unsupported branches
rut_tree <- collapseUnsupportedEdges(rut_tree, value = "node.label", 100)

# ... clean up some samples/taxa and collapse monophyletic genera with > 1 sample...
rut_tree$tip.label <- gsub("Murraya_koenigii", "Bergera_koenigii", rut_tree$tip.label)
rut_tree <- drop.tip(rut_tree, "Melicope_broadbentiana")
rut_tree$tip.label <- gsub("_.*|'", "", rut_tree$tip.label)
rut_tree <- drop.tip(rut_tree, "Boronella")
rut_tree$tip.label <- gsub("Chloroxylum", "Chloroxylon", rut_tree$tip.label)
rut_tree$tip.label <- gsub("Rhadinothamnus", "Chorilaena", rut_tree$tip.label)
sol <- AssessMonophyly(rut_tree, taxonomy = data.frame("tips" = rut_tree$tip.label, "genus" = rut_tree$tip.label))
rut_tree <- CollapseMonophyletics(solution = sol, rut_tree, taxlevels = 1)

# ... manually add genera not included in the Joyce tree based on positions in other studies - * denotes position in another phylogenetic study, ** denotes genera never including in a phylogeny and positioned based on taxonomic hypotheses
rut_tree <- AddTip(rut_tree, label = "Acradenia*", where = getMRCA(rut_tree, c("Bosistoa", "Bouchardatia")), edgeLength = 0.1) # From Appelhans et al. 2021
rut_tree <- AddTip(rut_tree, label = "Crossosperma*", where = "Acradenia*", edgeLength = 0.1) # From Appelhans et al. 2021
rut_tree <- AddTip(rut_tree, label = "Pitavia*", where = "Acradenia*", edgeLength = 0.1) # From Appelhans et al. 2021
rut_tree <- AddTip(rut_tree, label = "Maclurodendron*", where = "Acronychia", edgeLength = 0.1) # From Appelhans et al. 2021
rut_tree <- AddTip(rut_tree, label = "Cyanothamnus*", where = getMRCA(rut_tree, c("Tetractomia", "Dutaillyea")), edgeLength = 0.1) # From Appelhans et al. 2021
rut_tree <- AddTip(rut_tree, label = "Perryodendron*", where = getMRCA(rut_tree, c("Zieria", "Pitaviaster")), edgeLength = 0.1, lengthBelow = 0) # From Appelhans et al. 2021
rut_tree <- AddTip(rut_tree, label = "Eriostemon*", where = getMRCA(rut_tree, c("Crowea", "Leionema")), edgeLength = 0.1, lengthBelow = 0) # From Appelhans et al. 2021
rut_tree <- AddTip(rut_tree, label = "Diplolaena*", where = getMRCA(rut_tree, c("Nematolepis", "Leionema")), edgeLength = 0.1, lengthBelow = 0) # From Appelhans et al. 2021
rut_tree <- AddTip(rut_tree, label = "Neoschmidia*", where = getMRCA(rut_tree, c("Muiriantha", "Leionema")), edgeLength = 0.1) # From Appelhans et al. 2021
rut_tree <- AddTip(rut_tree, label = "Sigmatanthus*", where = "Erythrochiton", edgeLength = 0.1) # From Appelhans et al. 2021
rut_tree <- AddTip(rut_tree, label = "Dryades*", where = getMRCA(rut_tree, c("Neoraputia", "Galipea")), edgeLength = 0.1, lengthBelow = 0) # From Appelhans et al. 2021
rut_tree <- AddTip(rut_tree, label = "Desmotes*", where = "Toxosiphon", edgeLength = 0.1) # From Appelhans et al. 2021
rut_tree <- AddTip(rut_tree, label = "Esenbeckia*", where = "Helietta", edgeLength = 0.1) # From Appelhans et al. 2021
rut_tree <- AddTip(rut_tree, label = "Choisya*", where = "Peltostigma", edgeLength = 0.1) # From Appelhans et al. 2021
rut_tree <- AddTip(rut_tree, label = "Naringi*", where = "Citropsis", edgeLength = 0.1) # From Appelhans et al. 2021
rut_tree <- AddTip(rut_tree, label = "Wenzelia*", where = getMRCA(rut_tree, c("Merope", "Monanthocitrus")), edgeLength = 0.1) # From Appelhans et al. 2021
rut_tree <- AddTip(rut_tree, label = "Pamburus*", where = "Luvunga", edgeLength = 0.1) # From Appelhans et al. 2021
rut_tree <- AddTip(rut_tree, label = "Paramignya*", where = "Pamburus*", edgeLength = 0.1) # From Appelhans et al. 2021
rut_tree <- AddTip(rut_tree, label = "Micromelum*", where = "Bergera", edgeLength = 0.1) # From Appelhans et al. 2021
rut_tree <- AddTip(rut_tree, label = "Dutailliopsis**", where = "Dutaillyea", edgeLength = 0.1) # From Hartley 1995
rut_tree <- AddTip(rut_tree, label = "Polyaster**", where = getMRCA(rut_tree, c("Decazyx", "Peltostigma")), edgeLength = 0.1, lengthBelow = 0) # From Kubitzki et al. 2011
rut_tree <- AddTip(rut_tree, label = "Apocaulon**", where = getMRCA(rut_tree, c("Decagonocarpus", "Ertela")), edgeLength = 0.1, lengthBelow = 0) # From Kubitzki et al. 2011
rut_tree <- AddTip(rut_tree, label = "Euxylophora**", where = "Hortia", edgeLength = 0.1) # From Kubitzki et al. 2011
rut_tree <- AddTip(rut_tree, label = "Leptothyrsa**", where = "Adiscanthus", edgeLength = 0.1) # From Kubitzki et al. 2011
rut_tree <- AddTip(rut_tree, label = "Megastigma**", where = getMRCA(rut_tree, c("Decazyx", "Peltostigma")), edgeLength = 0.1, lengthBelow = 0) # From Kubitzki et al. 2011
rut_tree <- AddTip(rut_tree, label = "Naudinia**", where = getMRCA(rut_tree, c("Hortia", "Galipea")), edgeLength = 0.1, lengthBelow = 0) # From Kubitzki et al. 2011
rut_tree <- AddTip(rut_tree, label = "Raulinoa**", where = getMRCA(rut_tree, c("Metrodorea", "Balfourodendron")), edgeLength = 0.1, lengthBelow = 0) # From Kubitzki et al. 2011
rut_tree <- AddTip(rut_tree, label = "Rutaneblina**", where = getMRCA(rut_tree, c("Naudinia**", "Neoraputia")), edgeLength = 0.1, lengthBelow = 0) # From Kubitzki et al. 2011

# Convert the zero-length branches to real polytomies
rut_tree <-di2multi(rut_tree, tol = 0.00001)

# Rotate some sister branches for nicer plotting
rut_tree <- ape::rotate(rut_tree, node = c("Casimiroa", "Stauranthus"), polytom = c(2, 1))

# Test plot and save the tree
plot(rut_tree, cex = 0.7)
#write.tree(rut_tree, "02_data_preparation/Rutaceae_summary_tree_final.tre")
rut_tree <- ape::read.tree("02_data_preparation/Rutaceae_summary_tree_final.tre")


### Create different datasets for different plots ###

## Taxon distributions for TDWG levels 1 and 3, at ranks of genus and species
# Here I have made a best attempt to get as close to accurate as possible given the flaws in both GBIF and POWO

# How this works - Our very conservative data has been verified by comparing GBIF collections 
# against species distributions in POWO (i.e. it only contains species and distributions where a herb specimen
# has been collected from a location that is within the area that POWO says the species occurs in). The 
# downside of this is that many species from many locations are not in the GBIF data so they do not get 
# included. Here, we go back to the POWO distribution dictionary, select the species from genera which are
# have less species included there should be, and add the dict data for those extra species. We also add the dict
# data for species that are already in the dataset, in case they are only represented in a subset of the locations
# that they actually occur in. 

# Species at level 1
# Subset the data into unique combinations of "record_species_name" and "tdwg_lvl1_locality"
species_at_lvl1 <- unique(slim_gbif_data_regions_filtered[, c('record_species_name','tdwg_lvl1_locality')])
species_at_lvl1_dd <- unique(distribution_dictionary[, c('record_species_name','tdwg_lvl1_locality')])
a<-setdiff(species_at_lvl1_dd$record_species_name, species_at_lvl1$record_species_name)
powo_better_genera <- "Geleznowia|Phebalium|Cyanothamnus|Chorilaena|Zieria|Adenandra|Aegle|Aeglopsis|Agathosma|Bergera|Boronia|Choisya|Citrus|Correa|Coleonema|Conchocarpus|Drummondita|Dryades|Empleurum|Flindersia|Glycosmis|Haplophyllum|Ivodea|Leionema|Lubaria|Luvunga|Maclurodendron|Medicosma|Melicope|Micromelum|Murraya|Paramignya|Phebalium|Philotheca|Pilocarpus|Wenzelia|Vepris|Zanthoxylum"
to_add <- species_at_lvl1_dd[which(species_at_lvl1_dd$record_species_name %in% a[grep(powo_better_genera, a)]),] %>% filter(!grepl("×", record_species_name))
to_add2 <-species_at_lvl1_dd[which(species_at_lvl1_dd$record_species_name %in% species_at_lvl1$record_species_name),]

species_at_lvl1 <- rbind(species_at_lvl1, to_add, to_add2)
species_at_lvl1 <- unique(species_at_lvl1[, c('record_species_name','tdwg_lvl1_locality')]) %>% filter(!is.na(tdwg_lvl1_locality))

# Genera at level 1
genera_at_lvl1 <- species_at_lvl1
genera_at_lvl1$record_species_name <- gsub(" .*", "", genera_at_lvl1$record_species_name)
colnames(genera_at_lvl1)[1] <- 'record_genus_names'
genera_at_lvl1 <- unique(genera_at_lvl1[, c('record_genus_names','tdwg_lvl1_locality')])
# Some manual edits
genera_at_lvl1 <- rbind(genera_at_lvl1, c("Halfordia", "AUSTRALASIA"))

# Species at level 3 
# Subset the data into unique combinations of "record_species_name" and "tdwg_lvl3_locality"
species_at_lvl3 <- unique(slim_gbif_data_regions_filtered[, c('record_species_name','tdwg_lvl3_locality')])
species_at_lvl3_dd <- unique(distribution_dictionary[, c('record_species_name','tdwg_lvl3_locality')])
a<-setdiff(species_at_lvl3_dd$record_species_name, species_at_lvl3$record_species_name)
powo_better_genera <- "Geleznowia|Phebalium|Cyanothamnus|Chorilaena|Zieria|Adenandra|Aegle|Aeglopsis|Agathosma|Bergera|Boronia|Choisya|Citrus|Correa|Coleonema|Conchocarpus|Drummondita|Dryades|Empleurum|Flindersia|Glycosmis|Haplophyllum|Ivodea|Leionema|Lubaria|Luvunga|Maclurodendron|Medicosma|Melicope|Micromelum|Murraya|Paramignya|Phebalium|Philotheca|Pilocarpus|Wenzelia|Vepris|Zanthoxylum"
to_add <- species_at_lvl3_dd[which(species_at_lvl3_dd$record_species_name %in% a[grep(powo_better_genera, a)]),] %>% filter(!grepl("×", record_species_name))
to_add2 <-species_at_lvl3_dd[which(species_at_lvl3_dd$record_species_name %in% species_at_lvl3$record_species_name),]

species_at_lvl3 <- rbind(species_at_lvl3, to_add, to_add2)
species_at_lvl3 <- unique(species_at_lvl3[, c('record_species_name','tdwg_lvl3_locality')]) %>% filter(!is.na(tdwg_lvl3_locality))
species_at_lvl3$tdwg_lvl3_locality[species_at_lvl3$tdwg_lvl3_locality == "Panamá"] <- "Panama"
species_at_lvl3$tdwg_lvl3_locality[species_at_lvl3$tdwg_lvl3_locality == "Suriname"] <- "Surinam"
species_at_lvl3$tdwg_lvl3_locality[species_at_lvl3$tdwg_lvl3_locality == "Leeward Is."] <- "Leeward Is. AB Ant"
species_at_lvl3$tdwg_lvl3_locality[species_at_lvl3$tdwg_lvl3_locality == "Kirgizstan"] <- "Kirgizistan"
species_at_lvl3$tdwg_lvl3_locality[species_at_lvl3$tdwg_lvl3_locality == "Kirgizstan"] <- "Kirgizistan"
species_at_lvl3$tdwg_lvl3_locality[species_at_lvl3$tdwg_lvl3_locality == "Cocos (Keeling) Is."] <- "Cocos (Keeling) I."
species_at_lvl3$tdwg_lvl3_locality[species_at_lvl3$tdwg_lvl3_locality == "Gambia"] <- "Gambia, The"


# Genera at level 3
genera_at_lvl3 <- species_at_lvl3
genera_at_lvl3$record_species_name <- gsub(" .*", "", genera_at_lvl3$record_species_name)
colnames(genera_at_lvl3)[1] <- 'record_genus_names'
genera_at_lvl3 <- unique(genera_at_lvl3[, c('record_genus_names','tdwg_lvl3_locality')])
# Some manual edits
genera_at_lvl3 <- rbind(genera_at_lvl3, c("Halfordia", "Queensland"), c("Halfordia", "New South Wales"))



############################################
### Create the genus level 1 area matrix ###
############################################
genera_at_lvl1_matrix <- genera_at_lvl1 %>% 
  group_by(record_genus_names) %>%
  table()

# ... convert to dataframe, remove antarctic region
genera_at_lvl1_df <- as.data.frame.matrix(genera_at_lvl1_matrix) %>% select(!ANTARCTIC)
genera_at_lvl1_df[genera_at_lvl1_df == 0] <- NA
genera_at_lvl1_df$AFRICA[which(genera_at_lvl1_df$AFRICA == 1)] <- "AFRICA"
genera_at_lvl1_df$`ASIA-TEMPERATE`[which(genera_at_lvl1_df$`ASIA-TEMPERATE` == 1)] <- "ASIA-TEMPERATE"
genera_at_lvl1_df$`ASIA-TROPICAL`[which(genera_at_lvl1_df$`ASIA-TROPICAL` == 1)] <- "ASIA-TROPICAL"
genera_at_lvl1_df$AUSTRALASIA[which(genera_at_lvl1_df$AUSTRALASIA == 1)] <- "AUSTRALASIA"
genera_at_lvl1_df$EUROPE[which(genera_at_lvl1_df$EUROPE == 1)] <- "EUROPE"
genera_at_lvl1_df$`NORTHERN AMERICA`[which(genera_at_lvl1_df$`NORTHERN AMERICA` == 1)] <- "NORTHERN AMERICA"
genera_at_lvl1_df$PACIFIC[which(genera_at_lvl1_df$PACIFIC == 1)] <- "PACIFIC"
genera_at_lvl1_df$`SOUTHERN AMERICA`[which(genera_at_lvl1_df$`SOUTHERN AMERICA` == 1)] <- "SOUTHERN AMERICA"

# ... set row names same as tree (so that they include asterisks)
row.names(genera_at_lvl1_df) <- na.omit(rut_tree$tip.label[match(row.names(genera_at_lvl1_df), gsub("[*]" ,"", rut_tree$tip.label))])
# ... check overlap between df and tree tips
name.check(rut_tree, genera_at_lvl1_df)

# ... change to factors
genera_at_lvl1_df[sapply(genera_at_lvl1_df, is.character)] <- lapply(genera_at_lvl1_df[sapply(genera_at_lvl1_df, is.character)], as.factor)
# ... order by rownames
genera_at_lvl1_df <- genera_at_lvl1_df[order(row.names(genera_at_lvl1_df)),]
# ... final check overlap between df and tree tips
name.check(rut_tree,genera_at_lvl1_df)

##############################################
### Create the species level 1 area matrix ###
##############################################

# Get the area matrix
species_at_lvl1_matrix <- species_at_lvl1 %>% 
  group_by(record_species_name) %>%
  table()

# ... convert to dataframe, remove antarctic region
species_at_lvl1_df <- as.data.frame.matrix(species_at_lvl1_matrix) %>% select(!ANTARCTIC)
species_at_lvl1_df[species_at_lvl1_df == 0] <- NA
species_at_lvl1_df$AFRICA[which(species_at_lvl1_df$AFRICA == 1)] <- "AFRICA"
species_at_lvl1_df$`ASIA-TEMPERATE`[which(species_at_lvl1_df$`ASIA-TEMPERATE` == 1)] <- "ASIA-TEMPERATE"
species_at_lvl1_df$`ASIA-TROPICAL`[which(species_at_lvl1_df$`ASIA-TROPICAL` == 1)] <- "ASIA-TROPICAL"
species_at_lvl1_df$AUSTRALASIA[which(species_at_lvl1_df$AUSTRALASIA == 1)] <- "AUSTRALASIA"
species_at_lvl1_df$EUROPE[which(species_at_lvl1_df$EUROPE == 1)] <- "EUROPE"
species_at_lvl1_df$`NORTHERN AMERICA`[which(species_at_lvl1_df$`NORTHERN AMERICA` == 1)] <- "NORTHERN AMERICA"
species_at_lvl1_df$PACIFIC[which(species_at_lvl1_df$PACIFIC == 1)] <- "PACIFIC"
species_at_lvl1_df$`SOUTHERN AMERICA`[which(species_at_lvl1_df$`SOUTHERN AMERICA` == 1)] <- "SOUTHERN AMERICA"
# ... change to factors
species_at_lvl1_df[sapply(species_at_lvl1_df, is.character)] <- lapply(species_at_lvl1_df[sapply(species_at_lvl1_df, is.character)], as.factor)

############################################
### Create the genus level 3 area matrix ###
############################################

# Get the area matrix
genera_at_lvl3_matrix <- genera_at_lvl3 %>% 
  group_by(record_genus_names) %>%
  table()

# ... convert to dataframe
genera_at_lvl3_df <- as.data.frame.matrix(genera_at_lvl3_matrix)

##############################################
### Create the species level 3 area matrix ###
##############################################

# Get the area matrix
species_at_lvl3_matrix <- species_at_lvl3 %>% 
  group_by(record_species_name) %>%
  table()

# ... convert to dataframe
species_at_lvl3_df <- as.data.frame.matrix(species_at_lvl3_matrix)



##################################################################
### Create the species counts for tdwg lvl3 regions (for maps) ###
##################################################################

lvl3_species_counts <- sort(table(species_at_lvl3$tdwg_lvl3_locality), decreasing = TRUE)
all_lvl3_regions <- c(names(lvl3_species_counts), setdiff(tdwg_level3$LEVEL3_NAM, names(lvl3_species_counts)))
all_lvl3_spnums <- rep(NA, length(all_lvl3_regions))
all_lvl3_spnums[1:length(lvl3_species_counts)] <- lvl3_species_counts
tdwg_level3$species_number <- NA
tdwg_level3$species_number[match(all_lvl3_regions, tdwg_level3$LEVEL3_NAM)] <- all_lvl3_spnums

lvl3_endemicspecies_counts <- species_at_lvl3 %>% filter(!record_species_name %in% species_at_lvl3$record_species_name[which(duplicated(species_at_lvl3$record_species_name))])
lvl3_endemicspecies_counts <- sort(table(lvl3_endemicspecies_counts$tdwg_lvl3_locality), decreasing = TRUE)
all_lvl3_regions <- c(names(lvl3_endemicspecies_counts), setdiff(tdwg_level3$LEVEL3_NAM, names(lvl3_endemicspecies_counts)))
all_lvl3_spnums <- rep(NA, length(all_lvl3_regions))
all_lvl3_spnums[1:length(lvl3_endemicspecies_counts)] <- lvl3_endemicspecies_counts
tdwg_level3$endemic_species_number <- NA
tdwg_level3$endemic_species_number[match(all_lvl3_regions, tdwg_level3$LEVEL3_NAM)] <- all_lvl3_spnums

################################################################
### Create the genus counts for tdwg lvl3 regions (for maps) ###
################################################################

lvl3_genus_counts <- sort(table(genera_at_lvl3$tdwg_lvl3_locality), decreasing = TRUE)
all_lvl3_regions <- c(names(lvl3_genus_counts), setdiff(tdwg_level3$LEVEL3_NAM, names(lvl3_genus_counts)))
all_lvl3_gnums <- rep(NA, length(all_lvl3_regions))
all_lvl3_gnums[1:length(lvl3_genus_counts)] <- lvl3_genus_counts
tdwg_level3$genus_number <- NA
tdwg_level3$genus_number[match(all_lvl3_regions, tdwg_level3$LEVEL3_NAM)] <- all_lvl3_gnums

lvl3_endemicgenus_counts <- genera_at_lvl3 %>% filter(!record_genus_names %in% genera_at_lvl3$record_genus_names[which(duplicated(genera_at_lvl3$record_genus_names))])
lvl3_endemicgenus_counts <- sort(table(lvl3_endemicgenus_counts$tdwg_lvl3_locality), decreasing = TRUE)
all_lvl3_regions <- c(names(lvl3_endemicgenus_counts), setdiff(tdwg_level3$LEVEL3_NAM, names(lvl3_endemicgenus_counts)))
all_lvl3_gnums <- rep(NA, length(all_lvl3_regions))
all_lvl3_gnums[1:length(lvl3_endemicgenus_counts)] <- lvl3_endemicgenus_counts
tdwg_level3$endemic_genus_number <- NA
tdwg_level3$endemic_genus_number[match(all_lvl3_regions, tdwg_level3$LEVEL3_NAM)] <- all_lvl3_gnums


## Find the centroids for the tdwg lvl1 multipolygons, write them into the df
x<-terra::centroids(vect(tdwg_level1), inside=FALSE)
tdwg_level1$centroid_long <- crds(x)[,1]
tdwg_level1$centroid_lat <- crds(x)[,2]

## Find the centroids for the tdwg lvl3 multipolygons, write them into the df
x<-terra::centroids(vect(tdwg_level3), inside=FALSE)
tdwg_level3$centroid_long <- crds(x)[,1]
tdwg_level3$centroid_lat <- crds(x)[,2]



################
### Plotting ###
################


#### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
## Plot 1 - Plot of TDWG LVL 1 regions for showing alongside genus level phylogeny. ##
#### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###

lvl1_colours <- c('#77AADD', '#EEDD88', '#FFAABB', '#99DDFF',
                  '#7662E4', '#FF5E29', '#73EF6E', '#359E80') # from https://personal.sron.nl/~pault/data/tol_colors.py -'light'

# map plot
map_lvl1 <- ggplot() +
  geom_sf(data = tdwg_level1, aes(fill = LEVEL1_NAM), colour = "black") +
  scale_fill_manual(values = lvl1_colours) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(limits = c(-60, 90), expand = c(0, 0)) +
  theme(axis.text.y = element_text(size=8),
        axis.text.x = element_text(size=8),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.line.y = element_blank(),
        axis.line.x = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_line(colour = "#E1E1E1", linewidth = 0.5),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5),
        legend.position = "none")

# phylogeny plot
p <- ggtree(rut_tree, branch.length = "none") + 
  scale_y_reverse() + 
  geom_tiplab(size = 2, offset = 1.8, align = T, linesize = 0.25, fontface = 3)
phylo_lvl1 <- gheatmap(p, data = genera_at_lvl1_df, offset=0, width = 0.075, colnames=T) + scale_fill_manual(values = lvl1_colours, na.value = "white")

lvl1_phylo_plot <- grid.arrange(phylo_lvl1, map_lvl1, nrow = 2)


#### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
## Plot 2 - Map, network and upset plots of genus and species richness by TDWG level 1.  ##
#### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###

# Load packages
library(devtools)
library(ComplexHeatmap)
library(UpSetR)

names(lvl1_colours) <- c("AFRICA", "ASIA-TEMPERATE", "ASIA-TROPICAL", "AUSTRALASIA", "EUROPE", "NORTHERN AMERICA", "PACIFIC", "SOUTHERN AMERICA")

## First: genus level
# Reformat the df to suit the upset function
a <- genera_at_lvl1_df
a$taxon = row.names(a)
a[sapply(a, is.factor)] <- lapply(a[sapply(a, is.factor)], as.numeric)
a[is.na(a) == TRUE] <- 0

# First plot to get the set order, as we want to reverse it in the final plot
b <- upset(a, nsets = 8, nintersects = NA)

c <- upset(a, nsets = 8, nintersects = NA, 
      sets.bar.color = lvl1_colours[match(rev(b$Set_names), names(lvl1_colours))],
      keep.order = T,
      set_size.show = T,
      sets = rev(b$Set_names),
      mainbar.y.label = "Number of genera across regions (n = 150)",
      sets.x.label = "Number of genera in each region")

# Need to use cowplot to produce a ggplot-compatible version
upset_genus_lvl1 <- cowplot::plot_grid(NULL, c$Main_bar, c$Sizes, c$Matrix,
                           nrow=2, align='hv', rel_heights = c(3,1),
                           rel_widths = c(2,3))


## Second: species level
# Reformat the df to suit the upset function
a <- species_at_lvl1_df
a$taxon = row.names(a)
a[sapply(a, is.factor)] <- lapply(a[sapply(a, is.factor)], as.numeric)
a[is.na(a) == TRUE] <- 0

# First plot to get the set order, as we want to reverse it in the final plot
b <- upset(a, nsets = 8, nintersects = NA)

c <- upset(a, nsets = 8, nintersects = NA, 
      sets.bar.color = lvl1_colours[match(rev(b$Set_names), names(lvl1_colours))],
      keep.order = T,
      set_size.show = T,
      sets = rev(b$Set_names),
      mainbar.y.label = "Number of species across regions (n = 2126)",
      sets.x.label = "Number of species in each region")

# Need to use cowplot to produce a ggplot-compatible version
upset_species_lvl1 <- cowplot::plot_grid(NULL, c$Main_bar, c$Sizes, c$Matrix,
                                       nrow=2, align='hv', rel_heights = c(3,1),
                                       rel_widths = c(2,3))


lvl1_upset_plots <- grid.arrange(map_lvl1, as.grob(upset_genus_lvl1), upset_species_lvl1, nrow = 3)


## Third: genus level network
df <- !is.na(genera_at_lvl1_df)
co_mat <- t(df) %*% df
# Set diagonal values to 0
diag(co_mat) <- 0

# Assign dim names
dimnames(co_mat) <- list(colnames(df), colnames(df))

# Create graph from adjacency matrix
# ! edge weights are equal to frequency of co-occurrence
g <- graph_from_adjacency_matrix(co_mat, mode = "upper", weighted = TRUE)

# Assign nodes weight equal to species frequency
g <- set.vertex.attribute(g, "v_weight", value = colSums(df))

# Convert the graph to a format usable by ggplot
n <- ggnetwork(g)

# This step changes the x and y coordinates of vertices in the network to 
# those on an octagon (for 8 regions) - this is just a stylistic choice.
starting_angle = 1.9625
radius = 0.5
points = c(1:8)

df <- round(data.frame("x" = c(radius*cos(starting_angle + points*pi/4)),
                       "y" = c(radius*sin(starting_angle + points*pi/4))), digits = 2)

rownames(df) <- c("NORTHERN AMERICA", "SOUTHERN AMERICA", "AFRICA","ASIA-TEMPERATE", "ASIA-TROPICAL","AUSTRALASIA","PACIFIC","EUROPE")

df$region <- rownames(df)

# Swap the octagonal coords into the graph dataframe
for (i in 1:length(df$region)){
  region = df$region[i]
  region_x = unique(n$x[n$name == region])
  region_y = unique(n$y[n$name == region])
  
  n$x[n$x == region_x] <- df$x[i]
  n$y[n$y == region_y] <- df$y[i]
  
  n$xend[n$xend == region_x] <- df$x[i]
  n$yend[n$yend == region_y] <- df$y[i]
}

lvl1_genera_network <- ggplot(n, aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_edges(color = "grey50", curvature = 0.075, aes(linewidth = weight, alpha = weight)) +
  geom_nodes(aes(fill = name, size = v_weight), shape = 21, color = "black") +
  scale_fill_manual(values = lvl1_colours) +
  theme_blank(legend.position = "none")

co_mat[lower.tri(co_mat)] <- NA
genera_network_cooccurrence_matrix <- co_mat

## Fourth: species level network
df <- !is.na(species_at_lvl1_df)
co_mat <- t(df) %*% df
# Set diagonal values to 0
diag(co_mat) <- 0

# Assign dim names
dimnames(co_mat) <- list(colnames(df), colnames(df))

# Create graph from adjacency matrix
# ! edge weights are equal to frequency of co-occurrence
g <- graph_from_adjacency_matrix(co_mat, mode = "upper", weighted = TRUE)

# Assign nodes weight equal to species frequency
g <- set.vertex.attribute(g, "v_weight", value = colSums(df))

# Convert the graph to a format usable by ggplot
n <- ggnetwork(g)

# This step changes the x and y coordinates of vertices in the network to 
# those on an octagon (for 8 regions) - this is just a stylistic choice.
starting_angle = 1.9625
radius = 0.5
points = c(1:8)

df <- round(data.frame("x" = c(radius*cos(starting_angle + points*pi/4)),
                       "y" = c(radius*sin(starting_angle + points*pi/4))), digits = 2)

rownames(df) <- c("NORTHERN AMERICA", "SOUTHERN AMERICA", "AFRICA","ASIA-TEMPERATE", "ASIA-TROPICAL","AUSTRALASIA","PACIFIC","EUROPE")

df$region <- rownames(df)

# Swap the octagonal coords into the graph dataframe
for (i in 1:length(df$region)){
  region = df$region[i]
  region_x = unique(n$x[n$name == region])
  region_y = unique(n$y[n$name == region])
  
  n$x[n$x == region_x] <- df$x[i]
  n$y[n$y == region_y] <- df$y[i]
  
  n$xend[n$xend == region_x] <- df$x[i]
  n$yend[n$yend == region_y] <- df$y[i]
}

lvl1_species_network <- ggplot(n, aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_edges(color = "grey50", curvature = 0.075, aes(linewidth = weight, alpha = weight)) +
  geom_nodes(aes(fill = name, size = v_weight), shape = 21, color = "black") +
  scale_fill_manual(values = lvl1_colours) +
  theme_blank(legend.position = "none")

co_mat[lower.tri(co_mat)] <- NA
species_network_cooccurrence_matrix <- co_mat

lvl1_network_plots <- grid.arrange(lvl1_genera_network, lvl1_species_network, nrow = 2)


#### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
## Plot 3 - Heatmaps of genus & species richness by TDWG level 3. ##
#### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###

lvl3_genus_map <- ggplot() +
  geom_sf(data = tdwg_level3, aes(fill = log(genus_number)), colour = "black") +
  scale_fill_gradient2(low = "#A2E495", mid = "#E0FEFF", high = "#E672A5", midpoint = max(log(na.omit(tdwg_level3$genus_number)))/2) +
  geom_text(data = tdwg_level3, aes(x=centroid_long, y=centroid_lat, label = genus_number)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(limits = c(-60, 90), expand = c(0, 0)) +
  ggtitle(label = paste0("Number of genera by TDWG level 3 regions (n = ", length(rownames(genera_at_lvl1_df)), ")")) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.y = element_text(size=8),
        axis.text.x = element_text(size=8),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.line.y = element_blank(),
        axis.line.x = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_line(colour = "#E1E1E1", linewidth = 0.5),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5)
  )

lvl3_endemic_genus_map <- ggplot() +
  geom_sf(data = tdwg_level3, aes(fill = log(endemic_genus_number)), colour = "black") +
  scale_fill_gradient2(low = "#A2E495", mid = "#E0FEFF", high = "#E672A5", midpoint = max(log(na.omit(tdwg_level3$endemic_genus_number)))/2) +
  geom_text(data = tdwg_level3, aes(x=centroid_long, y=centroid_lat, label = endemic_genus_number)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(limits = c(-60, 90), expand = c(0, 0)) +
  ggtitle(label = paste0("Number of endemic genera by TDWG level 3 regions (n = ", sum(tdwg_level3$endemic_genus_number, na.rm = TRUE), ")")) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.y = element_text(size=8),
        axis.text.x = element_text(size=8),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.line.y = element_blank(),
        axis.line.x = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_line(colour = "#E1E1E1", linewidth = 0.5),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5)
  )

lvl3_species_map <- ggplot() +
  geom_sf(data = tdwg_level3, aes(fill = log(species_number)), colour = "black") +
  scale_fill_gradient2(low = "#A2E495", mid = "#E0FEFF", high = "#E672A5", midpoint = max(log(na.omit(tdwg_level3$species_number)))/2) +
  geom_text(data = tdwg_level3, aes(x=centroid_long, y=centroid_lat, label = species_number)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(limits = c(-60, 90), expand = c(0, 0)) +
  ggtitle(label = paste0("Number of species by TDWG level 3 regions (n = ", length(rownames(species_at_lvl1_df)), ")")) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.y = element_text(size=8),
        axis.text.x = element_text(size=8),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.line.y = element_blank(),
        axis.line.x = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_line(colour = "#E1E1E1", linewidth = 0.5),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5)
  )

lvl3_endemic_species_map <- ggplot() +
  geom_sf(data = tdwg_level3, aes(fill = log(endemic_species_number)), colour = "black") +
  scale_fill_gradient2(low = "#A2E495", mid = "#E0FEFF", high = "#E672A5", midpoint = max(log(na.omit(tdwg_level3$endemic_species_number)))/2) +
  geom_text(data = tdwg_level3, aes(x=centroid_long, y=centroid_lat, label = endemic_species_number)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(limits = c(-60, 90), expand = c(0, 0)) +
  ggtitle(label = paste0("Number of endemic species by TDWG level 3 regions (n = ", sum(tdwg_level3$endemic_species_number, na.rm = TRUE), ")")) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.y = element_text(size=8),
        axis.text.x = element_text(size=8),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.line.y = element_blank(),
        axis.line.x = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_line(colour = "#E1E1E1", linewidth = 0.5),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5)
  )

###############################################################
### Networks and graph clustering with the Leiden Algorithm ###
###############################################################

# First, create some plotting functions
plot_leiden_graph <- function(network, communities){
  pallete <- c("#FF0000", "#0000FF", "#FFFF00", "#008000", "#800080", "#FFA500", "#8B4513", "#FFC0CB", "#00FFFF", "#000000", "#FFFFFF")
  cols <- pallete[1:length(levels(communities))]
  
  ggplot(network, aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_edges(color = "grey50", curvature = 0.075, aes(linewidth = weight), alpha = 0.5) +
    geom_nodes(aes(fill = communities, size = v_weight), shape = 21, color = "black") +
    scale_fill_manual(values = cols) +
    scale_linewidth(range = c(0.02, 4)) +
    theme_blank()
}
plot_leiden_map <- function(map_data, leiden_graph){
  pallete <- c("#FF0000", "#0000FF", "#FFFF00", "#008000", "#800080", "#FFA500", "#8B4513", "#FFC0CB", "#00FFFF", "#000000", "#FFFFFF")
  cols <- pallete[1:leiden_graph$nb_clusters]
  
  map_df <- cbind(map_data, as.factor(leiden_graph$membership[match(map_data$LEVEL3_NAM, leiden_graph$names)]))
  names(map_df)[length(map_df) - 1] <- "membership"
  
  ggplot() +
    geom_sf(data = map_df, aes(fill = membership), colour = "black") +
    scale_fill_manual(values = cols) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(limits = c(-60, 90), expand = c(0, 0)) +
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.y = element_text(size=8),
          axis.text.x = element_text(size=8),
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          axis.line.y = element_blank(),
          axis.line.x = element_blank(),
          panel.background = element_blank(),
          panel.grid.major = element_line(colour = "#E1E1E1", linewidth = 0.5),
          panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5)
    )
}

# RUN FOR GENERA

genera_at_lvl3_df[genera_at_lvl3_df == 0] <- NA

# Remove any cells that are unique/singletons (i.e. they do not interact with other nodes in the network)
df <- genera_at_lvl3_df[which(rowSums(!is.na(genera_at_lvl3_df)) == 1), which(colSums(!is.na(genera_at_lvl3_df) == 1) == 1)]
df_r <- rownames(df)[which(rowSums(!is.na(df)) != 0)]
genera_at_lvl3_df_trimmed <- genera_at_lvl3_df[which(is.na(match(rownames(genera_at_lvl3_df), df_r))),]
df_c <- colnames(df)[which(colSums(!is.na(df)) != 0)]
genera_at_lvl3_df_trimmed <- genera_at_lvl3_df_trimmed[, which(is.na(match(colnames(genera_at_lvl3_df_trimmed), df_c)))]
df <- !is.na(genera_at_lvl3_df_trimmed)
co_mat <- t(df) %*% df
# Set diagonal values to 0
diag(co_mat) <- 0

# Assign dim names
dimnames(co_mat) <- list(colnames(df), colnames(df))

# Create graph from adjacency matrix
# ! edge weights are equal to frequency of co-occurrence
g <- graph_from_adjacency_matrix(co_mat, mode = "upper", weighted = TRUE)
g <- set.vertex.attribute(g, "v_weight", value = colSums(df))


# We'll use seven different resolution parameter values 
g_c_0.0 <- cluster_leiden(g,
                       objective_function = "modularity",
                       resolution_parameter = 0.0,
                       n_iterations = 100)

g_c_0.5 <- cluster_leiden(g,
                          objective_function = "modularity",
                          resolution_parameter = 0.5,
                          n_iterations = 100)

g_c_0.75 <- cluster_leiden(g,
                          objective_function = "modularity",
                          resolution_parameter = 0.75,
                          n_iterations = 100)

g_c_1.0 <- cluster_leiden(g,
                          objective_function = "modularity",
                          resolution_parameter = 1,
                          n_iterations = 100)

g_c_1.25 <- cluster_leiden(g,
                       objective_function = "modularity",
                       resolution_parameter = 1.25,
                       n_iterations = 100)

g_c_1.5 <- cluster_leiden(g,
                           objective_function = "modularity",
                           resolution_parameter = 1.5,
                           n_iterations = 100)
g_c_1.75 <- cluster_leiden(g,
                       objective_function = "modularity",
                       resolution_parameter = 1.75,
                       n_iterations = 100)

# Convert the graph to a format usable by ggplot
g_n <- ggnetwork(g)
g_n$community_g_c_0.0 <- as.factor(g_c_0.0$membership[match(g_n$name, g_c_0.0$names)])
g_n$community_g_c_0.5 <- as.factor(g_c_0.5$membership[match(g_n$name, g_c_0.5$names)])
g_n$community_g_c_0.75 <- as.factor(g_c_0.75$membership[match(g_n$name, g_c_0.75$names)])
g_n$community_g_c_1.0 <- as.factor(g_c_1.0$membership[match(g_n$name, g_c_1.0$names)])
g_n$community_g_c_1.25 <- as.factor(g_c_1.25$membership[match(g_n$name, g_c_1.25$names)])
g_n$community_g_c_1.5 <- as.factor(g_c_1.5$membership[match(g_n$name, g_c_1.5$names)])
g_n$community_g_c_1.75 <- as.factor(g_c_1.75$membership[match(g_n$name, g_c_1.75$names)])


# Plot them...
leiden_graph_genera_0.0 <- plot_leiden_graph(g_n, g_n$community_g_c_0.0)
leiden_map_genera_0.0 <- plot_leiden_map(tdwg_level3, g_c_0.0)

leiden_graph_genera_0.5 <- plot_leiden_graph(g_n, g_n$community_g_c_0.5)
leiden_map_genera_0.5 <- plot_leiden_map(tdwg_level3, g_c_0.5)

leiden_graph_genera_0.75 <- plot_leiden_graph(g_n, g_n$community_g_c_0.75)
leiden_map_genera_0.75 <- plot_leiden_map(tdwg_level3, g_c_0.75)

leiden_graph_genera_1.0 <- plot_leiden_graph(g_n, g_n$community_g_c_1.0)
leiden_map_genera_1.0 <- plot_leiden_map(tdwg_level3, g_c_1.0)

leiden_graph_genera_1.25 <- plot_leiden_graph(g_n, g_n$community_g_c_1.25)
leiden_map_genera_1.25 <- plot_leiden_map(tdwg_level3, g_c_1.25)

leiden_graph_genera_1.5 <- plot_leiden_graph(g_n, g_n$community_g_c_1.5)
leiden_map_genera_1.5 <- plot_leiden_map(tdwg_level3, g_c_1.5)

leiden_graph_genera_1.75 <- plot_leiden_graph(g_n, g_n$community_g_c_1.75)
leiden_map_genera_1.75 <- plot_leiden_map(tdwg_level3, g_c_1.75)


# RUN FOR SPECIES

species_at_lvl3_df[species_at_lvl3_df == 0] <- NA

# Remove any cells that are unique/singletons (i.e. they do not interact with other nodes in the network)
df <- species_at_lvl3_df[which(rowSums(!is.na(species_at_lvl3_df)) == 1), which(colSums(!is.na(species_at_lvl3_df) == 1) == 1)]
df_r <- rownames(df)[which(rowSums(!is.na(df)) != 0)]
species_at_lvl3_df_trimmed <- species_at_lvl3_df[which(is.na(match(rownames(species_at_lvl3_df), df_r))),]
df_c <- colnames(df)[which(colSums(!is.na(df)) != 0)]
species_at_lvl3_df_trimmed <- species_at_lvl3_df_trimmed[, which(is.na(match(colnames(species_at_lvl3_df_trimmed), df_c)))]
df <- !is.na(species_at_lvl3_df_trimmed)
co_mat <- t(df) %*% df
# Set diagonal values to 0
diag(co_mat) <- 0

# Assign dim names
dimnames(co_mat) <- list(colnames(df), colnames(df))

# Create graph from adjacency matrix
# ! edge weights are equal to frequency of co-occurrence
g <- graph_from_adjacency_matrix(co_mat, mode = "upper", weighted = TRUE)
g <- set.vertex.attribute(g, "v_weight", value = colSums(df))

# We'll use seven different resolution parameter values 
g_c_0.0 <- cluster_leiden(g,
                       objective_function = "modularity",
                       resolution_parameter = 0.0,
                       n_iterations = 100)

g_c_0.5 <- cluster_leiden(g,
                          objective_function = "modularity",
                          resolution_parameter = 0.5,
                          n_iterations = 100)

g_c_0.75 <- cluster_leiden(g,
                          objective_function = "modularity",
                          resolution_parameter = 0.75,
                          n_iterations = 100)

g_c_1.0 <- cluster_leiden(g,
                          objective_function = "modularity",
                          resolution_parameter = 1,
                          n_iterations = 100)

g_c_1.25 <- cluster_leiden(g,
                       objective_function = "modularity",
                       resolution_parameter = 1.25,
                       n_iterations = 100)

g_c_1.5 <- cluster_leiden(g,
                           objective_function = "modularity",
                           resolution_parameter = 1.5,
                           n_iterations = 100)
g_c_1.75 <- cluster_leiden(g,
                       objective_function = "modularity",
                       resolution_parameter = 1.75,
                       n_iterations = 100)

# Convert the graph to a format usable by ggplot
g_n <- ggnetwork(g)
g_n$community_g_c_0.0 <- as.factor(g_c_0.0$membership[match(g_n$name, g_c_0.0$names)])
g_n$community_g_c_0.5 <- as.factor(g_c_0.5$membership[match(g_n$name, g_c_0.5$names)])
g_n$community_g_c_0.75 <- as.factor(g_c_0.75$membership[match(g_n$name, g_c_0.75$names)])
g_n$community_g_c_1.0 <- as.factor(g_c_1.0$membership[match(g_n$name, g_c_1.0$names)])
g_n$community_g_c_1.25 <- as.factor(g_c_1.25$membership[match(g_n$name, g_c_1.25$names)])
g_n$community_g_c_1.5 <- as.factor(g_c_1.5$membership[match(g_n$name, g_c_1.5$names)])
g_n$community_g_c_1.75 <- as.factor(g_c_1.75$membership[match(g_n$name, g_c_1.75$names)])

# Plot them...
leiden_graph_species_0.0 <- plot_leiden_graph(g_n, g_n$community_g_c_0.0)
leiden_map_species_0.0 <- plot_leiden_map(tdwg_level3, g_c_0.0)

leiden_graph_species_0.5 <- plot_leiden_graph(g_n, g_n$community_g_c_0.5)
leiden_map_species_0.5 <- plot_leiden_map(tdwg_level3, g_c_0.5)

leiden_graph_species_0.75 <- plot_leiden_graph(g_n, g_n$community_g_c_0.75)
leiden_map_species_0.75 <- plot_leiden_map(tdwg_level3, g_c_0.75)

leiden_graph_species_1.0 <- plot_leiden_graph(g_n, g_n$community_g_c_1.0)
leiden_map_species_1.0 <- plot_leiden_map(tdwg_level3, g_c_1.0)

leiden_graph_species_1.25 <- plot_leiden_graph(g_n, g_n$community_g_c_1.25)
leiden_map_species_1.25 <- plot_leiden_map(tdwg_level3, g_c_1.25)

leiden_graph_species_1.5 <- plot_leiden_graph(g_n, g_n$community_g_c_1.5)
leiden_map_species_1.5 <- plot_leiden_map(tdwg_level3, g_c_1.5)

leiden_graph_species_1.75 <- plot_leiden_graph(g_n, g_n$community_g_c_1.75)
leiden_map_species_1.75 <- plot_leiden_map(tdwg_level3, g_c_1.75)


######################
### Save the plots ###
######################

ggsave(filename = "04_raw_plots/Fig_lvl1_phylo_plot.pdf", plot = lvl1_phylo_plot, device = "pdf", width = 16.4, height = 23.4)
ggsave(filename = "04_raw_plots/Fig_upset_plots.pdf", plot = lvl1_upset_plots, device = "pdf", width = 16.4, height = 23.4)
ggsave(filename = "04_raw_plots/Fig_network_plots.pdf", plot = lvl1_network_plots, device = "pdf", width = 6, height = 12)
ggsave(filename = "04_raw_plots/Fig_lvl3_genera_map.pdf", plot = lvl3_genus_map, device = "pdf", height = 16.4, width = 23.4)
ggsave(filename = "04_raw_plots/Fig_lvl3_endemic_genera_map.pdf", plot = lvl3_endemic_genus_map, device = "pdf", height = 16.4, width = 23.4)
ggsave(filename = "04_raw_plots/Fig_lvl3_species_map.pdf", plot = lvl3_species_map, device = "pdf", height = 16.4, width = 23.4)
ggsave(filename = "04_raw_plots/Fig_lvl3_endemic_species_map.pdf", plot = lvl3_endemic_species_map, device = "pdf", height = 16.4, width = 23.4)
ggsave(filename = "04_raw_plots/Fig_leiden_graph_genera_plot_0_75.pdf", plot = grid.arrange(leiden_graph_genera_0.75, leiden_map_genera_0.75, nrow = 2), device = "pdf", width = 8.2, height = 11.6)
ggsave(filename = "04_raw_plots/Fig_leiden_graph_genera_plot_1_75.pdf", plot = grid.arrange(leiden_graph_genera_1.75, leiden_map_genera_1.75, nrow = 2), device = "pdf", width = 8.2, height = 11.6)
ggsave(filename = "04_raw_plots/Fig_leiden_graph_species_plot_0_75.pdf", plot = grid.arrange(leiden_graph_species_0.75, leiden_map_species_0.75, nrow = 2), device = "pdf", width = 8.2, height = 11.6)
ggsave(filename = "04_raw_plots/Fig_leiden_graph_species_plot_1_75.pdf", plot = grid.arrange(leiden_graph_species_1.75, leiden_map_species_1.75, nrow = 2), device = "pdf", width = 8.2, height = 11.6)

# Save the different leiden threshold to a specific folder to make gifs
ifelse(!dir.exists(file.path("04_raw_plots/_leiden_gif_files/genera")), dir.create(file.path("04_raw_plots/_leiden_gif_files/genera")), FALSE)
ifelse(!dir.exists(file.path("04_raw_plots/_leiden_gif_files/species")), dir.create(file.path("04_raw_plots/_leiden_gif_files/species")), FALSE)

ggsave(filename = "04_raw_plots/_leiden_gif_files/genera/leiden_graph_genera_0.00.pdf", plot = grid.arrange(leiden_graph_genera_0.0, leiden_map_genera_0.0, nrow = 2), width = 8, height = 11)
ggsave(filename = "04_raw_plots/_leiden_gif_files/genera/leiden_graph_genera_0.50.pdf", plot = grid.arrange(leiden_graph_genera_0.5, leiden_map_genera_0.5, nrow = 2), width = 8, height = 11)
ggsave(filename = "04_raw_plots/_leiden_gif_files/genera/leiden_graph_genera_0.75.pdf", plot = grid.arrange(leiden_graph_genera_0.75, leiden_map_genera_0.75, nrow = 2), width = 8, height = 11)
ggsave(filename = "04_raw_plots/_leiden_gif_files/genera/leiden_graph_genera_1.00.pdf", plot = grid.arrange(leiden_graph_genera_1.0, leiden_map_genera_1.0, nrow = 2), width = 8, height = 11)
ggsave(filename = "04_raw_plots/_leiden_gif_files/genera/leiden_graph_genera_1.25.pdf", plot = grid.arrange(leiden_graph_genera_1.25, leiden_map_genera_1.25, nrow = 2), width = 8, height = 11)
ggsave(filename = "04_raw_plots/_leiden_gif_files/genera/leiden_graph_genera_1.50.pdf", plot = grid.arrange(leiden_graph_genera_1.5, leiden_map_genera_1.5, nrow = 2), width = 8, height = 11)
ggsave(filename = "04_raw_plots/_leiden_gif_files/genera/leiden_graph_genera_1.75.pdf", plot = grid.arrange(leiden_graph_genera_1.75, leiden_map_genera_1.75, nrow = 2), width = 8, height = 11)

ggsave(filename = "04_raw_plots/_leiden_gif_files/species/leiden_graph_species_0.00.pdf", plot = grid.arrange(leiden_graph_species_0.0, leiden_map_species_0.0, nrow = 2), width = 8, height = 11)
ggsave(filename = "04_raw_plots/_leiden_gif_files/species/leiden_graph_species_0.50.pdf", plot = grid.arrange(leiden_graph_species_0.5, leiden_map_species_0.5, nrow = 2), width = 8, height = 11)
ggsave(filename = "04_raw_plots/_leiden_gif_files/species/leiden_graph_species_0.75.pdf", plot = grid.arrange(leiden_graph_species_0.75, leiden_map_species_0.75, nrow = 2), width = 8, height = 11)
ggsave(filename = "04_raw_plots/_leiden_gif_files/species/leiden_graph_species_1.00.pdf", plot = grid.arrange(leiden_graph_species_1.0, leiden_map_species_1.0, nrow = 2), width = 8, height = 11)
ggsave(filename = "04_raw_plots/_leiden_gif_files/species/leiden_graph_species_1.25.pdf", plot = grid.arrange(leiden_graph_species_1.25, leiden_map_species_1.25, nrow = 2), width = 8, height = 11)
ggsave(filename = "04_raw_plots/_leiden_gif_files/species/leiden_graph_species_1.50.pdf", plot = grid.arrange(leiden_graph_species_1.5, leiden_map_species_1.5, nrow = 2), width = 8, height = 11)
ggsave(filename = "04_raw_plots/_leiden_gif_files/species/leiden_graph_species_1.75.pdf", plot = grid.arrange(leiden_graph_species_1.75, leiden_map_species_1.75, nrow = 2), width = 8, height = 11)


#########################################
### Save the data frames to csv files ###
#########################################
# For the phylogeny, upset plots + networks
species_at_lvl1_df_csv <- species_at_lvl1_df
species_at_lvl1_df_csv$taxon = row.names(species_at_lvl1_df_csv)
species_at_lvl1_df_csv[sapply(species_at_lvl1_df_csv, is.factor)] <- lapply(species_at_lvl1_df_csv[sapply(species_at_lvl1_df_csv, is.factor)], as.numeric)
species_at_lvl1_df_csv <- species_at_lvl1_df_csv[1:length(colnames(species_at_lvl1_df_csv))-1]
write.csv(species_at_lvl1_df_csv, "05_data_files/Rut_species_dists_TDWG_lvl1.csv")

genera_at_lvl1_df_csv <- genera_at_lvl1_df
genera_at_lvl1_df_csv$taxon = row.names(genera_at_lvl1_df_csv)
genera_at_lvl1_df_csv[sapply(genera_at_lvl1_df_csv, is.factor)] <- lapply(genera_at_lvl1_df_csv[sapply(genera_at_lvl1_df_csv, is.factor)], as.numeric)
genera_at_lvl1_df_csv <- genera_at_lvl1_df_csv[1:length(colnames(genera_at_lvl1_df_csv))-1]
write.csv(genera_at_lvl1_df_csv, "05_data_files/Rut_genera_dists_TDWG_lvl1.csv")

# Cooccurrence matrices for the networks
write.csv(genera_network_cooccurrence_matrix, "05_data_files/Rut_genera_network_cooccurrence_matrix_TDWG_lvl1.csv")
write.csv(species_network_cooccurrence_matrix, "05_data_files/Rut_species_network_cooccurrence_matrix_TDWG_lvl1.csv")

# For the lvl 3 genus and species maps, and Leiden clustering
write.csv(species_at_lvl3_df, "05_data_files/Rut_species_dists_TDWG_lvl3.csv")
write.csv(genera_at_lvl3_df, "05_data_files/Rut_genera_dists_TDWG_lvl3.csv")
write.csv(data.frame(LEVEL3_NAME = tdwg_level3$LEVEL3_NAM,
                     LEVEL3_CODE = tdwg_level3$LEVEL3_COD,
                     LEVEL1_CODE = tdwg_level3$LEVEL1_COD,
                     Genus_number = tdwg_level3$genus_number,
                     Endemic_genus_number = tdwg_level3$endemic_genus_number,
                     Species_number = tdwg_level3$species_number,
                     Endemic_species_number = tdwg_level3$endemic_species_number),
          "05_data_files/Rut_regional_diversity_and_endemism_TDWG_lvl3.csv")

# ... the names of TDWG in each Leiden cluster
gen_c1 <-  g_n %>% group_by(community_g_c1) %>% summarize(TDWG_lvl3_region = paste(sort(unique(name)),collapse=", "))
gen_c3 <-  g_n %>% group_by(community_g_c3) %>% summarize(TDWG_lvl3_region = paste(sort(unique(name)),collapse=", "))
spp_c1 <- n %>% group_by(community_g_c1) %>% summarize(TDWG_lvl3_region = paste(sort(unique(name)),collapse=", "))
spp_c3 <- n %>% group_by(community_g_c3) %>% summarize(TDWG_lvl3_region = paste(sort(unique(name)),collapse=", "))
write.csv(gen_c1, "05_data_files/Leiden_clusters_genera_0.75_TDWG_lvl3.csv", row.names = F)
write.csv(gen_c3, "05_data_files/Leiden_clusters_genera_1.75_TDWG_lvl3.csv", row.names = F)
write.csv(spp_c1, "05_data_files/Leiden_clusters_species_0.75_TDWG_lvl3.csv", row.names = F)
write.csv(spp_c3, "05_data_files/Leiden_clusters_species_1.75_TDWG_lvl3.csv", row.names = F)



