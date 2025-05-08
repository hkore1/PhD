library(rnaturalearth)
library(raster)
library(terra)
library(crsuggest)
library(sf)
library(ggplot2)
library(dplyr)
library(edgebundle)
library(gridExtra)
library(ragg)
library(rstudioapi)

# Set workdir to source file location
setwd(dirname(getActiveDocumentContext()$path))

########################################################################################################
### Set up map data + hex grid and the dataframe that contains info on what hex each record falls in ###
########################################################################################################

optimal_CRS = 3112

# Set up shapefiles for mapping
map_countries <- ne_countries(type = "countries", country = c("Australia", "Indonesia", "Papua New Guinea",
                                                        "New Caledonia", "New Zealand", 
                                                        "Vanuatu", "Solomon Islands", "Fiji"),
                        scale = "large", returnclass = "sf")

map_dependencies <- ne_countries(type = "map_units", geounit = c("Norfolk Island"),
                              scale = "large", returnclass = "sf")

map_all <- rbind(map_countries, map_dependencies)

map_all <- map_all %>% 
  st_transform(optimal_CRS) %>%
  summarise()

# Make region valid
region <- st_make_valid(map_all)

# Create a grid for the region
grid <- region %>%
  st_make_grid(cellsize = units::as_units(10000, "km^2"), square = FALSE, flat_topped = TRUE) %>% # Specify cells with 10,000 km2 area (https://github.com/r-spatial/sf/issues/1505)
  st_intersection(region)

# Read in point data from GBIF
point_data_GBIF <- readRDS("slim_gbif_data_regions_filtered.rds") %>%
  filter(tdwg_lvl1_locality %in% c("ASIA-TROPICAL", "AUSTRALASIA", "PACIFIC"))

# Read in point data from AVH
point_data_AVH <- read.csv("records-2025-04-03/records-2025-04-03.csv")

# Read in the distribution dictionary and filter out non-native species from the AVH data
dist_dict <- readRDS("distribution_dictionary.rds") %>% 
  filter(tdwg_lvl1_location %in% c("ASIA-TROPICAL", "AUSTRALASIA", "PACIFIC")) %>%
  distinct(species, tdwg_lvl3_location, .keep_all = TRUE)

point_data_AVH <-  point_data_AVH[which(point_data_AVH$species %in% dist_dict$species) ,] # Remove any species not in the dictionary

# Make variables the same in both datasets
x <- point_data_AVH[intersect(colnames(point_data_GBIF), colnames(point_data_AVH))]
y <- point_data_GBIF[intersect(colnames(point_data_GBIF), colnames(point_data_AVH))]

# Combine 'em
point_data <- rbind(x, y)

# Remove any records where stateProvince isn't in the dist_dict for AU collections
z <- point_data %>%
  filter(countryCode == "AU")

for (i in 1:length(rownames(z))){
  entry = z[i,]
  
  dist_locs <- subset(dist_dict, species == entry$species)$tdwg_lvl3_location
  
  if (entry$stateProvince %in% dist_locs == FALSE){
    z$stateProvince[i] <- "REMOVE"
  }
}

z <- z %>%
  filter(stateProvince != "REMOVE")

point_data <- rbind(z, filter(point_data, countryCode != "AU"))

# Attempt to filter out duplicates - this is not super robust but we'll only count the presence/absence of a taxon from a grid cell so the fact there could be duplicates doesn't really matter
point_data_no_dups <- point_data %>%
  distinct(catalogNumber, decimalLatitude, decimalLongitude, .keep_all = TRUE) %>% # Dups 
  filter(!decimalLatitude == "" | !decimalLongitude == "") # Records need coordinates

# Further cleaning
point_data_no_dups <- point_data_no_dups[grep("garden", point_data_no_dups$locality, ignore.case = TRUE, invert = TRUE) ,] # Remove any record with "garden" in the locality
point_data_no_dups <- point_data_no_dups[grep("cultivated", point_data_no_dups$locality, ignore.case = TRUE, invert = TRUE) ,] # Remove any record with "cultivated" in the locality
point_data_no_dups <- point_data_no_dups %>%
  filter(!(countryCode == "NZ" & !species %in% c("Leionema nudum","Melicope ternata","Melicope simplex"))) # Keep only NZs endemics, get rid of wrong things for NZ

# Specific wrong points manually detected in first pass
point_data_no_dups <- point_data_no_dups %>%
  filter(!(catalogNumber %in% c("AD 98585667", "CANB 252964.1", "PERTH 979791", "AD 97130016", "AK50234", "AD 97613138C", "AK50171", "MEL 0067588A", "SP102531", "MEL 0062097A", "NE 50492", "NSW393250", "MELUD112771a", "NSW10359", "NSW777654", "NSW779024", "CHR 688773", "QRS 6956.1", "QRS 120816.1", "QRS 120816.2")))

# Write data collection points to new variable
pnts <- data.frame(
  "x" = point_data_no_dups$decimalLongitude,
  "y" = point_data_no_dups$decimalLatitude)

# create a points collection
pnts_sf <- do.call("st_sfc",c(lapply(1:nrow(pnts), 
                                     function(i) {st_point(as.numeric(pnts[i, ]))}), list("crs" = 4326))) 

pnts_trans <- st_transform(pnts_sf, crs = optimal_CRS) # apply transformation to pnts sf
map_trans <- grid # apply transformation to polygons sf

# intersect and extract state name
pnts$name <- (st_intersects(pnts_trans, map_trans, sparse = T))
pnts$name <- unlist(lapply(pnts$name, function(x) if(identical(x, integer(0))) NA_character_ else x))

# append original data to new dataframe
pnts <- cbind(pnts, point_data_no_dups)

# remove points that do not intersect with the grid
pnts <- pnts[!is.na(pnts$name), ]





###############################################################
### Now get in to plotting different variations of the data ###
###############################################################

##----------------------------##
## First: plot all species... ##
##----------------------------##

# View number of collection records for each species
RecordCount <- pnts %>%
  group_by(genus) %>%
  count(species)

# Get hexes where each species occur
RegionOccurrences <- pnts %>%
  group_by(species) %>%
  distinct(name) 

# Get number of species in each hex
RegionFrequency <- RegionOccurrences %>%
  group_by(name) %>%
  count()

# Create a new variable for frequency for each hex
grid <- as.data.frame(grid)
grid$hexNumber <- as.integer(row.names(grid))
grid$freq <- NA

# Create list 'x', which matches hexNumber in 'grid' to their index in 'SubregionFrequency'
x <- lapply(grid$hexNumber, function(x) which(RegionFrequency$name %in% x))

## Converts int(0) to NA for values in 'x'
#... find zero-length values
idx <- !(sapply(x, length))
#... replace these values with NA
x[idx] <- NA

# Unlist 'x' to convert to vector
x <- unlist(x)

# Writes the frequency values from 'SubregionFrequency' to 'grid$freq' based on the indexes in 'x'
grid$freq <- RegionFrequency$n[x]

# Remove NA hexes, transform the data freq column a bit but I don't end up using these.
grid <- grid[which(!is.na(grid$freq)) ,]
grid$ln_freq <- log10(grid$freq)
grid$sqrt_freq <- sqrt(grid$freq)

# Colour palette for map plot
colScalefunc <- colorRampPalette(c("#E6FD06", "#E84E4E", "#572479"))
ColourPalette <- colScalefunc(max(grid$freq))

hexMap <- ggplot() +
  geom_sf(data = region, fill = "grey20") +
  geom_sf(data = grid, aes(geometry = geometry, fill = freq), colour = "grey20", linewidth = 0.05) +
  scale_fill_gradientn(colours = c(ColourPalette[1:max(RegionFrequency$n)]), name = "Number of \nspecies") +
  labs(x = "Longitude", y = "Latitude") +
  xlim(st_bbox(grid$geometry)[1], st_bbox(grid$geometry)[3]) +
  ylim(st_bbox(grid$geometry)[2], st_bbox(grid$geometry)[4]) +
  theme(axis.title = element_text(hjust = 0.5, face = "bold", family = "XCharter"),
        axis.text = element_text(family = "XCharter"),
        legend.title = element_text(face = "bold", family = "XCharter"),
        legend.text = element_text(family = "XCharter"),
        panel.background = element_rect(fill = "white"),
        panel.grid = element_line(colour = "#E6E6E6"),
        panel.border = element_rect(fill = "NA", colour = "black"))


##--------------------------------------##
## Second: plot network and communities ##
##--------------------------------------##

## Find the centroids for the grid hexes multipolygons, write them into the grid df
x<-st_coordinates(st_transform(st_centroid(grid$geometry), optimal_CRS))
grid$centroid_long <- x[,1]
grid$centroid_lat <- x[,2]

a <- RegionOccurrences %>% 
  group_by(species) %>%
  table()

# ... convert to dataframe, remove antarctic region
b <- as.data.frame.matrix(a)

b[b == 0] <- NA

# Remove any cells that are unique/singletons (i.e. they do not interact with other nodes in the network)
df <- b[which(rowSums(!is.na(b)) == 1), which(colSums(!is.na(b) == 1) == 1)]
df_r <- rownames(df)[which(rowSums(!is.na(df)) != 0)]
genera_at_lvl3_df_trimmed <- b[which(is.na(match(rownames(b), df_r))),]
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
library(igraph)
library(ggnetwork)
plot_leiden_graph <- function(network, communities){
  pallete <- c("#FFA500", "#0000FF", "#FF0000", "#008000", "#800080", "#FFA500", "#8B4513", "#FFC0CB", "#00FFFF", "#000000", "#FFFFFF")
  cols <- pallete[1:length(levels(communities))]
  
  ggplot(network, aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_edges(color = "grey50", curvature = 0.075, aes(linewidth = weight), alpha = 0.5, show.legend = FALSE) +
    geom_nodes(aes(fill = communities, size = v_weight), shape = 21, color = "black", stroke = 0.1, show.legend = FALSE) +
    scale_fill_manual(values = cols) +
    scale_linewidth(range = c(0.02, 4)) +
    theme_blank() +
    theme(panel.border = element_rect(fill = "NA", colour = "black"))
}

g <- igraph::graph_from_adjacency_matrix(co_mat, mode = "upper", weighted = TRUE)
g <- igraph::set.vertex.attribute(g, "v_weight", value = colSums(df))
g_n <- ggnetwork(g)

# We'll use X different resolution parameter values (0.15, 0.5, 1, 1.5)
g_c <- cluster_leiden(g,
                          objective_function = "modularity",
                          resolution_parameter = 0.5,
                          n_iterations = 100)
g_n$community_g_c<- as.factor(g_c$membership[match(g_n$name, g_c$names)])

communities <- g_n$community_g_c
names(communities) <- g_n$name
communities_dict <- communities[duplicated(names(communities)) == FALSE]

leiden_graph <- plot_leiden_graph(g_n, g_n$community_g_c)

# named geo coordinate vectors...
lat <- grid$centroid_lat
names(lat) <- grid$hexNumber
lon <- grid$centroid_long
names(lon) <- grid$hexNumber
# ... as a (x, y) coordinate matrix
geo <- cbind(lon[ V(g)$name ], lat[ V(g)$name ])

V(g)$community <- communities_dict[match(V(g)$name, names(communities_dict))]

V(g)$longitude <- geo[,1]
V(g)$latitude <- geo[,2]

xy <- cbind(V(g)$longitude, V(g)$latitude)
verts <- data.frame(x = V(g)$longitude, y = V(g)$latitude, community = V(g)$community, weight = V(g)$v_weight)

# Perform the edge-bundling. This bundles edges so the network is less hairball-y.
pbundle <- edge_bundle_path(g, xy, max_distortion = 20, weight_fac = 30, segments = 30)

# Plot the network
netMap <- ggplot() +
  geom_sf(data = region, fill = "grey20") +
  geom_path(data = pbundle, aes(x, y, group = group), col = alpha("#1E88E5", 0.4), linewidth = 0.05) +
  geom_path(data = pbundle, aes(x, y, group = group), col = alpha("red", 0.8), linewidth = 0.005 ) +
  geom_sf(data = region, fill = NA, colour = alpha("black", 0.5)) +
  geom_point(data = verts, aes(x, y), col = "#1E88E5", size = 0.25) +
  geom_point(data = verts, aes(x, y), col = "red", size = 0.25, alpha = 0.75) +
  geom_point(data = verts, aes(x, y), col = "white", size = 0.2, alpha = 0.5) +
  theme(axis.title = element_text(hjust = 0.5, face = "bold", family = "XCharter"),
        axis.text = element_text(family = "XCharter"),
        legend.title = element_text(face = "bold", family = "XCharter"),
        legend.text = element_text(family = "XCharter"),
        panel.background = element_rect(fill = "white"),
        panel.grid = element_line(colour = "#E6E6E6"),
        panel.border = element_rect(fill = "NA", colour = "black"))

# Plot the Leiden cluster communities
LeidenPlot <- ggplot() +
  geom_sf(data = region, colour = "grey20") +
  geom_point(
    data = verts, aes(x, y, fill = community, size = weight), shape = 21, colour = "grey20", stroke = 0.1
  ) +
  scale_fill_manual(values = c("#FFA500", "#0000FF", "#FF0000", "#008000")) +
  annotation_custom(ggplotGrob(leiden_graph), xmin = min(verts$x)*1.1, xmax = min(verts$x)*0.5, ymin = min(verts$y), ymax = min(verts$y)*0.5) +
  labs(x = "Longitude", y = "Latitude") +
  xlim(st_bbox(grid$geometry)[1], st_bbox(grid$geometry)[3]) +
  ylim(st_bbox(grid$geometry)[2], st_bbox(grid$geometry)[4]) +
  theme(axis.title = element_text(hjust = 0.5, face = "bold", family = "XCharter"),
        axis.text = element_text(family = "XCharter"),
        legend.title = element_text(face = "bold", family = "XCharter"),
        legend.text = element_text(family = "XCharter"),
        panel.background = element_rect(fill = "white"),
        panel.grid = element_line(colour = "#E6E6E6"),
        panel.border = element_rect(fill = "NA", colour = "black"))


ggsave(plot = hexMap, filename = "Fig_hexMap_Australasia_all_species.jpg", device = "jpg", width = 12, height = 10)
ggsave(plot = netMap, filename = "Fig_netMap_Australasia.jpg", device = "jpg", width = 12, height = 10)
ggsave(plot = LeidenPlot, filename = "Fig_LeidenHex_Australasia_all_species.jpg", device = "jpg", width = 14, height = 11.5)


##----------------------------------------------##
## Third: plot maps for each Australasian group ##
##----------------------------------------------##

# Specify the groups
Boronia_group <- c("Acronychia", "Boronia", "Brombya", "Comptonella", "Cyanothamnus", "Dutailliopsis", "Dutaillyea", "Euodia", "Maclurodendron", "Medicosma", "Melicope", "Neobyrnesia", "Perryodendron", "Picrella", "Pitaviaster", "Sarcomelicope", "Tetractomia", "Zieria")
Eriostemon_group <- c("Asterolasia", "Chorilaena", "Correa", "Crowea", "Diplolaena", "Drummondita", "Eriostemon", "Geleznowia", "Halfordia", "Leionema", "Muiriantha", "Myrtopsis", "Nematolepis", "Neoschmidia", "Phebalium", "Philotheca")
Flindersia_group <- c("Acradenia", "Bosistoa", "Bouchardatia", "Coatesia", "Crossosperma", "Dinosperma", "Flindersia", "Geijera", "Lunasia", "Pentaceras", "Pitavia")


# Create a grid for the region for each group
Boronia_grid <- region %>%
  st_make_grid(cellsize = units::as_units(10000, "km^2"), square = FALSE, flat_topped = TRUE) %>% # Specify cells with 10,000 km2 area (https://github.com/r-spatial/sf/issues/1505)
  st_intersection(region)
Boronia_grid <- as.data.frame(Boronia_grid)
Boronia_grid$hexNumber <- as.integer(row.names(Boronia_grid))
Boronia_grid$freq <- NA

Eriostemon_grid <- Boronia_grid
Flindersia_grid <- Boronia_grid

# Get the respective counts for each hex for each group
# Get hexes where each species occur
Boronia_RegionFrequency <- pnts %>%
  filter(grepl(paste(Boronia_group, collapse = "|"), species)) %>%
  group_by(species) %>%
  distinct(name) %>%
  group_by(name) %>%
  count()
  
Eriostemon_RegionFrequency <- pnts %>%
  filter(grepl(paste(Eriostemon_group, collapse = "|"), species)) %>%
  group_by(species) %>%
  distinct(name) %>%
  group_by(name) %>%
  count()

Flindersia_RegionFrequency <- pnts %>%
  filter(grepl(paste(Flindersia_group, collapse = "|"), species)) %>%
  group_by(species) %>%
  distinct(name) %>%
  group_by(name) %>%
  count()

# Define function to populate the grids and run it
calc_grid <- function(grid, freq_table){
  # Create list 'x', which matches hexNumber in 'grid' to their index in 'SubregionFrequency'
  x <- lapply(grid$hexNumber, function(x) which(freq_table$name %in% x))
  
  ## Converts int(0) to NA for values in 'x'
  #... find zero-length values
  idx <- !(sapply(x, length))
  #... replace these values with NA
  x[idx] <- NA
  
  # Unlist 'x' to convert to vector
  x <- unlist(x)
  
  # Writes the frequency values from 'SubregionFrequency' to 'grid$freq' based on the indexes in 'x'
  grid$freq <- freq_table$n[x]
  
  # Remove NA hexes
  grid <- grid[which(!is.na(grid$freq)) ,]
  
  return(grid)
}

Boronia_grid <- calc_grid(Boronia_grid, Boronia_RegionFrequency)
Eriostemon_grid <- calc_grid(Eriostemon_grid, Eriostemon_RegionFrequency)
Flindersia_grid <- calc_grid(Flindersia_grid, Flindersia_RegionFrequency)

# Colour palette for map plot
colScalefunc <- colorRampPalette(c("#E6FD06", "#E84E4E", "#572479"))
ColourPalette <- colScalefunc(max(grid$freq))

Boronia_plot <- ggplot() +
  geom_sf(data = region, fill = "grey20") +
  geom_sf(data = Boronia_grid, aes(geometry = geometry, fill = freq), colour = "grey20", linewidth = 0.05) +
  scale_fill_gradientn(colours = c(ColourPalette[1:max(RegionFrequency$n)]), name = "Number of \nspecies") +
  labs(x = "Longitude", y = "Latitude") +
  xlim(st_bbox(Boronia_grid$geometry)[1], st_bbox(Boronia_grid$geometry)[3]) +
  ylim(st_bbox(Boronia_grid$geometry)[2], st_bbox(Boronia_grid$geometry)[4]) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 6),
        panel.background = element_rect(fill = "white"),
        panel.grid = element_line(colour = "#E6E6E6"),
        panel.border = element_rect(fill = "NA", colour = "black"))

Eriostemon_plot <- ggplot() +
  geom_sf(data = region, fill = "grey20") +
  geom_sf(data = Eriostemon_grid, aes(geometry = geometry, fill = freq), colour = "grey20", linewidth = 0.05) +
  scale_fill_gradientn(colours = c(ColourPalette[1:max(RegionFrequency$n)]), name = "Number of \nspecies") +
  labs(x = "Longitude", y = "Latitude") +
  xlim(st_bbox(Boronia_grid$geometry)[1], st_bbox(Boronia_grid$geometry)[3]) +
  ylim(st_bbox(Boronia_grid$geometry)[2], st_bbox(Boronia_grid$geometry)[4]) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 6),
        panel.background = element_rect(fill = "white"),
        panel.grid = element_line(colour = "#E6E6E6"),
        panel.border = element_rect(fill = "NA", colour = "black"))

Flindersia_plot <- ggplot() +
  geom_sf(data = region, fill = "grey20") +
  geom_sf(data = Flindersia_grid, aes(geometry = geometry, fill = freq), colour = "grey20", linewidth = 0.05) +
  scale_fill_gradientn(colours = c(ColourPalette[1:max(RegionFrequency$n)]), name = "Number of \nspecies") +
  labs(x = "Longitude", y = "Latitude") +
  xlim(st_bbox(Boronia_grid$geometry)[1], st_bbox(Boronia_grid$geometry)[3]) +
  ylim(st_bbox(Boronia_grid$geometry)[2], st_bbox(Boronia_grid$geometry)[4]) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 6),
        panel.background = element_rect(fill = "white"),
        panel.grid = element_line(colour = "#E6E6E6"),
        panel.border = element_rect(fill = "NA", colour = "black"))

ggsave(plot = Boronia_plot, filename = "Fig_hexMap_Australasia_Boronia_group.pdf", device = "pdf", width = 12, height = 10)
ggsave(plot = Eriostemon_plot, filename = "Fig_hexMap_Australasia_Eriostemon_group.pdf", device = "pdf", width = 12, height = 10)
ggsave(plot = Flindersia_plot, filename = "Fig_hexMap_Australasia_Flindersia_group.pdf", device = "pdf", width = 12, height = 10)




######### '-'-'-'-'-' EXTRA BITS '-'-'-'-'-' ##########
## Plot the hexNumbers for reference
hexNumMap <- ggplot() +
  geom_sf(data = region, fill = "grey20") +
  geom_sf(data = grid, aes(geometry = geometry), fill = "white", colour = "black", linewidth = 0.1) +
  geom_text(data = grid, aes(label = hexNumber, x = centroid_long, y = centroid_lat), size = 1.2, fontface = 2)
ggsave(plot = hexNumMap, filename = "_useful_complementary_files/hexcodes_mapped.pdf", device = "pdf", width = 12, height = 10)

## save the underlying data to csv
write.csv(df, "_useful_complementary_files/RegionOccurrences_matrix.csv")
write.csv(RegionOccurrences[order(RegionOccurrences$species),], "_useful_complementary_files/RegionOccurrences.csv")


## I just used this bit to calculate numbers for section 1.4.1
mm<-dist_dict %>% filter(tdwg_lvl1_location %in% c("AUSTRALASIA", "PACIFIC") | tdwg_lvl3_location %in% c("New Guinea","Bismarck Archipelago","Solomon Is.","Vanuatu","Maluku"))
mm$genus <- gsub(" .*","", mm$species)

Aurantioideae_Australasian_genera <- c("Bergera", "Citrus", "Luvunga", "Micromelum", "Glycosmis", "Clausena", "Murraya")
Zanthoxyloideae_Australasian_genera <- c("Acradenia", "Acronychia", "Asterolasia", "Boronia", "Bosistoa", "Bouchardatia", "Brombya", "Chorilaena", "Coatesia", "Comptonella", "Correa", "Crossosperma", "Crowea", "Cyanothamnus", "Dinosperma", "Diplolaena", "Drummondita", "Dutailliopsis", "Dutaillyea", "Eriostemon", "Euodia", "Flindersia", "Geijera", "Geleznowia", "Halfordia", "Leionema", "Lunasia", "Medicosma", "Melicope", "Muiriantha", "Myrtopsis", "Nematolepis", "Neobyrnesia", "Neoschmidia", "Pentaceras", "Perryodendron", "Phebalium", "Philotheca", "Picrella", "Pitaviaster", "Sarcomelicope", "Tetractomia", "Zanthoxylum", "Zieria")

Aur_spp <- mm[grep(paste(Aurantioideae_Australasian_genera,collapse="|"), mm$species),] %>% distinct(species, .keep_all = T)
Zan_spp <- mm[grep(paste(Zanthoxyloideae_Australasian_genera,collapse="|"), mm$species),] %>% distinct(species, .keep_all = T)
