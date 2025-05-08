library(rstudioapi)
library(ape)
library(ggtree)
library(xlsx)
library(geiger)

# Set workdir to source file location
setwd(dirname(getActiveDocumentContext()$path))

# Load tree
rut_tree <- read.tree("../1.3_Phytogeography_of_Rutaceae/02_data_preparation/Rutaceae_summary_tree_final.tre")

# Load classification spreadsheet
classification_df <- read.xlsx(file = "Table1.1_Classification_of_Rutaceae_genera.xlsx", sheetIndex = 1)

# Set names in the spreadsheet to be the same as in the tree (i.e. including asterisks)
classification_df$Genus <- rut_tree$tip.label[match(classification_df$Genus, gsub("[*]" ,"", rut_tree$tip.label))]
classification_df <- classification_df[2:10]

name.check(rut_tree, classification_df, data.names=classification_df$Genus)

# Base tree plot
p1 <- ggtree(rut_tree, branch.length = "none")  %<+% classification_df +
  scale_y_reverse() +
  geom_tiplab(size = 2, linesize = 0.25, fontface = 3) 

# Add vertical strips for the Appelhans et al. 2021 subfamily classification
p2 <- p1 + geom_strip("Harrisonia", "Spathelia", label = "Cneoroideae", offset = 4.5, extend = 0.25) +
  geom_strip("Chloroxylon", "Thamnosma", label = "Rutoideae", offset = 4.5, extend = 0.25) +
  geom_strip("Clausena", "Burkillanthus", label = "Aurantioideae", offset = 4.5, extend = 0.25) +
  geom_strip("Casimiroa", "Neoraputia", label = "Zanthoxyloideae", offset = 4.5, extend = 0.25) +
  geom_strip("Decatropis", "Decatropis", label = "Zanthoxyloideae", offset = 4.5, extend = 0.25) +
  geom_strip("Amyris", "Cneoridium", label = "Amyridoideae", offset = 4.5, extend = 0.25) +
  geom_strip("Stauranthus", "Stauranthus", label = "Amyridoideae", offset = 4.5, extend = 0.25) +
  geom_strip("Haplophyllum", "Haplophyllum", label = "Haplophylloideae", offset = 4.5, extend = 0.25)

# Create a color palettes for the Kubitzksi et al. 2011 subfamilies
Kub_subfamily_colors <- c("#FF005D", "#4BABFF", "#FDC418")

# Plot the Kubitzki subfamilies
p3 <- p2 + geom_tippoint(aes(color = Kubitzki_Subfamily), position = position_nudge(x = 4.5, y = 0)) +
  scale_color_manual(values = Kub_subfamily_colors, na.translate = FALSE)

# Plot the Kubitzki alliances
p4 <- p3 + geom_tippoint(aes(shape = Kubitzki_Alliance), position = position_nudge(x = 4.25, y = 0)) +
  scale_shape_manual(values = c(1:12), na.translate = FALSE)

# Save the plot
ggsave(filename = "plots_and_figs/Fig1.1_phylo_plot.pdf", plot = p4, device = "pdf", width = 8.3, height = 11.7)



