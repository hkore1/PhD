library(rstudioapi)
library(dplyr)
library(ggsankey)
library(ggplot2)
library(xlsx)

# Set workdir to source file location
setwd(dirname(getActiveDocumentContext()$path))

# Load genera from classification spreadsheet
classification_df <- read.xlsx(file = "Table1.1_Classification_of_Rutaceae_genera.xlsx", sheetIndex = 1)

classification_df_subset <- classification_df %>% select(Genus, Engler_Subfamily, Kubitzki_Subfamily, Groppo_Subfamily, Morton_Subfamily, Appelhans_Subfamily)

# Create plot dataframe (https://stackoverflow.com/questions/74141362/how-to-skip-nodes-with-na-value-in-ggsankey)
plot_df <- do.call(rbind, apply(classification_df_subset, 1, function(x) {
  x <- na.omit(x[-1])
  data.frame(x = names(x), node = x, 
             next_x = dplyr::lead(names(x)), 
             next_node = dplyr::lead(x), row.names = NULL)
})) %>%
  mutate(x = factor(x, names(classification_df_subset)[-1]),
         next_x = factor(next_x, names(classification_df_subset)[-1]))

# Change the order of factor levels to plot nicely.
plot_df$node <- factor(plot_df$node, levels = c("Amyridoideae", "Haplophylloideae", "Zanthoxyloideae", "Rutoideae", "Flindersioideae", "Toddalioideae", "Aurantioideae", "Cneoroideae", "Dictyolomatoideae", "Spathelioideae"))


col_pallete <- c("#117733", "#372884", "#44AA99", "#88CCEE", "#DDCC77", "#CC6677", "#AA4499", "#882255", "#F5F1AE", "#8DE84C")

sankey_plot <- ggplot(plot_df, aes(x = x,
           next_x = next_x,
           node = node,
           next_node = next_node,
           fill = node,,
           label = node)) +
  geom_sankey(flow.alpha = 0.5,
              node.color = NA,
              show.legend = TRUE) +
  geom_sankey_text(size = 3, color = "black", fill = NA, hjust = 0, 
                   position = position_nudge(x = 0.1), aes(
                     x = as.numeric(x) + .05,
                     label = after_stat(paste0(node, "\nn = ", freq))
                   )) +
  scale_fill_manual(values = col_pallete) +
  theme_void()
sankey_plot

ggsave(filename = "Fig1.2_Sankey_diagram_unedited.pdf", plot = sankey_plot, device = "pdf")

