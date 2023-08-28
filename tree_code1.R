library(ggtree)
library(readr)
library(ape)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpp)
library(deeptime)
library(ggnewscale)

t <- read.tree('tree_sp.newick')
data <- read_csv('gram_stain_prediction_sp.csv', show_col_types = FALSE)
p <- ggtree(t, layout = 'circular', size = 0.02)

# add tip data for holdoutss ----------------------------------------------
data_holdout <- read_csv('gram_stain_holdout_sp.csv', show_col_types = FALSE)
data_holdout <- data_holdout |> 
    select(NCBI_ID, Attribute, Score) |> 
    distinct() |> 
    pivot_wider(
        names_from = 'Attribute', values_from = 'Score', values_fill = 0
    ) |> 
    as.data.frame() |> 
    rename(node = NCBI_ID) |> 
    relocate(node, .after = last_col()) |> 
    filter(node %in% t$tip.label) |> 
    tibble::column_to_rownames(var = 'node')
data_holdout <- data_holdout[t$tip.label,]
rownames(data_holdout) <- t$tip.label
# data_holdout <- data_holdout[!is.na(rowSums(data_holdout)),]
data_holdout[is.na(data_holdout)] <- 0
colnames(data_holdout) <- gsub(' ', '_', colnames(data_holdout))

# Add data for output of propagation --------------------------------------
tip_data <- data |> 
    filter(Rank == 'species') |> 
    mutate(
        NCBI_ID = paste0(sub('^(\\w).*$', '\\1', Rank), '__', NCBI_ID)
    ) |> 
    select(NCBI_ID, Attribute, Score) |> 
    distinct() |> 
    pivot_wider(
        names_from = 'Attribute', values_from = 'Score'
    ) |> 
    as.data.frame() |> 
    rename(node = NCBI_ID) |> 
    relocate(node, .after = last_col()) |> 
    filter(node %in% t$tip.label) |> 
    tibble::column_to_rownames(var = 'node')
tip_data <- tip_data[t$tip.label,]
rownames(tip_data) <- t$tip.label
tip_data[is.na(tip_data)] <- 0
# tip_data <- tip_data[!is.na(rowSums(tip_data)),]
colnames(tip_data) <- gsub(' ', '_', colnames(tip_data))
tip_data <- tip_data[, colnames(data_holdout)]

p_heatmap_ <- p |> 
    gheatmap(
        data = tip_data, offset = 0.1, width = 0.1, 
        # colnames_angle = 90, colnames_offset_y = .25,,
        colnames = FALSE,
        # legend_title = 'Predicted score',
        color = NA
        # low = 'white', high = 'firebrick'
    ) +
    scale_fill_viridis_c(option = 'D', name = 'Predicted score', na.value = 'white')

p_heatmap <- p_heatmap_ + new_scale_fill()
p_heatmap <- p_heatmap |> 
    gheatmap(
        data = data_holdout, offset = 1, width = 0.1, 
        # colnames_angle = 90, colnames_offset_y = .25,,
        colnames = FALSE,
        color = NA
    ) +
    scale_fill_gradient(
        name = 'Holdout score',
        low = 'white', high = 'dodgerblue4', na.value = 'white', 
    )

# ggsave(
#     filename = 'test_tree_plot.png', plot = p_heatmap, width = 10, height = 10,
#     units = 'in', dpi = 300
# )

# Add inset pies ----------------------------------------------------------
node_data <- data |> 
    filter(! Rank %in% c('species', 'strain')) |> 
    mutate(
        NCBI_ID = paste0(sub('^(\\w).*$', '\\1', Rank), '__', NCBI_ID)
    ) |> 
    select(NCBI_ID, Attribute, Score) |> 
    distinct() |> 
    pivot_wider(
        names_from = 'Attribute', values_from = 'Score'
    ) |> 
    as.data.frame() |> 
    rename(node = NCBI_ID) |> 
    relocate(node, .after = last_col()) |> 
    filter(node %in% t$node.label) |> 
    tibble::column_to_rownames(var = 'node') 
node_data <- node_data[t$node.label,] # I do this to make sure that I have the same order for nodes
rownames(node_data) <- length(t$tip.label) + 1:t$Nnode
node_data[is.na(node_data)] <- 0

# node_data <- node_data[!is.na(rowSums(node_data)),]
node_data$node <- as.integer(rownames(node_data))
colnames(node_data) <- gsub(' ', '_', colnames(node_data))
myColors <- c('red', 'blue', 'yellow')
names(myColors) <- colnames(node_data)[1:3]

# pies <- nodepie(node_data, cols = 1:3, color = myColors)
pies <- nodepie(node_data, cols = 1:3)
pies <- lapply(pies, function(g) g + scale_fill_manual(values = myColors))

p_with_pies <- p_heatmap +
    geom_inset(pies, width = .1, height = .1) 

# df <- tibble::tibble(node = as.numeric(node_data$node), pies = pies) # A tibble of pies (insets)
# p_with_pies <- p_heatmap %<+% 
#     df +
#     ggpp::geom_plot(
#         data = td_filter(!isTip), 
#         mapping = aes(x = x, y = y, label = pies),
#         vp.width = 0.01, vp.height = 0.01, 
#         hjust = 0.5, vjust = 0.5,
#         show.legend = TRUE
# )

# Save plot ---------------------------------------------------------------
ggsave(
    filename = 'test_tree_plot.png', plot = p_with_pies, width = 10, height = 10,
    units = 'in', dpi = 600
)
