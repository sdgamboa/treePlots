library(phytools)
library(taxPPro)
library(data.tree)
library(readr)
library(dplyr)
library(tidyr)

tree <- read.tree('tree_sp.newick')

data <- read_csv('gram_stain_prediction_sp.csv', show_col_types = FALSE)
data <- data |> 
    # filter(!grepl('^t__', NCBI_ID)) |> 
    # filter(!grepl('^s__', NCBI_ID)) |> 
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
    relocate(node, .after = last_col()) 
rownames(data) <- data$node
data$node <- NULL

new_data <- data[tree$node.label,]
new_row_names <- Ntip(tree) + 1:tree$Nnode
rownames(new_data) <- new_row_names
new_data <- new_data[which(!is.na(rowSums(new_data))),]

cols <- c(
    'red', 'blue', 'yellow'
)

tree$tip.label <- NULL
plotTree(tree, type = 'fan', fsize = 0.01, ftype = 'i')
nodelabels(
    node = rownames(new_data),
    pie = new_data, piecol = cols, cex  = 0.2
)
