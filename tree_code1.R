library(ggtree)
library(readr)
library(ape)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpp)
library(deeptime)


t <- read.tree('tree_sp.newick')
# tp <- ggtree(t, layout = 'circular')

data <- read_csv('gram_stain_prediction_sp.csv', show_col_types = FALSE)
data <- data |> 
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
    filter(node %in% t$node.label)

rownames(data) <- data$node
data$node <- NULL
data <- data[t$node.label,]
rownames(data) <- length(t$tip.label) + 1:t$Nnode
data <- data[!is.na(rowSums(data)),]
data$node <- as.integer(rownames(data))

myPies <- nodepie(data, cols = 1:3)
# t2 <- inset(tree_view = tp, insets = myPies[1:3], width = 0.03, height = 0.03)

df <- tibble::tibble(node=as.numeric(data$node), pies=myPies)
p1 <- ggtree(t, layout = 'circular', size = 0.025)
p2 <- p1 %<+% df
p3 <- p2 + 
    geom_plot(
    data = td_filter(!isTip), 
    mapping = aes(x = x, y = y, label = pies), 
    vp.width = 0.01, vp.height = 0.01, 
    hjust = 0.5, vjust = 0.5
) 
ggsave(
    filename = 'test_tree_plot.png', plot = p3, width = 10, height = 10,
    units = 'in', dpi = 300
)
