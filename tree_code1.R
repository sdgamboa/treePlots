library(ggtree)
library(readr)
library(ape)
library(dplyr)
library(tidyr)
t <- read.tree('tree_sp.newick')
tp <- ggtree(t, layout = 'circular')
data <- read_csv('gram_stain_prediction_sp.csv', show_col_types = FALSE)
data <- data |> 
    filter(!grepl('^t__', NCBI_ID)) |> 
    select(NCBI_ID, Attribute, Score) |> 
    distinct() |> 
    pivot_wider(
        names_from = 'Attribute', values_from = 'Score'
    ) |> 
    as.data.frame() |> 
    rename(node = NCBI_ID) |> 
    relocate(node, .after = last_col())
myPies <- nodepie(data, cols = 1:3)
myPlot <- inset(tp, myPies, width = 0.03, height = 0.03)
ggsave(
    filename = 'test_tree_plot.png', plot = myPlot, width = 20, height = 20, 
    units = 'in'
)
