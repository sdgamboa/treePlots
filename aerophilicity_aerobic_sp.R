
## In this script everything will be at the species level
library(dplyr)

tree <- ape::read.tree(file = 'tree_sp.newick')
nodes <-  unique(c(tree$node.label, tree$tip.label))

holdout <- read.csv('aerophilicity_holdout_sp.csv')
predicted <- read.csv('aerophilicity_prediction_sp.csv') |> 
    mutate(NCBI_ID = paste0(sub('^(\\w).*$', '\\1__', Rank), NCBI_ID)) |> 
    filter(Frequency != 'never')

ev1 <- c('igc', 'exp', 'nas', 'tas')
ev2 <- c('asr', 'inh')
    
## aerobic ####
l <- list(
    aerobic_holdout_ring_color = holdout |> 
        filter(Attribute == 'aerobic') |> 
        select(NCBI_ID) |> 
        mutate(
            x = 'ring_color',
            y = '2',
            z = 'r'
        ),
    aerobic_training_node_color = predicted |> 
        filter(Attribute == 'aerobic') |> 
        filter(Evidence %in% ev1) |> 
        select(NCBI_ID) |> 
        mutate(
            x = 'clade_marker_color',
            y = 'b'
        ),
    aerobic_predicted_node_color = predicted |> 
        filter(Attribute == 'aerobic') |> 
        filter(Evidence %in% ev2) |> 
        select(NCBI_ID) |> 
        mutate(
            x = 'clade_marker_color',
            y = 'y'
        ),
    ## ring
    aerobic_predicted_holdout_ring_color <-  predicted |> 
        filter(Attribute == 'aerobic') |> 
        filter(Evidence %in% ev2) |> # propagation 
        filter(NCBI_ID %in% unique(holdout$NCBI_ID)) |> 
        select(NCBI_ID) |> 
        mutate(
            x = 'ring_color', # option name
            y = '1', # ring level
            z = 'g' # ring color
        )
)






df <- bind_rows(l)
fname <- 'aerophilicity_aerobic_sp.annot'
write.table(
    x = df, file = fname, quote = FALSE, col.names = FALSE, 
    row.names = FALSE, sep = '\t'
)


