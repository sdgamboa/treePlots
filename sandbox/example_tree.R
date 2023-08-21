library(taxPPro)
library(dplyr)
library(data.tree)
data('tree_list')
data_tree <- as.Node(tree_list)
subtree <- data_tree$d__2$p__1224$c__1236$o__91347
Prune(subtree, function(x) !grepl("t__", x$name))
newick <- ToNewick(subtree, heightAttribute = NULL)
fileConn <- file('~/Projects/CUNY/treePlots/sandbox/example_tree.txt')
writeLines(newick, fileConn)
close(fileConn)

node_names <- unname(subtree$Get(function(node) node$name))
df <- data.frame(x = node_names)

## Add annotations for color
df <- df |> 
    mutate(
        y = 'clade_marker_color',
        z = ifelse(grepl('^s__', x), 'r', 'b')
    )

## Add annotations for size

sp <- df |> 
    filter(grepl('^s__', x)) |> 
    pull(x)


n_size1 <- round(length(sp) / 2)
n_size2 <- length(sp) - size1
df2 <- data.frame(
    x = sp,
    y = 'clade_marker_size',
    z = as.character(c(rep(60, n_size1), rep(300, n_size2)))
)


df <- bind_rows(df, df2)

write.table(
    x = df,
    file = '~/Projects/CUNY/treePlots/sandbox/example_tree.annot',
    sep = '\t', quote = FALSE, col.names = FALSE, row.names = FALSE
)







