
library(data.tree)
library(taxPPro)

data("tree_list")
data_tree <- as.Node(tree_list)
Prune(data_tree, function(x) !grepl("t__", x$name))
sp_newick <- ToNewick(data_tree, heightAttribute = NULL)

if (!interactive()) {
    fname_sp <- file.path(getwd(), 'tree_sp.txt')
} else {
    fname_sp <- 'tree_sp.txt'
}

fileConn <- file(fname_sp)
writeLines(sp_newick, fileConn)
close(fileConn)


