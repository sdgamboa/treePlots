
# setwd('/home/user/Projects/CUNY/treePlots')

library(data.tree)
library(taxPPro)
data("tree_list")
data_tree <- as.Node(tree_list)
full_newick <- ToNewick(data_tree, heightAttribute = NULL)
fileConn <- file('full_tree.txt')
writeLines(full_newick, fileConn)
close(fileConn)

Prune(data_tree, function(x) !grepl("t__", x$name))
sp_newick <- ToNewick(data_tree, heightAttribute = NULL)
fileConn <- file('sp_tree.txt')
writeLines(sp_newick, fileConn)
close(fileConn)


