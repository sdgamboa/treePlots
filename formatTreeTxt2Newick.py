from ete3 import Tree
t = Tree('sp_tree.txt', format = 1)
t.write(format = 1, outfile = 'sp_tree.newick')
