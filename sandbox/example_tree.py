from ete3 import Tree
t = Tree('example_tree.txt', format = 1)
t.write(format = 1, outfile = 'example_tree.newick')
