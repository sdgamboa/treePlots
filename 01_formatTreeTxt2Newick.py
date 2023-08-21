
## activate coda environment with ete3 for this
from ete3 import Tree
t = Tree('tree_sp.txt', format = 1)
t.write(format = 1, outfile = 'tree_sp.newick')
