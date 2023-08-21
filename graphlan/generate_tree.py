from ete3 import NCBITaxa
ncbi = NCBITaxa()
tree = ncbi.get_topology([9606, 9598, 100090, 7707, 8782])
print(tree)
tree.write(outfile = 'myTree.newick')
