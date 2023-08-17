
Generate a simple tree


I need to export to txt in R and then export to newick from within python.
Otherwise, the file cannot be processed correctly by GraPhlan.


```
graphlan.py sp_tree.newick newtree.png --size 100 # it seems that 100 is the limit for this tree (I get errors otherwise).
```

`--size 25` seems to be a good option

```
graphlan_annotate.py --annot sp_tree_annotation.txt sp_tree.newick sp_tree.xml
graphlan.py sp_tree.xml sp_tree.png --size 25
```


```
scp samuelgamboa@super:/mnt/STORE1/bighome/samuelgamboa/Projects/taxPProValidation/remove_auc.sh .
```

Check path with

```
readlink -f <filename>
```
