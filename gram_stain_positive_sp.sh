#! /bin/bash

## graphlan conda environment must be active

attrname=gram_stain_positive


sed -i -e "s/\tNA.*$//" "$attrname"_sp.annot

graphlan_annotate.py --annot general.annot tree_sp.newick general_sp_tree_"$attrname".xml
graphlan_annotate.py --annot "$attrname"_sp.annot general_sp_tree_"$attrname".xml "$attrname"_sp.xml
graphlan.py "$attrname"_sp.xml "$attrname"_sp.png --dpi 300 --size 25
