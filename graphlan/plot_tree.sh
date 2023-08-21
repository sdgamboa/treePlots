#! /bin/bash

# Basic script for annotating a tree and generating an image
graphlan_annotate.py --annot myAnnot.txt myTree.newick myTree.xml
graphlan.py myTree.xml myTree.png --dpi 300 --size 25 
