#! /bin/bash

# Basic script for annotating a tree and generating an image
graphlan_annotate.py --annot example_tree.annot example_tree.newick example_tree.xml
graphlan.py example_tree.xml example_tree.png --dpi 300 --size 25 
