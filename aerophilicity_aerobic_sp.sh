#! /bin/bash

## graphlan conda environment must be active

sed -i -e "s/\tNA.*$//" aerophilicity_aerobic_sp.annot

graphlan_annotate.py --annot aerophilicity_aerobic_sp.annot tree_sp.newick aerophilicity_aerobic_sp.xml
graphlan.py aerophilicity_aerobic_sp.xml aerophilicity_aerobic_sp.png --dpi 300 --size 25
