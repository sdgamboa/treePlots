#! /bin/bash

## graphlan conda environment must be active

sed -i -e "s/\tNA.*$//" aerophilicity_aerobic_sp.annot

graphlan_annotate.py --annot general.annot tree_sp.newick general_sp_tree_aer.xml
graphlan_annotate.py --annot aerophilicity_aerobic_sp.annot general_sp_tree_aer.xml aerophilicity_aerobic_sp.xml
graphlan.py aerophilicity_aerobic_sp.xml aerophilicity_aerobic_sp.png --dpi 300 --size 25
