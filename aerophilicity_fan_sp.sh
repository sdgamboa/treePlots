#! /bin/bash

## graphlan conda environment must be active

sed -i -e "s/\tNA.*$//" aerophilicity_fan_sp.annot

graphlan_annotate.py --annot general.annot tree_sp.newick general_sp_tree_fan.xml
graphlan_annotate.py --annot aerophilicity_fan_sp.annot general_sp_tree_fan.xml aerophilicity_fan_sp.xml
graphlan.py aerophilicity_fan_sp.xml aerophilicity_fan_sp.png --dpi 300 --size 25
