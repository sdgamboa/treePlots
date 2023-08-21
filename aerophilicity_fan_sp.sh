#! /bin/bash

## graphlan conda environment must be active

sed -i -e "s/\tNA.*$//" aerophilicity_fan_sp.annot

graphlan_annotate.py --annot aerophilicity_fan_sp.annot tree_sp.newick aerophilicity_fan_sp.xml
graphlan.py aerophilicity_fan_sp.xml aerophilicity_fan_sp.png --dpi 300 --size 25
