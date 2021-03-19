#!/bin/bash

input=${1?Error: no name given}
name=($(echo $input | tr "." "\n"))

# HYDROGENs
echo "Hydrogen addition in $input"
# Point out the path to phenix
source ~/phenix/phenix-1.18.2-3874/phenix_env.sh
phenix.reduce -NOFLIP -Quiet $input > tmp.pdb  
# > /dev/null 2>&1
mv tmp.pdb $name'_new.pdb'

# WATER PROTONATION
echo "Water protonation in $input"
if grep -q "HOH" $name'_new.pdb' # if any water exists in file
then
sed s/"AHOH"/" HOH"/g $name'_new.pdb' | sed /"BHOH"/d > $name'_tmp.pdb'
phenix.ready_set add_h_to_water=True output_file_name=$name"_H" $name'_tmp.pdb' 
mv $name"_H.pdb" $name'.pdb'
fi

# rm $name'_new.pdb'
# rm $name'_tmp.pdb'
# rm *.eff 

# SECONDARY STRUCTURE
echo "Making PHI file"
stride $name'.pdb' > $name'.phi'
