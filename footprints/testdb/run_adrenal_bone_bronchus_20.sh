#/bin/bash

# Run adrenal_gland first
cd /scratch/github/BDDS/footprints/testdb/adrenal_gland_20
./run_adrenal_gland_20.sh 

# Run bone_element second
cd /scratch/github/BDDS/footprints/testdb/bone_element_20
./run_bone_element_20.sh 


# Run bronchus last
cd /scratch/github/BDDS/footprints/testdb/bronchus_20
./run_bronchus_20.sh
