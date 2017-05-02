#/bin/bash

# Run brain first
cd /scratch/github/BDDS/footprints/testdb/brain_16
./run_brain_16.sh 

# Run skin second
cd /scratch/github/BDDS/footprints/testdb/skin_16
./run_skin_16.sh 


# Run lymphoblasts last
cd /scratch/github/BDDS/footprints/testdb/lymphoblast_16
./run_lymphoblast_16.sh
