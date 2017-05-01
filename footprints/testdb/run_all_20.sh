#/bin/bash

# Run brain first
cd /scratch/github/BDDS/footprints/testdb/brain_20
./run_brain_20.sh &

wait

# Run skin second
cd /scratch/github/BDDS/footprints/testdb/skin_20
./run_skin_20.sh &

wait

# Run lymphoblasts last
cd /scratch/github/BDDS/footprints/testdb/lymphoblast_20
./run_lymphoblast_20.sh &
