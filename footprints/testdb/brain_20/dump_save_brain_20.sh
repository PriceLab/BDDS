#/bin/bash

cd /scratch/db

pg_dump -Fc -h localhost -U postgres brain_hint_20 > ./brain_hint_20.dump &
pg_dump -Fc -h localhost -U postgres brain_wellington_20 > ./brain_wellington_20.dump

wait

aws s3 cp ./brain_hint_20.dump s3://marichards/completed_dbs/
aws s3 cp ./brain_wellington_20.dump s3://marichards/completed_dbs/
