#/bin/bash

sudo -u postgres psql brain_wellington_20 <<EOF

create index regions_index on regions (loc, start, endpos);
create index hits_index on hits (loc);

\connect brain_hint_20

create index regions_index on regions (loc, start, endpos);
create index hits_index on hits (loc);

EOF 
