#/bin/bash

sudo -u postgres psql skin_wellington_20 <<EOF

create index regions_index on regions (loc, start, endpos);
create index hits_index on hits (loc);

\connect skin_hint_20

create index regions_index on regions (loc, start, endpos);
create index hits_index on hits (loc);

EOF 
