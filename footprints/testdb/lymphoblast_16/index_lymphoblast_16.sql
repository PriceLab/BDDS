\connect lymphoblast_wellington_16;
create index regions_index on regions (loc, start, endpos);
create index hits_index on hits (loc);

\connect lymphoblast_hint_16;
create index regions_index on regions (loc, start, endpos);
create index hits_index on hits (loc);
