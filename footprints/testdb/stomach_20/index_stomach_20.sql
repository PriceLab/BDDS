\connect stomach_wellington_20;
create index regions_index on regions (loc, start, endpos);
create index hits_index on hits (loc);

\connect stomach_hint_20
create index regions_index on regions (loc, start, endpos);
create index hits_index on hits (loc);
