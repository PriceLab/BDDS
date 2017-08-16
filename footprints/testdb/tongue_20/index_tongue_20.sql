\connect tongue_wellington_20;
create index regions_index on regions (loc, start, endpos);
create index hits_index on hits (loc);

\connect tongue_hint_20
create index regions_index on regions (loc, start, endpos);
create index hits_index on hits (loc);
