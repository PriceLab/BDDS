\connect wholeBrain-wellington;
drop index hits_index;
create index hits_index on hits (loc);

