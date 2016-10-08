\connect wholeBrain-wellington;
drop index regions_index;
create index regions_index on regions (loc, chrom, start, endpos);


