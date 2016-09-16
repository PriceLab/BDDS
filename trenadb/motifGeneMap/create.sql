\connect trena;
drop table tfmotifs;
create table tfmotifs(motif varchar,
                      gene varchar);
grant all on table tfmotifs to trena;
