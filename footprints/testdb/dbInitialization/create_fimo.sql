/* 
this is legacy code used to set up using the trena role/table instead of the trenatest role

I'm not changing it to the trenatest table because fimo is big, so I'm using the current database (on whovian)
*/

\connect trena;
drop table fimo_hg38;
create table fimo_hg38(motifname varchar,
                       chrom varchar,
                       start int,
                       endpos int,
                       strand varchar,
                       motifscore float,
                       pval float,
                       empty char(1),
                       sequence varchar
                       );
grant all on table fimo_hg38 to trena;
