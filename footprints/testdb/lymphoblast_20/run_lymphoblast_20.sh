#/bin/bash

cd /scratch/db
sudo -u postgres psql postgres << EOF

create database lymphoblast_wellington_20;
grant all privileges on database lymphoblast_wellington_20 to trena;
create database lymphoblast_hint_20;
grant all privileges on database lymphoblast_hint_20 to trena;

\connect lymphoblast_wellington_20

create table regions(loc varchar primary key,
		     chrom varchar,
		     start int,
		     endpos int);

grant all on table "regions" to trena;

create table hits(loc varchar,
                  fp_start int,
                  fp_end int,
		  type varchar,
		  name varchar,
		  length int,
		  strand char(1),
		  sample_id varchar,
		  method varchar,
		  provenance varchar,
		  score1 real,
		  score2 real,
		  score3 real,
		  score4 real,
		  score5 real,
		  score6 real);

grant all on table "hits" to trena;

\connect lymphoblast_hint_20

create table regions(loc varchar primary key,
		     chrom varchar,
		     start int,
		     endpos int);

grant all on table "regions" to trena;

create table hits(loc varchar,
                  fp_start int,
                  fp_end int,
		  type varchar,
		  name varchar,
		  length int,
		  strand char(1),
		  sample_id varchar,
		  method varchar,
		  provenance varchar,
		  score1 real,
		  score2 real,
		  score3 real,
		  score4 real,
		  score5 real,
		  score6 real);

grant all on table "hits" to trena;

EOF

cd /scratch/github/BDDS/footprints/testdb/lymphoblast_20

R -f hint.R &
R -f wellington.R &
