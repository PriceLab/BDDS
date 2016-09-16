\connect trena;
\copy tfmotifs from './genesMotifs.tsv' delimiter E'\t' CSV header NULL as 'NULL';
