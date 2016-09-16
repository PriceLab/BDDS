\connect piq;
\copy regions from 'regions.tsv' delimiter E'\t' CSV NULL as 'NULL';
