\connect piq;
\copy hits from 'hits.tsv' delimiter E'\t' CSV NULL as 'NULL';
