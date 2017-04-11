# In ~20 of the names there was a space that created some unwanted parsing of the names.
# these are the postgresql commands I used to fix those errors. 
\connect hg38

# identify all the offending motifs
SELECT distinct(motif) from motifsgenes where motif like '%% %%';

UPDATE motifsgenes SET motif = 'gr09.v1_Phd1_deBruijn' WHERE motif = 'gr09.v1_Phd1_deBruijn UP00351';
UPDATE motifsgenes SET motif = 'Nkx2-6_3437.1' WHERE motif = 'UP00147_1 Nkx2-6_3437.1';
UPDATE motifsgenes SET motif = 'Nkx6-3_3446.1' WHERE motif = 'UP00238_1 Nkx6-3_3446.1';
UPDATE motifsgenes SET motif = 'sci09.v2_Plagl1_0972' WHERE motif = 'sci09.v2_Plagl1_0972 UP00088';
UPDATE motifsgenes SET motif = 'Nkx1-1_3856.3' WHERE motif = 'UP00220_1 Nkx1-1_3856.3';
UPDATE motifsgenes SET motif = 'UP00238_1 Nkx6-3_3446.1' WHERE motif = 'UP00238_1 Nkx6-3_3446.1';
UPDATE motifsgenes SET motif = 'sci09.v1_Sp100_2947' WHERE motif = 'sci09.v1_Sp100_2947 UP00049';
UPDATE motifsgenes SET motif = 'sci09.v2_Osr1_3033' WHERE motif = 'sci09.v2_Osr1_3033 UP00027';
UPDATE motifsgenes SET motif = 'sci09.v1_Bbx_3753' WHERE motif = 'sci09.v1_Bbx_3753 UP00012';
UPDATE motifsgenes SET motif = 'sci09.v2_Sp100_2947' WHERE motif = 'sci09.v2_Sp100_2947 UP00049';
UPDATE motifsgenes SET motif = 'sci09.v1_Osr1_3033' WHERE motif = 'sci09.v1_Osr1_3033 UP00027';
UPDATE motifsgenes SET motif = 'Nkx1-2_3214.1' WHERE motif = 'UP00139_1 Nkx1-2_3214.1';
UPDATE motifsgenes SET motif = 'sci09.v2_Zbtb3_1048' WHERE motif = 'sci09.v2_Zbtb3_1048 UP00031';
UPDATE motifsgenes SET motif = 'gr09.v1_Yap1_deBruijn' WHERE motif = 'gr09.v1_Yap1_deBruijn UP00327';
UPDATE motifsgenes SET motif = 'sci09.v1_Nr2f2_2192' WHERE motif = 'sci09.v1_Nr2f2_2192 UP00009';
UPDATE motifsgenes SET motif = 'sci09.v2_Bbx_3753' WHERE motif = 'sci09.v2_Bbx_3753 UP00012';
UPDATE motifsgenes SET motif = 'sci09.v1_Plagl1_0972' WHERE motif = 'sci09.v1_Plagl1_0972 UP00088';
UPDATE motifsgenes SET motif = 'sci09.v1_Zbtb12_2932' WHERE motif = 'sci09.v1_Zbtb12_2932 UP00019';
UPDATE motifsgenes SET motif = 'sci09.v1_Zbtb3_1048' WHERE motif = 'sci09.v1_Zbtb3_1048 UP00031';
UPDATE motifsgenes SET motif = 'sci09.v2_Nr2f2_2192' WHERE motif = 'sci09.v2_Nr2f2_2192 UP00009';
UPDATE motifsgenes SET motif = 'sci09.v2_Zbtb12_2932' WHERE motif = 'sci09.v2_Zbtb12_2932 UP00019';


\connect skin_hint

UPDATE hits SET name = 'Nkx2-6_3437.1' WHERE name = 'UP00147_1';
UPDATE hits SET name = 'Nkx6-3_3446.1' WHERE name = 'UP00238_1';
UPDATE hits SET name = 'Nkx1-1_3856.3' WHERE name = 'UP00220_1';
UPDATE hits SET name = 'Nkx1-2_3214.1' WHERE name = 'UP00139_1';
