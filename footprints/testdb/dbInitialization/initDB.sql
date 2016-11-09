--
-- A script to set up a new postgres database for footprinting work.
--
-- This script should be run with a superuser role to initialize a fresh postgres install
--
-- psql -U postgres -h localhost -f initDB.sql
--

-- create non-superuser for making new databases
CREATE USER databasemaker password 'databasemaker';
CREATE DATABASE databasemaker;
GRANT ALL PRIVILEGES on DATABASE databasemaker to databasemaker;
ALTER ROLE databasemaker CREATEROLE CREATEDB LOGIN;

-- switch to the databasemaker role
SET ROLE databasemaker;

-- create the trenatest user and grant login
CREATE USER trenatest password 'trenatest';
CREATE DATABASE trenatest;
GRANT ALL PRIVILEGES on DATABASE trenatest to trenatest;
ALTER ROLE trenatest LOGIN;
GRANT trenatest TO databasemaker;

-- now the database is ready for the trenatest user to create databases and run the create.sql scripts (except for FIMO)

-- create the trena user and grant login
CREATE USER trena password 'trena';
CREATE DATABASE trena;
GRANT ALL PRIVILEGES on DATABASE trena to trena;
ALTER ROLE trena LOGIN;
GRANT trena to databasemaker;

-- now the trena user can run createFimo.sql
