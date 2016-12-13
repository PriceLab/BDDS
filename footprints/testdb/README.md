The *testdb* directory includes the processing code to intersect hint or wellington ouptut with the FIMO database and put the results in a database. This code is in the `src` directory.

The code is based on examples in the `legacy_examples` directory, which can be considered legacy code, and aren't actively used any more.

Currently, my workflow is:
- [ ] create a database (using code in the `dbInitialization` folder)
- [ ] make a new folder for running each batch of data
- [ ] copy the `hint.R` and `wellington.R` master scripts to the newly created folder from a previous run (such as the `lymphoblast` folder)
- [ ] process the data using the master scripts (in R) to fill the database
- [ ] index the database 
- [ ] make database read only

I intend to write more details about that workflow below.
