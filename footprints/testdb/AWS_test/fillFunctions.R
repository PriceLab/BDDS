simpleLoopFill <- function(dbConnection, tbls){

    print("Beginning simple for loop on regions")
    # Loop through the table, construct each query,
    # and run it in series
    for(i in 1:length(tbls$regions$loc)){

        my.query <- sprintf("insert into regions values ('%s','%s',%d,%d) on conflict (loc) do nothing;",
	                    tbls$regions$loc[i],
			    tbls$regions$chrom[i],
			    tbls$regions$start[i],
			    tbls$regions$end[i])
	dbSendQuery(dbConnection,my.query)
    }
    print("Continuing with hits")

    for(i in 1:length(tbls$regions$loc)){

        my.query <- sprintf("insert into regions values ('%s','%s',%d,%d) on conflict (loc) do nothing;",
	                            tbls$regions$loc[i],				    
								                            tbls$regions$start[i],
											                                tbls$regions$end[i])
															        dbSendQuery(dbConnection,my.query)
																    }


    print("Simple for loop completed")

} #simpleLoopFill
#---------------------------------------------------------------------------
tempTableFill <- function(dbConnection, tbl,ID,tableID){

    # Create the temporary table using the passed ID
    my.query <- sprintf("create table %s (like regions);",ID)
    dbSendQuery(dbConnection, my.query)
    
    # Fill the table all at once
    dbWriteTable(dbConnection, ID, tbl, row.names = FALSE, append=TRUE)

    # Use the temp table to fill the real one, then delete it
    my.query <- sprintf("insert into %s select * from %s on conflict (loc) do nothing;",
                        tableID, ID)
    dbSendQuery(dbConnection, my.query)
    dbRemoveTable(dbConnection, name = ID)
    
} #tempTableFill
#---------------------------------------------------------------------------
batchFill <- function(dbConnection, tbl){

    # Create a temporary text file of queries
    




} #batchFill