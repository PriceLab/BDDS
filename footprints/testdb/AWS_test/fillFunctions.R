simpleLoopFill <- function(dbConnection, tbl){


    print("Beginning simple for loop")
    # Loop through the table, construct each query,
    # and run it in series
    for(i in 1:length(tbl$loc)){

        my.query <- sprintf("insert into regions values ('%s','%s',%d,%d) on conflict (loc) do nothing;",
	                    tbl$loc[i],tbl$chrom[i],tbl$start[i], tbl$end[i])
	dbSendQuery(dbConnection,my.query)
    }

    print("Simple for loop completed")

} #simpleLoopFill
#---------------------------------------------------------------------------
tempTableFill <- function(dbConnection, tbl,ID){

    # Create the temporary table using the passed ID
    my.query <- sprintf("create table %s (like regions);",ID)
    dbSendQuery(dbConnection, my.query)
    
    # Fill the table all at once
    dbWriteTable(dbConnection, ID, tbl, row.names = FALSE, append=TRUE)

    # Use the temp table to fill the real one, then delete it
    my.query <- sprintf("insert into regions select * from %s on conflict (loc) do nothing;",
                        ID)
    dbSendQuery(dbConnection, my.query)
    dbSendQuery(dbConnection, sprintf("drop table %s;",ID))
    
} #tempTableFill
#---------------------------------------------------------------------------
