# ------------------------------------------------------------------------------
# main program.R
#
# analysis code to estimate vaccine averted AMR health burden
# ------------------------------------------------------------------------------

# remove all objects from workspace
rm (list = ls ())

# start time
start_time <- Sys.time ()
print (paste0 ("start time = ", start_time))

# source functions
source ("functions.R")

# move to base directory (run code from source directory)
source_wd <- getwd ()
setwd ("../")

# # create data table of (194) WHO countries classified by (6) WHO regions 
# and (4) World Bank income levels
create_country_table ()

# return to source directory
setwd (source_wd)


# end time
end_time <- Sys.time ()
print (paste0 ("end time = ", start_time))
print (Sys.time () - start_time)
# ------------------------------------------------------------------------------
