# ------------------------------------------------------------------------------
# main program.R
#
# analysis code to estimate vaccine averted AMR health burden
# ------------------------------------------------------------------------------

# load libraries
library (countrycode)         # country codes
library (data.table)          # data table

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

# ------------------------------------------------------------------------------

# create data table of (194) WHO countries classified by (6) WHO regions 
# and (4) World Bank income levels
regional_classification_file       <- file.path ("data", "regional_classification.csv")
country_income_classification_file <- file.path ("data", "country_income_classification.csv")
create_country_table (regional_classification_file,
                      country_income_classification_file)

# ------------------------------------------------------------------------------

# return to source directory
setwd (source_wd)

# end time
end_time <- Sys.time ()
print (paste0 ("end time = ", end_time))
print (Sys.time () - start_time)
# ------------------------------------------------------------------------------
