# ------------------------------------------------------------------------------
# main anlaysis.R
#
# analysis code to estimate vaccine averted AMR health burden
# ------------------------------------------------------------------------------

# load libraries
library (readr)
library (readxl)
library (dplyr)
library (ggplot2)
library (reshape2)
library (tidyverse)
library (data.table)
library (ggplot2)
library (rriskDistributions)

# remove all objects from workspace
rm (list = ls ())

# start time
start_time <- Sys.time ()
print (paste0 ("start time = ", start_time))

# source functions
source ("function.R")

# move to base directory (run code from source directory)
source_wd <- getwd ()
setwd("~/GitHub/vaccine_amr/data")

# ------------------------------------------------------------------------------
## IHME data on the AMR burden ##
# WHO regions: Africa, Americas, Eastern Mediterranean, Europe, South-East Asia, 
#             Western Pacific, Unclassified 
# year: 2019
# sex: both
# ------------------------------------------------------------------------------
# create data table of AMR burden (deaths) classified by pathogen, 
# and disease presentation, age groups

AMR_death_burden_file <- read_excel("IHME AMR burden.xlsx", 
                              col_names = FALSE)             

death_burden_dt <- create_death_burden_table(AMR_death_burden_file)

# create data table of vaccine profile
vaccine_profile_file <- read_excel("Vaccine profile for IHME burden.xlsx", 
                              sheet = "Vaccine profile assumptions")

vaccine_profile_dt <- create_vaccine_profile_table(vaccine_profile_file)

# ------------------------------------------------------------------------------
# Death trend by pathogen across all age groups
# the outputs were used to decide the age of vaccination

pathogenlist <- unique(death_burden_dt$Pathogen) 

lapply(pathogenlist, create_death_by_pathogen_graph)

# ------------------------------------------------------------------------------
# estimate vaccine averted AMR deaths

BurdenAverted_attributable_dt <- estimate_vaccine_impact_a()

# estimate vaccine avertable associated AMR deaths
BurdenAverted_associated_dt <- estimate_vaccine_impact_b()


# ------------------------------------------------------------------------------




# ------------------------------------------------------------------------------






# ------------------------------------------------------------------------------
# return to source directory
setwd (source_wd)

# end time
end_time <- Sys.time ()
print (paste0 ("end time = ", end_time))
print (Sys.time () - start_time)
# ------------------------------------------------------------------------------
