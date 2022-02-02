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
source ("function.ua.R")

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

# create combined table: vaccine profile + disease burden
# create separate burden file for health burdens attributable to AMR and associated with AMR  

combined_dt <- create_combined_table(death_burden_dt    = death_burden_dt, 
                         vaccine_profile_dt = vaccine_profile_dt)

# ------------------------------------------------------------------------------
# Death trend by pathogen across all age groups
# the outputs were used to decide the age of vaccination

pathogenlist <- unique(death_burden_dt$Pathogen) 

lapply(pathogenlist, create_death_by_pathogen_graph)

# ------------------------------------------------------------------------------
# estimate vaccine averted AMR deaths & uncertainty analysis

run <- 100

deaths_attributable_psa <- uncertainty_analysis_baseline(psa   = run, 
                                                         data  = read_csv("attributable_burden.csv"))


deaths_associated_psa <- uncertainty_analysis_baseline(psa   = run, 
                                                       data  = read_csv("associated_burden.csv"))

# ------------------------------------------------------------------------------
# vaccine impact by region for deaths attributable to AMR  
impact_by_region_attributable <- data.table(WHO_region     = character(), 
                                            averted_burden = numeric(), 
                                            run_id         = numeric())

for(i in 1:run){
  dt <- estimate_vaccine_impact(i, data=deaths_attributable_psa)
  dt <- dt %>%
    group_by(WHO_region, run_id) %>%
    summarise(averted_burden=sum(va_health_burden), .groups = 'drop')
  impact_by_region_attributable <- rbindlist (list (impact_by_region_attributable, dt),
                                         use.names = TRUE) 
}

# vaccine impact by region for deaths associated to AMR  
impact_by_region_associated <- data.table(WHO_region     = character(), 
                                          averted_burden = numeric(), 
                                          run_id         = numeric())

for(i in 1:run){
  dt <- estimate_vaccine_impact(i, data=deaths_associated_psa)
  dt <-  dt %>%
    group_by(WHO_region, run_id) %>%
    summarise(averted_burden=sum(va_health_burden), .groups = 'drop')
  impact_by_region_associated <- rbindlist (list (impact_by_region_associated, dt),
                                         use.names = TRUE) 
}

# ------------------------------------------------------------------------------
# [table 2] Deaths and DALYs associated with and attributable to bacterial antimicrobial resistance
# globally and by WHO_region, 2019

Attributable_death_averted <- aggregate_impact_by_region(impact_by_region=impact_by_region_attributable)

Associated_death_averted <- aggregate_impact_by_region(impact_by_region=impact_by_region_associated)

# combine into one table
Deaths_Averted <- left_join(Attributable_death_averted, Associated_death_averted, 
                            by=c("Counts" = "Counts"))

fwrite (x    = Deaths_Averted,
        file = "Deaths_Averted.csv")
# ------------------------------------------------------------------------------






# ------------------------------------------------------------------------------
# return to source directory
setwd (source_wd)

# end time
end_time <- Sys.time ()
print (paste0 ("end time = ", end_time))
print (Sys.time () - start_time)
# ------------------------------------------------------------------------------
