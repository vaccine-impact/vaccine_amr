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
library (formattable)

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
## IHME data on the AMR burden ##
# WHO regions: Africa, Americas, Eastern Mediterranean, Europe, South-East Asia, 
#             Western Pacific, Unclassified 
# year: 2019
# sex: both
# ------------------------------------------------------------------------------
# create data table of AMR burden (deaths) classified by pathogen, 
# and disease presentation, age groups

AMR_death_burden_file <- read_excel(file.path("data", "IHME AMR burden.xlsx"),
                                    col_names = FALSE)             

death_burden_dt <- create_death_burden_table(AMR_death_burden = AMR_death_burden_file,
                                             death_burden_file = file.path ("tables", "AMR_death_burden.csv"))

# create data table of vaccine profile
vaccine_profile_file <- read_excel(file.path("data", "Vaccine profile for IHME burden.xlsx"), 
                                   sheet = "Vaccine profile assumptions")

vaccine_profile_dt <- create_vaccine_profile_table(vaccine_profile = vaccine_profile_file,
                                                   vaccine_profile_file = file.path("tables", "vaccine_profile.csv"))

# create combined table: vaccine profile + disease burden
# create separate burden file for health burdens attributable to AMR and associated with AMR  

combined_dt <- create_combined_table(death_burden_dt          = death_burden_dt, 
                                     vaccine_profile_dt       = vaccine_profile_dt,
                                     attributable_burden_file = file.path("tables", "attributable_burden.csv"),
                                     associated_burden_file   = file.path("tables", "associated_burden.csv"))

# ------------------------------------------------------------------------------
# Death trend by pathogen across all age groups
# the outputs were used to decide the age of vaccination

pathogenlist <- unique(death_burden_dt$Pathogen)

lapply(pathogenlist, create_death_by_pathogen_graph)

# ------------------------------------------------------------------------------
# estimate vaccine averted AMR deaths & uncertainty analysis

set.seed (3)  # seed for random number generator
run <- 200 # number of runs for probabilistic sensitivity analysis

deaths_attributable_psa <- uncertainty_analysis_baseline(psa   = run, 
                                                         data  = read_csv(file.path("tables", "attributable_burden.csv")))


deaths_associated_psa <- uncertainty_analysis_baseline(psa   = run, 
                                                       data  = read_csv(file.path("tables", "associated_burden.csv")))

# ------------------------------------------------------------------------------
# Deaths and DALYs associated with and attributable to 
# AMR globally and by WHO region, 2019

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

Attributable_death_averted <- aggregate_impact_by_region(impact_by_region=impact_by_region_attributable)

Associated_death_averted   <- aggregate_impact_by_region(impact_by_region=impact_by_region_associated)

# combine into one table
Death_Averted <- left_join(Associated_death_averted, Attributable_death_averted,
                           by=c("Counts" = "Counts"))

# ------------------------------------------------------------------------------
# [table 2] Deaths and DALYs associated with and attributable to 
# bacterial antimicrobial resistance globally and by WHO region, 2019

Death_Averted$`2.5%.x`  <- comma(Death_Averted$`2.5%.x`,  format = "d")
Death_Averted$`50%.x`   <- comma(Death_Averted$`50%.x`,   format = "d")
Death_Averted$`97.5%.x` <- comma(Death_Averted$`97.5%.x`, format = "d")
Death_Averted$`2.5%.y`  <- comma(Death_Averted$`2.5%.y`,  format = "d")
Death_Averted$`50%.y`   <- comma(Death_Averted$`50%.y`,   format = "d")
Death_Averted$`97.5%.y` <- comma(Death_Averted$`97.5%.y`, format = "d")

Death_Averted [,"Associated with resistance":= paste(Death_Averted$"50%.x","(",Death_Averted$"2.5%.x","-",Death_Averted$"97.5%.x",")")]
                                                     
Death_Averted [,"Attributable to resistance":= paste(Death_Averted$"50%.y","(",Death_Averted$"2.5%.y","-",Death_Averted$"97.5%.y",")")]

Death_Averted <- Death_Averted [, c("Counts", "Associated with resistance", "Attributable to resistance")]

Death_Averted <- Death_Averted [c(2,3,4,5,6,7,1),]

fwrite (x    = Death_Averted,
        file = file.path ("tables", "Deaths_Averted.csv"))
# ------------------------------------------------------------------------------
# Figure 1: vaccine avertable attributable to and associated with 
# bacterial antimicrobial resistance by GBD  region, 2019

create_burden_averted_by_region_graph(Attributable_burden_averted = Attributable_death_averted,
                                      Associated_burden_averted   = Associated_death_averted)

# ------------------------------------------------------------------------------
# Figure 2: Global vaccine avertable deaths (counts) attributable to and associated with 
# bacterial antimicrobial resistance by infectious syndrome, 2019

# create table for avertable deaths attributable to AMR by disease presentation
impact_by_dp_attributable <- data.table(Disease_presentation = character(), 
                                        averted_burden       = numeric(), 
                                        run_id               = numeric())

for(i in 1:run){
  dt <- estimate_vaccine_impact(i, data=deaths_attributable_psa)
  dt <- dt %>%
    group_by(Disease_presentation, run_id) %>%
    summarise(averted_burden = sum(va_health_burden), .groups = 'drop')
    impact_by_dp_attributable <- rbindlist (list (impact_by_dp_attributable, dt),
                                              use.names = TRUE) 
}

# create table for avertable deaths associated with AMR by disease presentation 
impact_by_dp_associated <- data.table(Disease_presentation = character(), 
                                      averted_burden       = numeric(), 
                                      run_id               = numeric())

for(i in 1:run){
  dt <- estimate_vaccine_impact(i, data=deaths_associated_psa)
  dt <-  dt %>%
    group_by(Disease_presentation, run_id) %>%
    summarise(averted_burden=sum(va_health_burden), .groups = 'drop')
    impact_by_dp_associated <- rbindlist (list (impact_by_dp_associated, dt),
                                            use.names = TRUE)
}

# aggregate the data
Attributable_death_averted_dp <- aggregate_impact_by_dp(impact_by_dp=impact_by_dp_attributable)

Associated_death_averted_dp <- aggregate_impact_by_dp(impact_by_dp=impact_by_dp_associated)

# create Figure 2
create_burden_averted_by_dp_graph(Attributable_burden_averted = Attributable_death_averted_dp,
                                  Associated_burden_averted   = Associated_death_averted_dp)

# ------------------------------------------------------------------------------
# Figure 3: Global vaccine avertable deaths (counts) attributable to and associated with 
# bacterial antimicrobial resistance by pathogen, 2019

# create table for avertable deaths attributable to AMR by pathogen
impact_by_pathogen_attributable <- data.table(Pathogen             = character(), 
                                              averted_burden       = numeric(), 
                                              run_id               = numeric())

for(i in 1:run){
  dt <- estimate_vaccine_impact(i, data=deaths_attributable_psa)
  dt <- dt %>%
    group_by(Pathogen, run_id) %>%
    summarise(averted_burden = sum(va_health_burden), .groups = 'drop')
  impact_by_pathogen_attributable <- rbindlist (list (impact_by_pathogen_attributable, dt),
                                                use.names = TRUE) 
}

# create table for avertable deaths associated with AMR by pathogen
impact_by_pathogen_associated <- data.table(Pathogen         = character(), 
                                            averted_burden     = numeric(), 
                                            run_id           = numeric())

for(i in 1:run){
  dt <- estimate_vaccine_impact(i, data=deaths_associated_psa)
  dt <-  dt %>%
    group_by(Pathogen, run_id) %>%
    summarise(averted_burden=sum(va_health_burden), .groups = 'drop')
  impact_by_pathogen_associated <- rbindlist (list (impact_by_pathogen_associated, dt),
                                              use.names = TRUE) 
}

# aggregate the data
Attributable_death_averted_pathogen <- aggregate_impact_by_pathogen(impact_by_pathogen=impact_by_pathogen_attributable)

Associated_death_averted_pathogen <- aggregate_impact_by_pathogen(impact_by_pathogen=impact_by_pathogen_associated)

# create Figure 3
create_burden_averted_by_pathogen_graph(Attributable_burden_averted = Attributable_death_averted_pathogen,
                                        Associated_burden_averted   = Associated_death_averted_pathogen)
# ------------------------------------------------------------------------------
# [table in appendix] vaccine avertable health burdens associated with and attributable
# to AMR by WHO region, pathogen, disease presentation, and age group
AMR_burden_data_updated <- update_death_burden(combined_dt     = combined_dt,
                                               death_burden_dt = death_burden_dt,
                                               AMR_burden_data_updated_file = file.path ("tables", "AMR burden data_updated.csv"))

# ------------------------------------------------------------------------------
# return to source directory
setwd (source_wd)

# end time
end_time <- Sys.time ()
print (paste0 ("end time = ", end_time))
print (Sys.time () - start_time)
# ------------------------------------------------------------------------------
