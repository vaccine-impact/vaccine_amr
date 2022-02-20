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
library (patchwork)

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

# 2019 vaccine coverage for existing vaccines: HIB vaccine, PCV
existing_vaccine_coverage (hib_coverage_file = file.path("tables", "hib coverage.csv"),
                           pcv_coverage_file = file.path("tables", "pcv coverage.csv"))

# ------------------------------------------------------------------------------

# create data table of AMR burden (deaths) classified by pathogen, 
# and disease presentation, age groups

death_burden_dt <- read_excel(file.path("data", "IHME AMR burden.xlsx"),
                                    col_names = FALSE)             

death_burden_dt <- create_burden_table(AMR_burden = death_burden_dt,
                                       burden_file = file.path ("tables", "AMR_death_burden.csv"))


# create data table of AMR burden (DALYs) classified by pathogen, 
# and disease presentation, age groups

daly_burden_dt <- read_excel(file.path("data", "IHME AMR burden_DALYs.xlsx"),
                             col_names = FALSE)

daly_burden_dt <- create_burden_table(AMR_burden = daly_burden_dt,
                                      burden_file = file.path ("tables", "AMR_daly_burden.csv"))

# will add analysis for Neisseria gonorrhoeae later
daly_burden_dt <- daly_burden_dt[!(daly_burden_dt$Pathogen == "Neisseria gonorrhoeae"),]

# ------------------------------------------------------------------------------
# create data table of vaccine profile
vaccine_profile_file <- read_excel(file.path("data", "Vaccine profile for IHME burden.xlsx"), 
                                   sheet = "Vaccine profile assumptions")

vaccine_profile_dt <- create_vaccine_profile_table(vaccine_profile = vaccine_profile_file,
                                                   vaccine_profile_file = file.path("tables", "vaccine_profile.csv"))

# ------------------------------------------------------------------------------
# create combined table: vaccine profile + disease burden (deaths)
# create separate burden file for health burdens attributable to AMR and associated with AMR  

combined_dt <- create_combined_table(death_burden_dt          = death_burden_dt, 
                                     vaccine_profile_dt       = vaccine_profile_dt,
                                     attributable_burden_file = file.path("tables", "attributable_burden.csv"),
                                     associated_burden_file   = file.path("tables", "associated_burden.csv"))

# create combined table: vaccine profile + disease burden (DALYs)
# create separate burden file for health burdens attributable to AMR and associated with AMR  
daly_combined_dt <- create_combined_table(death_burden_dt     = daly_burden_dt, 
                                     vaccine_profile_dt       = vaccine_profile_dt,
                                     attributable_burden_file = file.path("tables", "daly_attributable_burden.csv"),
                                     associated_burden_file   = file.path("tables", "daly_associated_burden.csv"))
# ------------------------------------------------------------------------------



# ------------------------------------------------------------------------------
# Death trend by pathogen across all age groups
# the outputs were used to decide the age of vaccination

pathogenlist <- unique(death_burden_dt$Pathogen)

lapply(pathogenlist, create_death_by_pathogen_graph)

# ------------------------------------------------------------------------------



# ------------------------------------------------------------------------------
# estimate vaccine averted AMR deaths & uncertainty analysis

set.seed (3)  # seed for random number generator
run <- 200 # number of runs for probabilistic sensitivity analysis

# Deaths
deaths_attributable_psa <- uncertainty_analysis_baseline(psa       = run, 
                                                         tolerance = 0.0016,
                                                         data      = read_csv(file.path("tables", "attributable_burden.csv")))

deaths_associated_psa <- uncertainty_analysis_baseline(psa   = run, 
                                                       tolerance = 0.0016,
                                                       data  = read_csv(file.path("tables", "associated_burden.csv")))

# DALYs
daly_attributable_psa <- uncertainty_analysis_baseline(psa   = run, 
                                                       tolerance = 0.028,
                                                       data  = read_csv(file.path("tables", "daly_attributable_burden.csv")))

daly_associated_psa <- uncertainty_analysis_baseline(psa   = run, 
                                                     tolerance = 0.028,
                                                     data  = read_csv(file.path("tables", "daly_associated_burden.csv")))
# ------------------------------------------------------------------------------
# Deaths and DALYs associated with and attributable to 
# AMR globally and by WHO region, 2019

# vaccine impact on deaths by region

Attributable_death_averted <- aggregate_impact_by_region(input_data=deaths_attributable_psa)

Associated_death_averted   <- aggregate_impact_by_region(input_data=deaths_associated_psa)

Death_Averted <- left_join(Associated_death_averted, Attributable_death_averted,
                           by=c("Counts" = "Counts"))

# vaccine impact on DALYs by region

Attributable_daly_averted <- aggregate_impact_by_region(input_data = daly_attributable_psa)

Associated_daly_averted <- aggregate_impact_by_region(input_data = daly_associated_psa)

DALY_Averted <- left_join(Associated_daly_averted, Attributable_daly_averted,
                          by=c("Counts" = "Counts"))

# ------------------------------------------------------------------------------
# [table 2] Vaccine avertable deaths and DALYs associated with and attributable to 
# bacterial antimicrobial resistance globally and by WHO region, 2019

Death_Averted <- estimate_bruden_averted_by_region(input_data = Death_Averted)

DALY_Averted <- estimate_bruden_averted_by_region(input_data = DALY_Averted)

Table2_avertible_burden_by_region <- left_join(Death_Averted, DALY_Averted,
                                               by=c("Counts" = "Counts"))

fwrite (x    = Table2_avertible_burden_by_region, 
        file = file.path("tables", "Table2_avertible_burden_by_region.csv"))

# ------------------------------------------------------------------------------
# Figure 1: vaccine avertable burden attributable to and associated with 
# bacterial antimicrobial resistance by GBD  region, 2019

death_averted_by_region_graph <- create_burden_averted_by_region_graph(
                                    Attributable_burden_averted = Attributable_death_averted,
                                    Associated_burden_averted   = Associated_death_averted,
                                    ylim_max                    = 80000,
                                    ylabel                      = "Vaccine Avertable Deaths")

options(scipen=999)

daly_averted_by_region_graph <- create_burden_averted_by_region_graph(
                                    Attributable_burden_averted = Attributable_daly_averted,
                                    Associated_burden_averted   = Associated_daly_averted,
                                    ylim_max                    = 6000000,
                                    ylabel                      = "Vaccine Avertable DALYs")

death_averted_by_region_graph + daly_averted_by_region_graph

ggsave (filename = "Figure1_burden_averted_by_region.png",
        path = "figures",
        width = 15, 
        height = 6, 
        dpi = 600)
# ------------------------------------------------------------------------------
# Figure 2: Global vaccine avertable deaths (counts) attributable to and associated with 
# bacterial antimicrobial resistance by infectious syndrome, 2019

# create table for avertable deaths attributable to AMR by disease presentation

Attributable_death_averted_dp <- aggregate_impact_by_dp(data_input = deaths_attributable_psa)

Associated_death_averted_dp   <- aggregate_impact_by_dp(data_input = deaths_associated_psa)

death_averted_by_dp_graph <- create_burden_averted_by_dp_graph(
                              Attributable_burden_averted = Attributable_death_averted_dp,
                              Associated_burden_averted   = Associated_death_averted_dp,
                              ylim_max = 80000,
                              ylabel = "Vaccine Avertable Deaths")

# create table for avertable DALYs attributable to AMR by disease presentation

Attributable_daly_averted_dp <- aggregate_impact_by_dp(data_input = daly_attributable_psa)

Associated_daly_averted_dp   <- aggregate_impact_by_dp(data_input = daly_associated_psa)

daly_averted_by_dp_graph <- create_burden_averted_by_dp_graph(
                              Attributable_burden_averted = Attributable_daly_averted_dp,
                              Associated_burden_averted   = Associated_daly_averted_dp,
                              ylim_max = 6500000,
                              ylabel = "Vaccine Avertable DALYs")

death_averted_by_dp_graph / daly_averted_by_dp_graph

ggsave (filename = "Figure 2_burden_averted_by_dp.png",
        path = "figures",
        width = 15, 
        height = 10, 
        dpi = 600)

# ------------------------------------------------------------------------------
# Figure 3: Global vaccine avertable deaths (counts) attributable to and associated with 
# bacterial antimicrobial resistance by pathogen, 2019

# create table for avertable deaths attributable to AMR by pathogen
Attributable_death_averted_pathogen <- aggregate_impact_by_pathogen(data_input=deaths_attributable_psa)

Associated_death_averted_pathogen <- aggregate_impact_by_pathogen(data_input=deaths_associated_psa)

death_averted_by_pathogen_graph <- create_burden_averted_by_pathogen_graph(
                                    Attributable_burden_averted = Attributable_death_averted_pathogen,
                                    Associated_burden_averted   = Associated_death_averted_pathogen,
                                    ylim_max = 60000,
                                    ylabel = "Vaccine Avertable Deaths")

# create table for avertable DALYs attributable to AMR by pathogen
Attributable_daly_averted_pathogen <- aggregate_impact_by_pathogen(data_input=daly_attributable_psa)

Associated_daly_averted_pathogen <- aggregate_impact_by_pathogen(data_input=daly_associated_psa)

daly_averted_by_pathogen_graph <- create_burden_averted_by_pathogen_graph(
                                    Attributable_burden_averted = Attributable_daly_averted_pathogen,
                                    Associated_burden_averted   = Associated_daly_averted_pathogen,
                                    ylim_max = 5100000,
                                    ylabel = "Vaccine Avertable DALYs")

death_averted_by_pathogen_graph / daly_averted_by_pathogen_graph

ggsave (filename = "Figure3_burden_averted_by_pathogen.png",
        path = "figures",
        width = 15, 
        height = 10, 
        dpi = 600)
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
