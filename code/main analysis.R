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
vaccine_profile_file <- read_excel(file.path("data", "Vaccine profile assumptions.xlsx"), 
                                   sheet = "Vaccine profile _ input")

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

# Conservative scenario - vaccine impact on deaths by region

Attributable_death_averted <- aggregate_impact_by_region(input_data=deaths_attributable_psa)

Associated_death_averted   <- aggregate_impact_by_region(input_data=deaths_associated_psa)

Death_Averted              <- left_join(Associated_death_averted, Attributable_death_averted,
                                        by=c("Counts" = "Counts"))

# Conservative scenario - vaccine impact on DALYs by region

Attributable_daly_averted <- aggregate_impact_by_region(input_data = daly_attributable_psa)

Associated_daly_averted <- aggregate_impact_by_region(input_data = daly_associated_psa)

DALY_Averted <- left_join(Associated_daly_averted, Attributable_daly_averted,
                          by=c("Counts" = "Counts"))

# Optimistic scenario - vaccine impact on deaths by region

Attributable_death_averted_opt <- aggregate_impact_by_region(input_data=deaths_attributable_psa,
                                                         input_scenario="optimistic")

Associated_death_averted_opt   <- aggregate_impact_by_region(input_data=deaths_associated_psa,
                                                         input_scenario="optimistic")

Death_Averted_opt              <- left_join(Associated_death_averted_opt, Attributable_death_averted_opt,
                                            by=c("Counts" = "Counts"))

# Optimistic scenario - vaccine impact on DALYs by region
Attributable_daly_averted_opt <- aggregate_impact_by_region(input_data = daly_attributable_psa,
                                                            input_scenario = "optimistic")

Associated_daly_averted_opt   <- aggregate_impact_by_region(input_data = daly_associated_psa,
                                                            input_scenario = "optimistic")

DALY_Averted_opt              <- left_join(Associated_daly_averted_opt, Attributable_daly_averted_opt,
                                           by=c("Counts" = "Counts"))

# ------------------------------------------------------------------------------
# [table 2] Vaccine avertable deaths and DALYs associated with and attributable to 
# bacterial antimicrobial resistance globally and by WHO region, 2019

# Conservative scenario
Death_Averted <- estimate_bruden_averted_by_region(input_data = Death_Averted)

DALY_Averted <-  estimate_bruden_averted_by_region(input_data = DALY_Averted)

Table2_avertible_burden_by_region <- left_join(Death_Averted, DALY_Averted,
                                               by=c("Counts" = "Counts"))

fwrite (x    = Table2_avertible_burden_by_region, 
        file = file.path("tables", "Table2_avertible_burden_by_region.csv"))

# Optimal scenario
Death_Averted_opt <- estimate_bruden_averted_by_region(input_data = Death_Averted_opt)

DALY_Averted_opt <-  estimate_bruden_averted_by_region(input_data = DALY_Averted_opt)

Table2_avertible_burden_by_region_opt <- left_join(Death_Averted_opt, DALY_Averted_opt,
                                               by=c("Counts" = "Counts"))

fwrite (x    = Table2_avertible_burden_by_region_opt, 
        file = file.path("tables", "Table2_avertible_burden_by_region_opt.csv"))

# ------------------------------------------------------------------------------
# Figure 1: vaccine avertable burden attributable to and associated with 
# bacterial antimicrobial resistance by GBD  region, 2019
options(scipen=999)

# Conservative scenario
death_averted_by_region_graph <- create_burden_averted_by_region_graph(
                                    Attributable_burden_averted = Attributable_death_averted,
                                    Associated_burden_averted   = Associated_death_averted,
                                    ylim_max                    = 130000,
                                    ylabel                      = "Vaccine Avertable Deaths")


daly_averted_by_region_graph <- create_burden_averted_by_region_graph(
                                    Attributable_burden_averted = Attributable_daly_averted,
                                    Associated_burden_averted   = Associated_daly_averted,
                                    ylim_max                    = 11000000,
                                    ylabel                      = "Vaccine Avertable DALYs")

# Figure 1 - (1)
death_averted_by_region_graph + daly_averted_by_region_graph

ggsave (filename = "Figure1_burden_averted_by_region.png",
        path = "figures",
        width = 15, 
        height = 6, 
        dpi = 600)

# Optimistic scenario
death_averted_by_region_graph_opt <- create_burden_averted_by_region_graph(
  Attributable_burden_averted = Attributable_death_averted_opt,
  Associated_burden_averted   = Associated_death_averted_opt,
  ylim_max                    = 280000,
  ylabel                      = "Vaccine Avertable Deaths")


daly_averted_by_region_graph_opt <- create_burden_averted_by_region_graph(
  Attributable_burden_averted = Attributable_daly_averted_opt,
  Associated_burden_averted   = Associated_daly_averted_opt,
  ylim_max                    = 15000000,
  ylabel                      = "Vaccine Avertable DALYs")

# Figure 1 - (2)
death_averted_by_region_graph_opt + daly_averted_by_region_graph_opt

ggsave (filename = "Figure1_burden_averted_by_region_opt.png",
        path = "figures",
        width = 15, 
        height = 6, 
        dpi = 600)

# ------------------------------------------------------------------------------
# Figure 2: Global vaccine avertable deaths (counts) attributable to and associated with 
# bacterial antimicrobial resistance by infectious syndrome, 2019

# Conservative scenario

# create table for avertable deaths attributable to AMR by disease presentation
Attributable_death_averted_dp <- aggregate_impact_by_dp(data_input = deaths_attributable_psa)

Associated_death_averted_dp   <- aggregate_impact_by_dp(data_input = deaths_associated_psa)

death_averted_by_dp_graph <- create_burden_averted_by_dp_graph(
                              Attributable_burden_averted = Attributable_death_averted_dp,
                              Associated_burden_averted   = Associated_death_averted_dp,
                              ylim_max = 150000,
                              ylabel = "Vaccine Avertable Deaths")

# create table for avertable DALYs attributable to AMR by disease presentation

Attributable_daly_averted_dp <- aggregate_impact_by_dp(data_input = daly_attributable_psa)

Associated_daly_averted_dp   <- aggregate_impact_by_dp(data_input = daly_associated_psa)

daly_averted_by_dp_graph <- create_burden_averted_by_dp_graph(
                              Attributable_burden_averted = Attributable_daly_averted_dp,
                              Associated_burden_averted   = Associated_daly_averted_dp,
                              ylim_max = 12500000,
                              ylabel = "Vaccine Avertable DALYs")

# Figure 2 - (1)
death_averted_by_dp_graph / daly_averted_by_dp_graph

ggsave (filename = "Figure 2_burden_averted_by_dp.png",
        path = "figures",
        width = 15, 
        height = 10, 
        dpi = 600)


# Optimistic scenario

# create table for avertable deaths attributable to AMR by disease presentation
Attributable_death_averted_dp_opt <- aggregate_impact_by_dp(data_input = deaths_attributable_psa,
                                                            input_scenario="optimistic")

Associated_death_averted_dp_opt   <- aggregate_impact_by_dp(data_input = deaths_associated_psa,
                                                            input_scenario="optimistic")

death_averted_by_dp_graph_opt <- create_burden_averted_by_dp_graph(
  Attributable_burden_averted = Attributable_death_averted_dp_opt,
  Associated_burden_averted   = Associated_death_averted_dp_opt,
  ylim_max = 360000,
  ylabel = "Vaccine Avertable Deaths")

# create table for avertable DALYs attributable to AMR by disease presentation

Attributable_daly_averted_dp_opt <- aggregate_impact_by_dp(data_input = daly_attributable_psa,
                                                           input_scenario="optimistic")

Associated_daly_averted_dp_opt   <- aggregate_impact_by_dp(data_input = daly_associated_psa,
                                                           input_scenario="optimistic")

daly_averted_by_dp_graph_opt <- create_burden_averted_by_dp_graph(
  Attributable_burden_averted = Attributable_daly_averted_dp_opt,
  Associated_burden_averted   = Associated_daly_averted_dp_opt,
  ylim_max = 18000000,
  ylabel = "Vaccine Avertable DALYs")

# Figure 2 - (2)
death_averted_by_dp_graph_opt / daly_averted_by_dp_graph_opt

ggsave (filename = "Figure 2_burden_averted_by_dp_opt.png",
        path = "figures",
        width = 15, 
        height = 10, 
        dpi = 600)

# ------------------------------------------------------------------------------
# Figure 3: Global vaccine avertable deaths (counts) attributable to and associated with 
# bacterial antimicrobial resistance by pathogen, 2019

# Conservative scenario

# create table for avertable deaths attributable to AMR by pathogen
Attributable_death_averted_pathogen <- aggregate_impact_by_pathogen(data_input=deaths_attributable_psa)

Associated_death_averted_pathogen <- aggregate_impact_by_pathogen(data_input=deaths_associated_psa)

death_averted_by_pathogen_graph <- create_burden_averted_by_pathogen_graph(
                                    Attributable_burden_averted = Attributable_death_averted_pathogen,
                                    Associated_burden_averted   = Associated_death_averted_pathogen,
                                    ylim_max = 100000,
                                    ylabel = "Vaccine Avertable Deaths")

# create table for avertable DALYs attributable to AMR by pathogen
Attributable_daly_averted_pathogen <- aggregate_impact_by_pathogen(data_input=daly_attributable_psa)

Associated_daly_averted_pathogen <- aggregate_impact_by_pathogen(data_input=daly_associated_psa)

daly_averted_by_pathogen_graph <- create_burden_averted_by_pathogen_graph(
                                    Attributable_burden_averted = Attributable_daly_averted_pathogen,
                                    Associated_burden_averted   = Associated_daly_averted_pathogen,
                                    ylim_max = 9000000,
                                    ylabel = "Vaccine Avertable DALYs")

# Figure 3 - (1)
death_averted_by_pathogen_graph / daly_averted_by_pathogen_graph

ggsave (filename = "Figure3_burden_averted_by_pathogen.png",
        path = "figures",
        width = 15, 
        height = 10, 
        dpi = 600)

# Optimistic scenario

# create table for avertable deaths attributable to AMR by pathogen
Attributable_death_averted_pathogen_opt <- aggregate_impact_by_pathogen(data_input = deaths_attributable_psa,
                                                                        input_scenario = "optimistic")

Associated_death_averted_pathogen_opt <- aggregate_impact_by_pathogen(data_input = deaths_associated_psa,
                                                                      input_scenario = "optimistic")


death_averted_by_pathogen_graph_opt <- create_burden_averted_by_pathogen_graph(
  Attributable_burden_averted = Attributable_death_averted_pathogen_opt,
  Associated_burden_averted   = Associated_death_averted_pathogen_opt,
  ylim_max = 250000,
  ylabel = "Vaccine Avertable Deaths")

# create table for avertable DALYs attributable to AMR by pathogen
Attributable_daly_averted_pathogen_opt <- aggregate_impact_by_pathogen(data_input=daly_attributable_psa,
                                                                       input_scenario = "optimistic")

Associated_daly_averted_pathogen_opt <- aggregate_impact_by_pathogen(data_input=daly_associated_psa,
                                                                     input_scenario = "optimistic")

daly_averted_by_pathogen_graph_opt <- create_burden_averted_by_pathogen_graph(
  Attributable_burden_averted = Attributable_daly_averted_pathogen_opt,
  Associated_burden_averted   = Associated_daly_averted_pathogen_opt,
  ylim_max = 12500000,
  ylabel = "Vaccine Avertable DALYs")

# Figure 3 - (2)
death_averted_by_pathogen_graph_opt / daly_averted_by_pathogen_graph_opt

ggsave (filename = "Figure3_burden_averted_by_pathogen_opt.png",
        path = "figures",
        width = 15, 
        height = 10, 
        dpi = 600)

# ------------------------------------------------------------------------------
# Further analysis for pathogen with multiple vaccine options

# vaccine profile with multiple options
vaccine_profile_dt_add <- read_csv(file.path("tables", "vaccine_profile.csv"))

vaccine_profile_dt_add <- vaccine_profile_dt_add %>%
  filter(Pathogen == "Acinetobacter baumannii"|
           Pathogen == "Escherichia coli"|
           Pathogen == "Klebsiella pneumoniae"|
           Pathogen == "Mycobacterium tuberculosis"|
           Pathogen == "Streptococcus pneumoniae")

# Generate vaccine impact table for pathogen with multiple vaccine impact options
vaccine_impact_add <- bind_rows(list(
  
  create_burden_table_add(pathogen_input = "Acinetobacter baumannii",
                          vaccine_type_input = vaccine_profile_dt_add[1,]),
  
  create_burden_table_add(pathogen_input = "Acinetobacter baumannii",
                          vaccine_type_input = vaccine_profile_dt_add[2,]),
  
  create_burden_table_add(pathogen_input = "Escherichia coli",
                          vaccine_type_input = vaccine_profile_dt_add[3,]),
  
  create_burden_table_add(pathogen_input = "Escherichia coli",
                          vaccine_type_input = vaccine_profile_dt_add[4,]),
  
  create_burden_table_add(pathogen_input = "Escherichia coli",
                          vaccine_type_input = vaccine_profile_dt_add[5,]),
  
  create_burden_table_add(pathogen_input = "Klebsiella pneumoniae",
                          vaccine_type_input = vaccine_profile_dt_add[6,]),
  
  create_burden_table_add(pathogen_input = "Klebsiella pneumoniae",
                          vaccine_type_input = vaccine_profile_dt_add[7,]),
  
  create_burden_table_add(pathogen_input = "Mycobacterium tuberculosis",
                          vaccine_type_input = vaccine_profile_dt_add[8,]),
  
  create_burden_table_add(pathogen_input = "Mycobacterium tuberculosis",
                          vaccine_type_input = vaccine_profile_dt_add[9,]),
  
  create_burden_table_add(pathogen_input = "Streptococcus pneumoniae",
                          vaccine_type_input = vaccine_profile_dt_add[10,]),
  
  create_burden_table_add(pathogen_input = "Streptococcus pneumoniae",
                          vaccine_type_input = vaccine_profile_dt_add[11,])))

vaccine_impact_table_add <- cbind(vaccine_profile_dt_add, vaccine_impact_add)

fwrite (x    = vaccine_impact_table_add, 
        file = file.path("tables", "Table_vaccine_impact_for_multiple_vaccine_profiles.csv"))

# ------------------------------------------------------------------------------
# [table in appendix] vaccine avertable health burdens associated with and attributable
# to AMR by WHO region, pathogen, disease presentation, and age group
 AMR_death_data_updated <- update_death_burden(combined_dt     = combined_dt,
                                               death_burden_dt = death_burden_dt,
                                               AMR_burden_data_updated_file = file.path ("tables", "AMR deaths data_updated.csv"))

 AMR_daly_data_updated <- update_death_burden(combined_dt     = daly_combined_dt,
                                               death_burden_dt = daly_burden_dt,
                                               AMR_burden_data_updated_file = file.path ("tables", "AMR dalys data_updated.csv"))

# ------------------------------------------------------------------------------
# return to source directory
setwd (source_wd)

# end time
end_time <- Sys.time ()
print (paste0 ("end time = ", end_time))
print (Sys.time () - start_time)
# ------------------------------------------------------------------------------
