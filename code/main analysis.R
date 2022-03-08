# ------------------------------------------------------------------------------
# main analysis.R
#
# main analysis
# to estimate the vaccine avertable bacterial antimicrobial resistance burden 
# based on profiles of existing and future vaccines 
# at the global and regional levels
------------------------------------------------------------------------------
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
#              Western Pacific, Unclassified
# year: 2019
# sex: both
# ------------------------------------------------------------------------------
# estimate regional vaccine coverage for existing vaccines in 2018 and 2019
# - HIB vaccine, PCV
existing_vaccine_coverage (year = "2019",
                           hib_coverage_file = file.path("tables", "hib coverage 2019.csv"),
                           pcv_coverage_file = file.path("tables", "pcv coverage 2019.csv"))

existing_vaccine_coverage (year = "2018",
                           hib_coverage_file = file.path("tables", "hib coverage 2018.csv"),
                           pcv_coverage_file = file.path("tables", "pcv coverage 2018.csv"))

# ------------------------------------------------------------------------------
# create data table of AMR burden (deaths) classified 
# by pathogen, WHO region, disease presentation, and age groups

death_burden_dt <- read_excel(file.path("data", "IHME AMR burden.xlsx"),
                                    col_names = FALSE)             

death_burden_dt <- create_burden_table(AMR_burden = death_burden_dt,
                                       burden_file = file.path ("tables", "AMR_death_burden.csv"))

# create data table of AMR burden (DALYs) classified
# by pathogen, WHO region, disease presentation, and age groups

daly_burden_dt <- read_excel(file.path("data", "IHME AMR burden_DALYs.xlsx"),
                             col_names = FALSE)

daly_burden_dt <- create_burden_table(AMR_burden = daly_burden_dt,
                                      burden_file = file.path ("tables", "AMR_daly_burden.csv"))

# ------------------------------------------------------------------------------
# create data table of vaccine profile
vaccine_profile_file <- 
  read_excel(file.path("data", "Vaccine profile assumptions.xlsx"), 
             sheet = "Vaccine profile _ input")

vaccine_profile_dt <- 
  create_vaccine_profile_table(vaccine_profile = vaccine_profile_file,
                               vaccine_profile_file = file.path("tables", "vaccine_profile.csv"))

# ------------------------------------------------------------------------------
# create combined table: vaccine profile + disease burden (deaths)
# create separate burden file for health burdens attributable to AMR and associated with AMR  

combined_dt <- 
  create_combined_table(death_burden_dt          = death_burden_dt, 
                        vaccine_profile_dt       = vaccine_profile_dt,
                        attributable_burden_file = file.path("tables", "attributable_burden.csv"),
                        associated_burden_file   = file.path("tables", "associated_burden.csv"))

# create combined table: vaccine profile + disease burden (DALYs)
# create separate burden file for health burdens attributable to AMR and associated with AMR  
daly_combined_dt <- 
  create_combined_table(death_burden_dt          = daly_burden_dt, 
                        vaccine_profile_dt       = vaccine_profile_dt,
                        attributable_burden_file = file.path("tables", "daly_attributable_burden.csv"),
                        associated_burden_file   = file.path("tables", "daly_associated_burden.csv"))
# ------------------------------------------------------------------------------



# ------------------------------------------------------------------------------
# AMR burden by pathogen across all age groups:
# the outputs were used to find out the age group with the highest burdens
options(scipen=999)

pathogenlist_death <- unique(death_burden_dt$Pathogen)
pathogenlist_daly <- unique(daly_burden_dt$Pathogen)

lapply(pathogenlist_death, create_burden_by_pathogen_graph, 
       input_data = combined_dt, ylabel = "Number of Death Associated with AMR",
       burden_type = "deaths")

lapply(pathogenlist_daly, create_burden_by_pathogen_graph, 
       input_data = daly_combined_dt, ylabel = "Number of DALY Associated with AMR",
       burden_type = "dalys")

# ------------------------------------------------------------------------------



# ------------------------------------------------------------------------------
# estimate vaccine averted AMR deaths & uncertainty analysis

set.seed (3)  # seed for random number generator
run <- 400 # number of runs for probabilistic sensitivity analysis

# psa for deaths attributable to and associated with AMR
deaths_associated_psa <- uncertainty_analysis_baseline(
  psa       = run,
  tolerance = 0.001,
  data      = read_csv(file.path("tables", "associated_burden.csv")),
  mode      = "death")

deaths_attributable_psa <- uncertainty_analysis_baseline(
  psa       = run,
  tolerance = 0.0016,
  data      = read_csv(file.path("tables", "attributable_burden.csv")),
  mode      = "death")

# psa for DALYs attributable to and associated with AMR
daly_associated_psa <- uncertainty_analysis_baseline(
  psa       = run, 
  tolerance = 0.001,
  data      = read_csv(file.path("tables", "daly_associated_burden.csv")),
  mode      = "daly")

daly_attributable_psa <- uncertainty_analysis_baseline(
  psa       = run,
  tolerance = 0.001,
  data      = read_csv(file.path("tables", "daly_attributable_burden.csv")),
  mode      = "daly")

# ------------------------------------------------------------------------------
# deaths and DALYs associated with and attributable to AMR
# globally and by WHO region, 2019

# create table for avertable death estimates by region
# -- conservative scenario

Associated_death_averted_re <- 
  aggregate_impact_by_region(input_data=deaths_associated_psa)

Attributable_death_averted_re <- 
  aggregate_impact_by_region(input_data=deaths_attributable_psa)

# create table for avertable DALY estimates by region
# -- conservative scenario

Associated_daly_averted_re <- 
  aggregate_impact_by_region(input_data = daly_associated_psa)

Attributable_daly_averted_re <- 
  aggregate_impact_by_region(input_data = daly_attributable_psa)

# create table for avertable death estimates by region
# -- optimistic scenario

Associated_death_averted_re_opt <- 
  aggregate_impact_by_region(input_data=deaths_associated_psa,
                             input_scenario="optimistic")

Attributable_death_averted_re_opt <- 
  aggregate_impact_by_region(input_data=deaths_attributable_psa,
                             input_scenario="optimistic")

# create table for avertable DALY estimates by region
# -- optimistic scenario

Associated_daly_averted_re_opt    <- 
  aggregate_impact_by_region(input_data = daly_associated_psa,
                             input_scenario = "optimistic")

Attributable_daly_averted_re_opt  <- 
  aggregate_impact_by_region(input_data = daly_attributable_psa,
                             input_scenario = "optimistic")

# ------------------------------------------------------------------------------
# Table 2: vaccine avertable AMR health burden globally and by WHO region, 2019

Table_avertible_burden_by_region <- create_avertable_burden_table(
  Associated_death_averted       = Associated_death_averted_re,
  Attributable_death_averted     = Attributable_death_averted_re,
  Attributable_daly_averted      = Attributable_daly_averted_re,
  Associated_daly_averted        = Associated_daly_averted_re,
  Attributable_death_averted_opt = Attributable_death_averted_re_opt, 
  Associated_death_averted_opt   = Associated_death_averted_re_opt,
  Attributable_daly_averted_opt  = Attributable_daly_averted_re_opt, 
  Associated_daly_averted_opt    = Associated_daly_averted_re_opt)

fwrite (x    = Table_avertible_burden_by_region,
        file = file.path("tables", "Table2_avertible_burden_by_region.csv"))

# ------------------------------------------------------------------------------
# Figure 1: Vaccine impact by WHO region, 2019

# create graph with avertable AMR burden by region -- conservative scenario
death_averted_by_region_graph <- create_burden_averted_by_region_graph(
  Attributable_burden_averted = Attributable_death_averted_re,
  Associated_burden_averted   = Associated_death_averted_re,
  ylim_max                    = 130000,
  ylabel                      = "Vaccine Avertable Deaths",
  title_name                  = "Conservative Scenario")

daly_averted_by_region_graph <- create_burden_averted_by_region_graph(
  Attributable_burden_averted = Attributable_daly_averted_re,
  Associated_burden_averted   = Associated_daly_averted_re,
  ylim_max                    = 11000000,
  ylabel                      = "Vaccine Avertable DALYs",
  title_name                  = " ")

# create graph with avertable AMR burden by region -- optimistic scenario
death_averted_by_region_graph_opt <- create_burden_averted_by_region_graph(
  Attributable_burden_averted = Attributable_death_averted_re_opt,
  Associated_burden_averted   = Associated_death_averted_re_opt,
  ylim_max                    = 280000,
  ylabel                      = "Vaccine Avertable Deaths",
  title_name                  = "Optimistic Scenario")


daly_averted_by_region_graph_opt <- create_burden_averted_by_region_graph(
  Attributable_burden_averted = Attributable_daly_averted_re_opt,
  Associated_burden_averted   = Associated_daly_averted_re_opt,
  ylim_max                    = 15000000,
  ylabel                      = "Vaccine Avertable DALYs",
  title_name                  = " ")

# create [Figure 1]
(death_averted_by_region_graph + daly_averted_by_region_graph) /
(death_averted_by_region_graph_opt + daly_averted_by_region_graph_opt)

# save the image file
ggsave (filename = "Figure1_burden_averted_by_region.png",
        path = "figures",
        width = 15, 
        height = 12, 
        dpi = 600)

# ------------------------------------------------------------------------------
# global vaccine avertable deaths and DALYs attributable to and associated with 
# bacterial antimicrobial resistance by infectious syndrome, 2019

# disease presentation -- Neisseria gonorrhoeae data is only available for DALYs
DiseasePresentation_death <- unique(death_burden_dt$Disease_presentation)
DiseasePresentation_daly  <- unique(daly_burden_dt$Disease_presentation)

# create table for avertable death estimates by disease presentation
# -- conservative scenario
Associated_death_averted_dp   <- 
  aggregate_impact_by_dp(data_input          = deaths_associated_psa,
                         DiseasePresentation = DiseasePresentation_death)

Attributable_death_averted_dp <- 
  aggregate_impact_by_dp(data_input          = deaths_attributable_psa,
                         DiseasePresentation = DiseasePresentation_death)

# create table for avertable DALY estimates by disease presentation
# -- conservative scenario

Associated_daly_averted_dp   <- 
  aggregate_impact_by_dp(data_input          = daly_associated_psa,
                         DiseasePresentation = DiseasePresentation_daly)

Attributable_daly_averted_dp <- 
  aggregate_impact_by_dp(data_input          = daly_attributable_psa,
                         DiseasePresentation = DiseasePresentation_daly)

# create table for avertable death estimates by disease presentation
# -- optimistic scenario

Associated_death_averted_dp_opt   <- aggregate_impact_by_dp(
  data_input          = deaths_associated_psa,
  input_scenario      = "optimistic",
  DiseasePresentation = DiseasePresentation_death)

Attributable_death_averted_dp_opt <- aggregate_impact_by_dp(
  data_input          = deaths_attributable_psa,
  input_scenario      = "optimistic",
  DiseasePresentation = DiseasePresentation_death)

# create table for avertable DALY estimates by disease presentation
# -- optimistic scenario

Associated_daly_averted_dp_opt   <- aggregate_impact_by_dp(
  data_input          = daly_associated_psa,
  input_scenario      = "optimistic",
  DiseasePresentation = DiseasePresentation_daly)

Attributable_daly_averted_dp_opt <- aggregate_impact_by_dp(
  data_input          = daly_attributable_psa,
  input_scenario      = "optimistic",
  DiseasePresentation = DiseasePresentation_daly)

# ------------------------------------------------------------------------------
# Table: vaccine avertable AMR health burden by infectious syndrome, 2019

Table_avertible_burden_by_dp <- create_avertable_burden_table(
  Associated_death_averted       = Associated_death_averted_dp,
  Attributable_death_averted     = Attributable_death_averted_dp,
  Attributable_daly_averted      = Attributable_daly_averted_dp,
  Associated_daly_averted        = Associated_daly_averted_dp,
  Attributable_death_averted_opt = Attributable_death_averted_dp_opt, 
  Associated_death_averted_opt   = Associated_death_averted_dp_opt,
  Attributable_daly_averted_opt  = Attributable_daly_averted_dp_opt, 
  Associated_daly_averted_opt    = Associated_daly_averted_dp_opt)

Table_avertible_burden_by_dp <- 
  Table_avertible_burden_by_dp %>% arrange(Counts)

fwrite (x    = Table_avertible_burden_by_dp,
        file = file.path("tables", "Table_avertible_burden_by_infectious_syndrome.csv"))

# ------------------------------------------------------------------------------
# Figure 2: Vaccine impact by infectious syndrome, 2019

# create graph with avertable AMR burden by infectious syndrome 
# -- conservative scenario
death_averted_by_dp_graph <- create_burden_averted_by_dp_graph(
  Attributable_burden_averted = Attributable_death_averted_dp,
  Associated_burden_averted   = Associated_death_averted_dp,
  ylim_max = 150000,
  ylabel = "Vaccine Avertable Deaths",
  title_name = "Conservative Scenario")

daly_averted_by_dp_graph <- create_burden_averted_by_dp_graph(
  Attributable_burden_averted = Attributable_daly_averted_dp,
  Associated_burden_averted   = Associated_daly_averted_dp,
  ylim_max = 12500000,
  ylabel = "Vaccine Avertable DALYs",
  title_name = " ")

# create graph with avertable AMR burden by infectious syndrome 
# -- optimistic scenario
death_averted_by_dp_graph_opt <- create_burden_averted_by_dp_graph(
  Attributable_burden_averted = Attributable_death_averted_dp_opt,
  Associated_burden_averted   = Associated_death_averted_dp_opt,
  ylim_max = 360000,
  ylabel = "Vaccine Avertable Deaths",
  title_name = "Optimistic Scenario")

daly_averted_by_dp_graph_opt <- create_burden_averted_by_dp_graph(
  Attributable_burden_averted = Attributable_daly_averted_dp_opt,
  Associated_burden_averted   = Associated_daly_averted_dp_opt,
  ylim_max = 18000000,
  ylabel = "Vaccine Avertable DALYs",
  title_name = " ")

# create [Figure 2]
death_averted_by_dp_graph / daly_averted_by_dp_graph/
death_averted_by_dp_graph_opt / daly_averted_by_dp_graph_opt

# save the image file
ggsave (filename = "Figure 2_burden_averted_by_dp.png",
        path = "figures",
        width = 15, 
        height = 17, 
        dpi = 600)

# ------------------------------------------------------------------------------
# global vaccine avertable deaths and DALYs attributable to and associated with 
# bacterial antimicrobial resistance by pathogen, 2019

# create table for avertable death estimates by pathogen
# -- conservative scenario
Associated_death_averted_pathogen <- 
  aggregate_impact_by_pathogen(data_input=deaths_associated_psa,
                               pathogenlist=pathogenlist_death)

Attributable_death_averted_pathogen <- 
  aggregate_impact_by_pathogen(data_input=deaths_attributable_psa,
                               pathogenlist=pathogenlist_death)

# create table for avertable DALY estimates by pathogen
# -- conservative scenario
Attributable_daly_averted_pathogen <- 
  aggregate_impact_by_pathogen(data_input=daly_attributable_psa,
                               pathogenlist=pathogenlist_daly)

Associated_daly_averted_pathogen <- 
  aggregate_impact_by_pathogen(data_input=daly_associated_psa,
                               pathogenlist=pathogenlist_daly)

# create table for avertable death estimates by pathogen
# -- optimistic scenario
Attributable_death_averted_pathogen_opt <- 
  aggregate_impact_by_pathogen(data_input     = deaths_attributable_psa,
                               input_scenario = "optimistic",
                               pathogenlist   = pathogenlist_death)

Associated_death_averted_pathogen_opt <- 
  aggregate_impact_by_pathogen(data_input     = deaths_associated_psa,
                               input_scenario = "optimistic",
                               pathogenlist   = pathogenlist_death)

# create table for avertable DALY estimates by pathogen
# -- optimistic scenario
Attributable_daly_averted_pathogen_opt <- 
  aggregate_impact_by_pathogen(data_input     = daly_attributable_psa,
                               input_scenario = "optimistic",
                               pathogenlist   = pathogenlist_daly)

Associated_daly_averted_pathogen_opt  <- 
  aggregate_impact_by_pathogen(data_input     = daly_associated_psa,
                               input_scenario = "optimistic",
                               pathogenlist   = pathogenlist_daly)

# -------------------------------------------------------------------------
# Table: vaccine avertable AMR health burden by pathogen, 2019

Table_avertible_burden_by_pathogen <- create_avertable_burden_table(
  Associated_death_averted       = Associated_death_averted_pathogen,
  Attributable_death_averted     = Attributable_death_averted_pathogen,
  Attributable_daly_averted      = Attributable_daly_averted_pathogen,
  Associated_daly_averted        = Associated_daly_averted_pathogen,
  Attributable_death_averted_opt = Attributable_death_averted_pathogen_opt, 
  Associated_death_averted_opt   = Associated_death_averted_pathogen_opt,
  Attributable_daly_averted_opt  = Attributable_daly_averted_pathogen_opt, 
  Associated_daly_averted_opt    = Associated_daly_averted_pathogen_opt)

Table_avertible_burden_by_pathogen <- 
  Table_avertible_burden_by_pathogen %>% arrange(Counts)

fwrite (x    = Table_avertible_burden_by_pathogen,
        file = file.path("tables", "Table_avertible_burden_by_pathogen.csv"))

# -------------------------------------------------------------------------
# Figure 3: vaccine impact by pathogen, 2019

# create graph with avertable AMR burden by pathogen -- conservative scenario
death_averted_by_pathogen_graph <- create_burden_averted_by_pathogen_graph(
  Attributable_burden_averted = Attributable_death_averted_pathogen,
  Associated_burden_averted   = Associated_death_averted_pathogen,
  ylim_max = 100000,
  ylabel = "Vaccine Avertable Deaths",
  title_name = "Conservative Scenario")

daly_averted_by_pathogen_graph <- create_burden_averted_by_pathogen_graph(
  Attributable_burden_averted = Attributable_daly_averted_pathogen,
  Associated_burden_averted   = Associated_daly_averted_pathogen,
  ylim_max = 9000000,
  ylabel = "Vaccine Avertable DALYs",
  title_name = " ")

# create graph with avertable AMR burden by pathogen -- optimistic scenario
death_averted_by_pathogen_graph_opt <- create_burden_averted_by_pathogen_graph(
  Attributable_burden_averted = Attributable_death_averted_pathogen_opt,
  Associated_burden_averted   = Associated_death_averted_pathogen_opt,
  ylim_max = 250000,
  ylabel = "Vaccine Avertable Deaths",
  title_name = "Optimistic Scenario")

daly_averted_by_pathogen_graph_opt <- create_burden_averted_by_pathogen_graph(
  Attributable_burden_averted = Attributable_daly_averted_pathogen_opt,
  Associated_burden_averted   = Associated_daly_averted_pathogen_opt,
  ylim_max = 12500000,
  ylabel = "Vaccine Avertable DALYs",
  title_name = " ")

# create [Figure 3]
death_averted_by_pathogen_graph/daly_averted_by_pathogen_graph/
death_averted_by_pathogen_graph_opt/daly_averted_by_pathogen_graph_opt

# save the image file
ggsave (filename = "Figure3_burden_averted_by_pathogen.png",
        path = "figures",
        width = 15, 
        height = 18, 
        dpi = 600)

# ------------------------------------------------------------------------------
# further analysis for pathogen with multiple vaccine options

# vaccine profile with multiple options
vaccine_profile_dt_add <- read_csv(file.path("tables", "vaccine_profile.csv"))

vaccine_profile_dt_add <- vaccine_profile_dt_add %>%
  filter(Pathogen == "Acinetobacter baumannii"|
           Pathogen == "Escherichia coli"|
           Pathogen == "Klebsiella pneumoniae"|
           Pathogen == "Mycobacterium tuberculosis"|
           Pathogen == "Streptococcus pneumoniae")

# generate vaccine impact table for pathogen with multiple vaccine impact options
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

# vaccine avertable AMR burden estimates using the burden estimate median value

# create vaccine avertable burden file -- conservative scenario
AMR_death_data_updated <- update_death_burden(
  combined_dt                  = combined_dt,
  burden_dt                    = death_burden_dt,
  AMR_burden_data_updated_file = file.path ("tables", "AMR deaths data_conservative.csv"),
  input_scenario               = "conservative")

AMR_daly_data_updated <- update_death_burden(
  combined_dt                  = daly_combined_dt,
  burden_dt                    = daly_burden_dt,
  AMR_burden_data_updated_file = file.path ("tables", "AMR dalys data_conservative.csv"),
  input_scenario               = "conservative")

# create vaccine avertable burden file -- optimistic scenario
 AMR_death_data_updated_opt <- update_death_burden(
   combined_dt                  = combined_dt,
   burden_dt                    = death_burden_dt,
   AMR_burden_data_updated_file = file.path ("tables", "AMR deaths data_optimistic.csv"),
   input_scenario               = "optimistic")

 AMR_daly_data_updated_opt <- update_death_burden(
   combined_dt                  = daly_combined_dt,
   burden_dt                    = daly_burden_dt,
   AMR_burden_data_updated_file = file.path ("tables", "AMR dalys data_optimistic.csv"),
   input_scenario               = "optimistic")
 
 
# ------------------------------------------------------------------------------
# return to source directory
setwd (source_wd)

# end time
end_time <- Sys.time ()
print (paste0 ("end time = ", end_time))
print (Sys.time () - start_time)
# ------------------------------------------------------------------------------
