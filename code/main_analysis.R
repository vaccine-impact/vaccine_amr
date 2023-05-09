# ----------------------------------------------------------------------------
# main analysis.R
#
# main analysis
# to estimate the vaccine avertable bacterial antimicrobial resistance burden
# based on profiles of existing and future vaccines
# at the global and regional levels
# -----------------------------------------------------------------------------
# load libraries
library (readr)
library (readxl)
library (dplyr)
library (ggplot2)
library (reshape2)
library (tidyverse)
library (scales)
library (data.table)
library (ggplot2)
library (rriskDistributions)
library (msm)
library (formattable)
library (patchwork)
library (Cairo)
library (webr)

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
                           hib_coverage_file = file.path("tables", "hib_coverage_2019.csv"),
                           pcv_coverage_file = file.path("tables", "pcv_coverage_2019.csv"))

existing_vaccine_coverage (year = "2018",
                           hib_coverage_file = file.path("tables", "hib_coverage_2018.csv"),
                           pcv_coverage_file = file.path("tables", "pcv_coverage_2018.csv"))

# ------------------------------------------------------------------------------
# create data table of AMR burden (deaths) classified 
# by pathogen, WHO region, disease presentation, and age groups
options(scipen=999)

death_burden_dt <- read_excel(file.path("data", "IHME_AMR_burden.xlsx"),
                              col_names = FALSE)             

death_burden_dt <- create_burden_table(AMR_burden  = death_burden_dt,
                                       burden_file = file.path ("tables", "AMR_death_burden.csv"))

# create data table of AMR burden (DALYs) classified
# by pathogen, WHO region, disease presentation, and age groups

daly_burden_dt <- read_excel(file.path("data", "IHME_AMR_burden_DALYs.xlsx"),
                             col_names = FALSE)

daly_burden_dt <- create_burden_table(AMR_burden = daly_burden_dt,
                                      burden_file = file.path ("tables", "AMR_daly_burden.csv"))

# ------------------------------------------------------------------------------
# create data table of vaccine profile
vaccine_profile_file <-
  read_excel(file.path("data", "Vaccine_profile_assumptions.xlsx"),
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
                        susceptible_burden_file  = file.path("tables", "susceptible_burden.csv"),
                        attributable_burden_file = file.path("tables", "attributable_burden.csv"),
                        associated_burden_file   = file.path("tables", "associated_burden.csv"))

# create combined table: vaccine profile + disease burden (DALYs)
# create separate burden file for health burdens attributable to AMR and associated with AMR  
daly_combined_dt <-
  create_combined_table(death_burden_dt          = daly_burden_dt, 
                        vaccine_profile_dt       = vaccine_profile_dt,
                        susceptible_burden_file  = file.path("tables", "daly_susceptible_burden.csv"),
                        attributable_burden_file = file.path("tables", "daly_attributable_burden.csv"),
                        associated_burden_file   = file.path("tables", "daly_associated_burden.csv"))
# ------------------------------------------------------------------------------
# estimate pre-vaccine burden for existing vaccines -- HIB vaccine & PCV
estimate_prevaccination_burden(
  burden_input = read_csv(file.path("tables", "associated_burden.csv")),
  burden_file  = file.path("tables", "associated_burden.csv"))

estimate_prevaccination_burden(
  burden_input = read_csv(file.path("tables", "attributable_burden.csv")),
  burden_file  = file.path("tables", "attributable_burden.csv"))

estimate_prevaccination_burden(
  burden_input = read_csv(file.path("tables", "susceptible_burden.csv")),
  burden_file  = file.path("tables", "susceptible_burden.csv"))

estimate_prevaccination_burden(
  burden_input = read_csv(file.path("tables", "daly_associated_burden.csv")),
  burden_file  = file.path("tables", "daly_associated_burden.csv"))

estimate_prevaccination_burden(
  burden_input = read_csv(file.path("tables", "daly_attributable_burden.csv")),
  burden_file  = file.path("tables", "daly_attributable_burden.csv"))

estimate_prevaccination_burden(
  burden_input = read_csv(file.path("tables", "daly_susceptible_burden.csv")),
  burden_file  = file.path("tables", "daly_susceptible_burden.csv"))
# ------------------------------------------------------------------------------



# ------------------------------------------------------------------------------
# AMR burden by pathogen across all age groups:
# the outputs were used to find out the age group with the highest burden

# pathogenlist_death <- unique(death_burden_dt$Pathogen)
# pathogenlist_daly  <- unique(daly_burden_dt$Pathogen)

# lapply(pathogenlist_death, create_burden_by_pathogen_graph, 
#        input_data = combined_dt, ylabel = "Number of Death Associated with AMR",
#        burden_type = "deaths")

# lapply(pathogenlist_daly, create_burden_by_pathogen_graph, 
#        input_data = daly_combined_dt, ylabel = "Number of DALY Associated with AMR",
#        burden_type = "dalys")
# ------------------------------------------------------------------------------



# ------------------------------------------------------------------------------
# estimate vaccine averted AMR burden & uncertainty analysis baseline
memory.limit(size = 20000)

set.seed (3)  # seed for random number generator
run <- 400 # number of runs for probabilistic sensitivity analysis

# baseline for uncertainty analysis -- deaths associated with AMR
deaths_associated_psa <- uncertainty_analysis_baseline(
  psa       = run,
  tolerance = 0.001,
  data      = read_csv(file.path("tables", "associated_burden.csv")),
  mode      = "death")

fwrite (x    = deaths_associated_psa,
        file = file.path("tables", "deaths_associated_psa.csv"))

deaths_associated_psa$va_base <- 
  estimate_vaccine_impact(data     = read_csv(file.path("tables", "deaths_associated_psa.csv")),
                          scenario = "conservative")[, va_health_burden]

deaths_associated_psa$va_high <- 
  estimate_vaccine_impact(data     = read_csv(file.path("tables", "deaths_associated_psa.csv")),
                          scenario = "optimistic")[, va_health_burden]

deaths_associated_psa$va_incre <- 
  deaths_associated_psa$va_high - deaths_associated_psa$va_base

fwrite (x    = deaths_associated_psa,
        file = file.path("tables", "deaths_associated_psa.csv"))

# baseline for uncertainty analysis -- deaths attributable to AMR

deaths_attributable_psa <- uncertainty_analysis_baseline(
  psa       = run,
  tolerance = 0.0016,
  data      = read_csv(file.path("tables", "attributable_burden.csv")),
  mode      = "death")

fwrite (x    = deaths_attributable_psa,
        file = file.path("tables", "deaths_attributable_psa.csv"))

deaths_attributable_psa$va_base <- 
  estimate_vaccine_impact(data     = read_csv(file.path("tables", "deaths_attributable_psa.csv")),
                          scenario = "conservative")[, va_health_burden]

deaths_attributable_psa$va_high <- 
  estimate_vaccine_impact(data     = read_csv(file.path("tables", "deaths_attributable_psa.csv")),
                          scenario = "optimistic")[, va_health_burden]

deaths_attributable_psa$va_incre <- 
  deaths_attributable_psa$va_high - deaths_attributable_psa$va_base

fwrite (x    = deaths_attributable_psa,
        file = file.path("tables", "deaths_attributable_psa.csv"))

# baseline for uncertainty analysis -- DALYs associated with AMR

daly_associated_psa <- uncertainty_analysis_baseline(
  psa       = run, 
  tolerance = 0.001,
  data      = read_csv(file.path("tables", "daly_associated_burden.csv")),
  mode      = "daly")

fwrite (x    = daly_associated_psa,
        file = file.path("tables", "daly_associated_psa.csv"))

daly_associated_psa$va_base <- 
  estimate_vaccine_impact(data     = read_csv(file.path("tables", "daly_associated_psa.csv")),
                          scenario = "conservative")[, va_health_burden]

daly_associated_psa$va_high <- 
  estimate_vaccine_impact(data     = read_csv(file.path("tables", "daly_associated_psa.csv")),
                          scenario = "optimistic")[, va_health_burden]

daly_associated_psa$va_incre <- 
  daly_associated_psa$va_high - daly_associated_psa$va_base

fwrite (x    = daly_associated_psa,
        file = file.path("tables", "daly_associated_psa.csv"))

# baseline for uncertainty analysis -- DALYs attributable to AMR

daly_attributable_psa <- uncertainty_analysis_baseline(
  psa       = run,
  tolerance = 0.001,
  data      = read_csv(file.path("tables", "daly_attributable_burden.csv")),
  mode      = "daly")

fwrite (x    = daly_attributable_psa,
        file = file.path("tables", "daly_attributable_psa.csv"))

daly_attributable_psa$va_base <- 
  estimate_vaccine_impact(data     = read_csv(file.path("tables", "daly_attributable_psa.csv")),
                          scenario = "conservative")[, va_health_burden]

daly_attributable_psa$va_high <- 
  estimate_vaccine_impact(data     = read_csv(file.path("tables", "daly_attributable_psa.csv")),
                          scenario = "optimistic")[, va_health_burden]

daly_attributable_psa$va_incre <- 
  daly_attributable_psa$va_high - daly_attributable_psa$va_base

fwrite (x    = daly_attributable_psa,
        file = file.path("tables", "daly_attributable_psa.csv"))

# baseline for uncertainty analysis -- deaths associated with susceptible
 deaths_susceptible_psa <- uncertainty_analysis_baseline(
  psa       = run,
  tolerance = 0.001,
  data      = read_csv(file.path("tables", "susceptible_burden.csv")),
  mode      = "death")

 fwrite (x    = deaths_susceptible_psa,
        file = file.path("tables", "deaths_susceptible_psa.csv"))

 deaths_susceptible_psa$va_base <- 
  estimate_vaccine_impact(data     = read_csv(file.path("tables", "deaths_susceptible_psa.csv")),
                          scenario = "conservative")[, va_health_burden]

 deaths_susceptible_psa$va_high <- 
  estimate_vaccine_impact(data     = read_csv(file.path("tables", "deaths_susceptible_psa.csv")),
                          scenario = "optimistic")[, va_health_burden]

 deaths_susceptible_psa$va_incre <- 
  deaths_susceptible_psa$va_high - deaths_susceptible_psa$va_base

 fwrite (x    = deaths_susceptible_psa,
        file = file.path("tables", "deaths_susceptible_psa.csv"))

# baseline for uncertainty analysis -- DALYs associated with susceptible
 daly_susceptible_psa <- uncertainty_analysis_baseline(
   psa       = run, 
   tolerance = 0.001,
   data      = read_csv(file.path("tables", "daly_susceptible_burden.csv")),
   mode      = "daly")

 fwrite (x    = daly_susceptible_psa,
         file = file.path("tables", "daly_susceptible_psa.csv"))

 daly_susceptible_psa$va_base <- 
   estimate_vaccine_impact(data     = read_csv(file.path("tables", "daly_susceptible_psa.csv")),
                           scenario = "conservative")[, va_health_burden]

 daly_susceptible_psa$va_high <- 
   estimate_vaccine_impact(data     = read_csv(file.path("tables", "daly_susceptible_psa.csv")),
                           scenario = "optimistic")[, va_health_burden]

 daly_susceptible_psa$va_incre <- 
   daly_susceptible_psa$va_high - daly_susceptible_psa$va_base

 fwrite (x    = daly_susceptible_psa,
         file = file.path("tables", "daly_susceptible_psa.csv"))
# ------------------------------------------------------------------------------



# ------------------------------------------------------------------------------
# deaths and DALYs associated with and attributable to AMR
# globally and by WHO region, 2019

# create table for avertable death estimates by region
# -- Baseline Scenario

Associated_death_averted_re <- aggregate_impact_by_region(
  input_data = read_csv(file.path("tables", "deaths_associated_psa.csv")),
  mode       = "baseline")

Attributable_death_averted_re <- aggregate_impact_by_region(
  input_data = read_csv(file.path("tables", "deaths_attributable_psa.csv")),
  mode       = "baseline")

Susceptible_death_averted_re <- aggregate_impact_by_region(
  input_data = read_csv(file.path("tables", "deaths_susceptible_psa.csv")),
  mode       = "baseline")

# create table for avertable DALY estimates by region
# -- Baseline Scenario

Associated_daly_averted_re <- aggregate_impact_by_region(
  input_data = read_csv(file.path("tables", "daly_associated_psa.csv")),
  mode       = "baseline")

Attributable_daly_averted_re <- aggregate_impact_by_region(
  input_data = read_csv(file.path("tables", "daly_attributable_psa.csv")),
  mode       = "baseline")

Susceptible_daly_averted_re <- aggregate_impact_by_region(
  input_data = read_csv(file.path("tables", "daly_susceptible_psa.csv")),
  mode       = "baseline")

# create table for avertable death estimates by region
# -- High-potential Scenario

Associated_death_averted_re_opt <- aggregate_impact_by_region(
  input_data = read_csv(file.path("tables", "deaths_associated_psa.csv")),
  mode       = "high_potential")

Attributable_death_averted_re_opt <- aggregate_impact_by_region(
  input_data = read_csv(file.path("tables", "deaths_attributable_psa.csv")),
  mode       = "high_potential")

Susceptible_death_averted_re_opt <- aggregate_impact_by_region(
  input_data = read_csv(file.path("tables", "deaths_susceptible_psa.csv")),
  mode       = "high_potential")

# create table for avertable DALY estimates by region
# -- High-potential Scenario

Associated_daly_averted_re_opt <- aggregate_impact_by_region(
  input_data = read_csv(file.path("tables", "daly_associated_psa.csv")),
  mode       = "high_potential")

Attributable_daly_averted_re_opt <- aggregate_impact_by_region(
  input_data = read_csv(file.path("tables", "daly_attributable_psa.csv")),
  mode       = "high_potential")

Susceptible_daly_averted_re_opt <- aggregate_impact_by_region(
  input_data = read_csv(file.path("tables", "daly_susceptible_psa.csv")),
  mode       = "high_potential")

# create table for incremental avertable deaths estimated by region
# -- (High-potential Scenario - Baseline Scenario)
Associated_death_averted_re_inc <- aggregate_impact_by_region(
  input_data = read_csv(file.path("tables", "deaths_associated_psa.csv")),
  mode       = "incremental")

Attributable_death_averted_re_inc <- aggregate_impact_by_region(
  input_data = read_csv(file.path("tables", "deaths_attributable_psa.csv")),
  mode       = "incremental")

Susceptible_death_averted_re_inc <- aggregate_impact_by_region(
  input_data = read_csv(file.path("tables", "deaths_susceptible_psa.csv")),
  mode       = "incremental")

# create table for incremental avertable DALY estimated by region
# -- (High-potential Scenario - Baseline Scenario)

Associated_daly_averted_re_inc <- aggregate_impact_by_region(
  input_data = read_csv(file.path("tables", "daly_associated_psa.csv")),
  mode       = "incremental")

Attributable_daly_averted_re_inc <- aggregate_impact_by_region(
  input_data = read_csv(file.path("tables", "daly_attributable_psa.csv")),
  mode       = "incremental")

Susceptible_daly_averted_re_inc <- aggregate_impact_by_region(
  input_data = read_csv(file.path("tables", "daly_susceptible_psa.csv")),
  mode       = "incremental")
# ------------------------------------------------------------------------------
# Table: vaccine avertable AMR health burden globally and by WHO region, 2019

Table_avertable_burden_by_region <- create_avertable_burden_table(
  Associated_death_averted       = Associated_death_averted_re,
  Associated_death_averted_opt   = Associated_death_averted_re_opt,
  Associated_death_averted_inc   = Associated_death_averted_re_inc,
  Associated_daly_averted        = Associated_daly_averted_re,
  Associated_daly_averted_opt    = Associated_daly_averted_re_opt,
  Associated_daly_averted_inc    = Associated_daly_averted_re_inc,
  Attributable_death_averted     = Attributable_death_averted_re,
  Attributable_death_averted_opt = Attributable_death_averted_re_opt, 
  Attributable_death_averted_inc = Attributable_death_averted_re_inc, 
  Attributable_daly_averted      = Attributable_daly_averted_re,
  Attributable_daly_averted_opt  = Attributable_daly_averted_re_opt,
  Attributable_daly_averted_inc  = Attributable_daly_averted_re_inc,
  Susceptible_death_averted      = Susceptible_death_averted_re,
  Susceptible_death_averted_opt  = Susceptible_death_averted_re_opt,
  Susceptible_death_averted_inc  = Susceptible_death_averted_re_inc,
  Susceptible_daly_averted       = Susceptible_daly_averted_re,
  Susceptible_daly_averted_opt   = Susceptible_daly_averted_re_opt,
  Susceptible_daly_averted_inc   = Susceptible_daly_averted_re_inc)

fwrite (x    = Table_avertable_burden_by_region,
        file = file.path("tables", "Table_avertable_burden_by_region.csv"))

# ------------------------------------------------------------------------------
# Figure: vaccine impact by WHO region, 2019

# create graph with avertable AMR burden by region -- Baseline Scenario
death_averted_by_region_graph <- create_burden_averted_by_region_graph(
  Attributable_burden_averted = Attributable_death_averted_re,
  Associated_burden_averted   = Associated_death_averted_re,
  ylabel                      = "Vaccine Avertable Deaths")

daly_averted_by_region_graph <- create_burden_averted_by_region_graph(
  Attributable_burden_averted = Attributable_daly_averted_re,
  Associated_burden_averted   = Associated_daly_averted_re,
  ylabel                      = "Vaccine Avertable DALYs")

# create figure
death_averted_by_region_graph / daly_averted_by_region_graph

# save the image file
ggsave (filename = "Figure_avertable_burden_by_region.png",
        path = "figures",
        width = 6, 
        height = 8, 
        dpi = 600)

ggsave (filename = "Figure_avertable_burden_by_region.eps",
        path = "figures",
        device = "eps",
        width = 7, 
        height = 7)

# ------------------------------------------------------------------------------
# global vaccine avertable deaths and DALYs attributable to and associated with 
# bacterial antimicrobial resistance by infectious syndrome, 2019

# disease presentation -- Neisseria gonorrhoeae data is only available for DALYs
DiseasePresentation_death <- unique(death_burden_dt$Disease_presentation)
DiseasePresentation_daly  <- unique(daly_burden_dt$Disease_presentation)

# create table for avertable death estimates by disease presentation
# -- Baseline Scenario
Associated_death_averted_dp   <- aggregate_impact_by_dp(
  input_data          = read_csv(file.path("tables", "deaths_associated_psa.csv")),
  DiseasePresentation = DiseasePresentation_death,
  mode                = "baseline")

Attributable_death_averted_dp <- aggregate_impact_by_dp(
  input_data          = read_csv(file.path("tables", "deaths_attributable_psa.csv")),
  DiseasePresentation = DiseasePresentation_death,
  mode                = "baseline")

Susceptible_death_averted_dp <- aggregate_impact_by_dp(
  input_data          = read_csv(file.path("tables", "deaths_susceptible_psa.csv")),
  DiseasePresentation = DiseasePresentation_death,
  mode                = "baseline")

# create table for avertable DALY estimates by disease presentation
# -- Baseline Scenario
Associated_daly_averted_dp   <- aggregate_impact_by_dp(
  input_data          = read_csv(file.path("tables", "daly_associated_psa.csv")),
  DiseasePresentation = DiseasePresentation_daly,
  mode                = "baseline")

Attributable_daly_averted_dp <- aggregate_impact_by_dp(
  input_data          = read_csv(file.path("tables", "daly_attributable_psa.csv")),
  DiseasePresentation = DiseasePresentation_daly,
  mode                = "baseline")

Susceptible_daly_averted_dp <- aggregate_impact_by_dp(
  input_data          = read_csv(file.path("tables", "daly_susceptible_psa.csv")),
  DiseasePresentation = DiseasePresentation_daly,
  mode                = "baseline")

# create table for avertable death estimates by disease presentation
# -- High-potential Scenario
Associated_death_averted_dp_opt   <- aggregate_impact_by_dp(
  input_data          = read_csv(file.path("tables", "deaths_associated_psa.csv")),
  DiseasePresentation = DiseasePresentation_death,
  mode                = "high_potential")

Attributable_death_averted_dp_opt <- aggregate_impact_by_dp(
  input_data          = read_csv(file.path("tables", "deaths_attributable_psa.csv")),
  DiseasePresentation = DiseasePresentation_death,
  mode                = "high_potential")

Susceptible_death_averted_dp_opt <- aggregate_impact_by_dp(
  input_data          = read_csv(file.path("tables", "deaths_susceptible_psa.csv")),
  DiseasePresentation = DiseasePresentation_death,
  mode                = "high_potential")

# create table for avertable DALY estimates by disease presentation
# -- High-potential Scenario
Associated_daly_averted_dp_opt   <- aggregate_impact_by_dp(
  input_data          = read_csv(file.path("tables", "daly_associated_psa.csv")),
  DiseasePresentation = DiseasePresentation_daly,
  mode                = "high_potential")

Attributable_daly_averted_dp_opt <- aggregate_impact_by_dp(
  input_data          = read_csv(file.path("tables", "daly_attributable_psa.csv")),
  DiseasePresentation = DiseasePresentation_daly,
  mode                = "high_potential")

Susceptible_daly_averted_dp_opt <- aggregate_impact_by_dp(
  input_data          = read_csv(file.path("tables", "daly_susceptible_psa.csv")),
  DiseasePresentation = DiseasePresentation_daly,
  mode                = "high_potential")

# create table for incremental avertable death estimates by disease presentation
# -- (High-potential Scenario - Baseline Scenario)
Associated_death_averted_dp_inc   <- aggregate_impact_by_dp(
  input_data          = read_csv(file.path("tables", "deaths_associated_psa.csv")),
  DiseasePresentation = DiseasePresentation_death,
  mode                = "incremental")

Attributable_death_averted_dp_inc <- aggregate_impact_by_dp(
  input_data          = read_csv(file.path("tables", "daly_attributable_psa.csv")),
  DiseasePresentation = DiseasePresentation_death,
  mode                = "incremental")

Susceptible_death_averted_dp_inc <- aggregate_impact_by_dp(
  input_data          = read_csv(file.path("tables", "deaths_susceptible_psa.csv")),
  DiseasePresentation = DiseasePresentation_death,
  mode                = "incremental")

# create table for incremental avertable DALY estimates by disease presentation
# -- (High-potential Scenario - Baseline Scenario)
Associated_daly_averted_dp_inc   <- aggregate_impact_by_dp(
  input_data          = read_csv(file.path("tables", "daly_associated_psa.csv")),
  DiseasePresentation = DiseasePresentation_daly,
  mode                = "incremental")

Attributable_daly_averted_dp_inc <- aggregate_impact_by_dp(
  input_data          = read_csv(file.path("tables", "daly_attributable_psa.csv")),
  DiseasePresentation = DiseasePresentation_daly,
  mode                = "incremental")

Susceptible_daly_averted_dp_inc <- aggregate_impact_by_dp(
  input_data          = read_csv(file.path("tables", "daly_susceptible_psa.csv")),
  DiseasePresentation = DiseasePresentation_daly,
  mode                = "incremental")
# ------------------------------------------------------------------------------
# Table: vaccine avertable AMR health burden by infectious syndrome, 2019

Table_avertable_burden_by_dp <- create_avertable_burden_table(
  Associated_death_averted       = Associated_death_averted_dp,
  Associated_death_averted_opt   = Associated_death_averted_dp_opt,
  Associated_death_averted_inc   = Associated_death_averted_dp_inc, 
  Associated_daly_averted        = Associated_daly_averted_dp,
  Associated_daly_averted_opt    = Associated_daly_averted_dp_opt,
  Associated_daly_averted_inc    = Associated_daly_averted_dp_inc,
  Attributable_death_averted     = Attributable_death_averted_dp,
  Attributable_death_averted_opt = Attributable_death_averted_dp_opt, 
  Attributable_death_averted_inc = Attributable_death_averted_dp_inc, 
  Attributable_daly_averted      = Attributable_daly_averted_dp,
  Attributable_daly_averted_opt  = Attributable_daly_averted_dp_opt,
  Attributable_daly_averted_inc  = Attributable_daly_averted_dp_inc,
  Susceptible_death_averted      = Susceptible_death_averted_dp,
  Susceptible_death_averted_opt  = Susceptible_death_averted_dp_opt,
  Susceptible_death_averted_inc  = Susceptible_death_averted_dp_inc,
  Susceptible_daly_averted       = Susceptible_daly_averted_dp,
  Susceptible_daly_averted_opt   = Susceptible_daly_averted_dp_opt,
  Susceptible_daly_averted_inc   = Susceptible_daly_averted_dp_inc)

Table_avertable_burden_by_dp <- 
  Table_avertable_burden_by_dp %>% arrange(Counts)

fwrite (x    = Table_avertable_burden_by_dp,
        file = file.path("tables", "Table_avertable_burden_by_infectious_syndrome.csv"))

# ------------------------------------------------------------------------------
# Figure: Vaccine impact by infectious syndrome, 2019

# create graph with avertable AMR burden by infectious syndrome
# -- Baseline Scenario
death_averted_by_dp_graph <- create_burden_averted_by_dp_graph(
  Attributable_burden_averted = Attributable_death_averted_dp,
  Associated_burden_averted   = Associated_death_averted_dp,
  ylabel = "Vaccine Avertable Deaths")

daly_averted_by_dp_graph <- create_burden_averted_by_dp_graph(
  Attributable_burden_averted = Attributable_daly_averted_dp,
  Associated_burden_averted   = Associated_daly_averted_dp,
  ylabel = "Vaccine Avertable DALYs")

# create Figure
death_averted_by_dp_graph / daly_averted_by_dp_graph

# save the image file
ggsave (filename = "Figure_avertable_burden_by_dp.png",
        path = "figures",
        width = 6, 
        height = 8, 
        dpi = 600)

ggsave (filename = "Figure_avertable_burden_by_dp.eps",
        path = "figures",
        device = "eps",
        width = 7, 
        height = 7)

# ------------------------------------------------------------------------------
# global vaccine avertable deaths and DALYs attributable to and associated with 
# bacterial antimicrobial resistance by vaccine profiles, 2019

# vaccine profile with multiple options
vaccine_profile_dt_add <- read_csv(file.path("tables", "vaccine_profile.csv"))

# create table for avertable death estimates by vaccine profile
# -- Baseline Scenario
Associated_death_averted_vp <- vaccine_imapct_by_vaccine(
  data_in = read_csv(file.path("tables", "deaths_associated_psa.csv")),
  mode_in = "baseline")

Attributable_death_averted_vp <- vaccine_imapct_by_vaccine(
  data_in = read_csv(file.path("tables", "deaths_attributable_psa.csv")),
  mode_in = "baseline")

Susceptible_death_averted_vp <- vaccine_imapct_by_vaccine(
  data_in = read_csv(file.path("tables", "deaths_susceptible_psa.csv")),
  mode_in = "baseline")

# create table for avertable DALY estimates by vaccine profile
# -- Baseline Scenario
Associated_daly_averted_vp <- vaccine_imapct_by_vaccine(
  data_in = read_csv(file.path("tables", "daly_associated_psa.csv")),
  mode_in = "baseline")

Attributable_daly_averted_vp <- vaccine_imapct_by_vaccine(
  data_in = read_csv(file.path("tables", "daly_attributable_psa.csv")),
  mode_in = "baseline")

Susceptible_daly_averted_vp <- vaccine_imapct_by_vaccine(
  data_in = read_csv(file.path("tables", "daly_susceptible_psa.csv")),
  mode_in = "baseline")

# create table for avertable death estimates by vaccine profile
# -- High-potential Scenario
Associated_death_averted_vp_opt <- vaccine_imapct_by_vaccine(
  data_in = read_csv(file.path("tables", "deaths_associated_psa.csv")),
  mode_in = "high_potential")

Attributable_death_averted_vp_opt <- vaccine_imapct_by_vaccine(
  data_in = read_csv(file.path("tables", "deaths_attributable_psa.csv")),
  mode_in = "high_potential")

Susceptible_death_averted_vp_opt <- vaccine_imapct_by_vaccine(
  data_in = read_csv(file.path("tables", "deaths_susceptible_psa.csv")),
  mode_in = "high_potential")

# create table for avertable DALY estimates by vaccine profile
# -- High-potential Scenario
Associated_daly_averted_vp_opt <- vaccine_imapct_by_vaccine(
  data_in = read_csv(file.path("tables", "daly_associated_psa.csv")),
  mode_in = "high_potential")

Attributable_daly_averted_vp_opt <- vaccine_imapct_by_vaccine(
  data_in = read_csv(file.path("tables", "daly_attributable_psa.csv")),
  mode_in = "high_potential")

Susceptible_daly_averted_vp_opt <- vaccine_imapct_by_vaccine(
  data_in = read_csv(file.path("tables", "daly_susceptible_psa.csv")),
  mode_in = "high_potential")

# create table for avertable incremental death estimates by vaccine profile
# -- (High-potential Scenario - Baseline Scenario)
Associated_death_averted_vp_inc <- vaccine_imapct_by_vaccine(
  data_in = read_csv(file.path("tables", "deaths_associated_psa.csv")),
  mode_in = "incremental")

Attributable_death_averted_vp_inc <- vaccine_imapct_by_vaccine(
  data_in = read_csv(file.path("tables", "deaths_attributable_psa.csv")),
  mode_in = "incremental")

Susceptible_death_averted_vp_inc <- vaccine_imapct_by_vaccine(
  data_in = read_csv(file.path("tables", "deaths_susceptible_psa.csv")),
  mode_in = "incremental")

# create table for avertable incremental DALY estimates by vaccine profile
# -- (High-potential Scenario - Baseline Scenario)
Associated_daly_averted_vp_inc <- vaccine_imapct_by_vaccine(
  data_in = read_csv(file.path("tables", "daly_associated_psa.csv")),
  mode_in = "incremental")

Attributable_daly_averted_vp_inc <- vaccine_imapct_by_vaccine(
  data_in = read_csv(file.path("tables", "daly_attributable_psa.csv")),
  mode_in = "incremental")

Susceptible_daly_averted_vp_inc <- vaccine_imapct_by_vaccine(
  data_in = read_csv(file.path("tables", "daly_susceptible_psa.csv")),
  mode_in = "incremental")
# ------------------------------------------------------------------------------
# Table: vaccine avertable AMR health burden by vaccine profile, 2019

Table_avertable_burden_by_vp <- create_avertable_burden_table(
  Associated_death_averted       = Associated_death_averted_vp,
  Associated_death_averted_opt   = Associated_death_averted_vp_opt,
  Associated_death_averted_inc   = Associated_death_averted_vp_inc,
  Associated_daly_averted        = Associated_daly_averted_vp,
  Associated_daly_averted_opt    = Associated_daly_averted_vp_opt,
  Associated_daly_averted_inc    = Associated_daly_averted_vp_inc,
  Attributable_death_averted     = Attributable_death_averted_vp,
  Attributable_death_averted_opt = Attributable_death_averted_vp_opt, 
  Attributable_death_averted_inc = Attributable_death_averted_vp_inc, 
  Attributable_daly_averted      = Attributable_daly_averted_vp,
  Attributable_daly_averted_opt  = Attributable_daly_averted_vp_opt,
  Attributable_daly_averted_inc  = Attributable_daly_averted_vp_inc,
  Susceptible_death_averted      = Susceptible_death_averted_vp,
  Susceptible_death_averted_opt  = Susceptible_death_averted_vp_opt,
  Susceptible_death_averted_inc  = Susceptible_death_averted_vp_inc,
  Susceptible_daly_averted       = Susceptible_daly_averted_vp,
  Susceptible_daly_averted_opt   = Susceptible_daly_averted_vp_opt,
  Susceptible_daly_averted_inc   = Susceptible_daly_averted_vp_inc)

Table_avertable_burden_by_vp <- left_join(vaccine_profile_dt_add, 
                                          Table_avertable_burden_by_vp,
                                          by = c("Vaccine_pathogen" = "Counts"))

fwrite (x    = Table_avertable_burden_by_vp,
        file = file.path("tables", "Table_avertable_burden_by_vaccine_profile.csv"))

# ------------------------------------------------------------------------------
# Figure: vaccine impact by vaccine profile, 2019

# create graph with avertable AMR burden by vaccine profile -- Baseline Scenario
death_averted_by_vp_graph <- create_burden_averted_by_vp_graph(
  Attributable_burden_averted = Attributable_death_averted_vp,
  Associated_burden_averted   = Associated_death_averted_vp,
  ylabel = "Vaccine Avertable Deaths",
  ylim_max = 150000)

daly_averted_by_vp_graph <- create_burden_averted_by_vp_graph(
  Attributable_burden_averted = Attributable_daly_averted_vp,
  Associated_burden_averted   = Associated_daly_averted_vp,
  ylabel = "Vaccine Avertable DALYs",
  ylim_max = 11000000)

# create Figure
death_averted_by_vp_graph / daly_averted_by_vp_graph

# save the image file
ggsave (filename = "Figure_burden_averted_by_vaccine_profile.png",
        path = "figures",
        width = 6, 
        height = 9,
        dpi = 600)

ggsave (filename = "Figure_avertable_burden_by_vaccine_profile.eps",
        path = "figures",
        device = "eps",
        width = 6, 
        height = 10)

# ------------------------------------------------------------------------------



# ------------------------------------------------------------------------------
# PCV and Hib vaccine impact with current coverage
# -- baseline scenario

vaccine_impact_current <- bind_rows(list(
# impact of current coverage of Hib vaccine on deaths associated with AMR
edit_table(estimate_existing_vaccine_impact(
    input_data = read_csv(file.path("tables", "deaths_associated_psa.csv")))[Counts == "Haemophilus influenzae",]),

# impact of current coverage of Hib vaccine on deaths attributable to AMR
edit_table(estimate_existing_vaccine_impact(
    input_data = read_csv(file.path("tables", "deaths_attributable_psa.csv")))[Counts == "Haemophilus influenzae",]),

# impact of current coverage of Hib vaccine on dalys associated with AMR
edit_table(estimate_existing_vaccine_impact(
    input_data = read_csv(file.path("tables", "daly_associated_psa.csv")))[Counts == "Haemophilus influenzae",]),

# impact of current coverage of Hib vaccine on dalys attributable to AMR
edit_table(estimate_existing_vaccine_impact(
    input_data = read_csv(file.path("tables", "daly_attributable_psa.csv")))[Counts == "Haemophilus influenzae",]),

# impact of current coverage of PCV on deaths associated with AMR
edit_table(estimate_existing_vaccine_impact(
    input_data = read_csv(file.path("tables", "deaths_associated_psa.csv")))[Counts == "Streptococcus pneumoniae",]),

# impact of current coverage of PCV on deaths attributable to AMR
edit_table(estimate_existing_vaccine_impact(
    input_data = read_csv(file.path("tables", "deaths_attributable_psa.csv")))[Counts == "Streptococcus pneumoniae",]),

# impact of current coverage of PCV on dalys associated with AMR
edit_table(estimate_existing_vaccine_impact(
    input_data = read_csv(file.path("tables", "daly_associated_psa.csv")))[Counts == "Streptococcus pneumoniae",]),

# impact of current coverage of PCV on dalys attributable to AMR
edit_table(estimate_existing_vaccine_impact(
    input_data = read_csv(file.path("tables", "daly_attributable_psa.csv")))[Counts == "Streptococcus pneumoniae",])))

vaccine_impact_current$type <- rep(c("death_associated", "death_attributable", "daly_associated", "daly_attributable"), 2)

fwrite (x    = vaccine_impact_current,
        file = file.path("tables", "Table_avertable_burden_with_current_coverage.csv"))

# ------------------------------------------------------------------------------
# Vaccine avertable deaths by infectious syndrome and pathogen

deaths_associated_dp_pathogen   <- 
  aggregate_impact_by_dp_pathogen(input_data = read_csv(file.path("tables", "deaths_associated_psa.csv")),
                                  input_rep  = 1:62)

deaths_attributable_dp_pathogen <- 
  aggregate_impact_by_dp_pathogen(input_data = read_csv(file.path("tables", "deaths_attributable_psa.csv")),
                                  input_rep  = 1:62)

daly_associated_dp_pathogen     <- 
  aggregate_impact_by_dp_pathogen(input_data = read_csv(file.path("tables", "daly_associated_psa.csv")),
                                  input_rep  = 1:63)

daly_attributable_dp_pathogen   <- 
  aggregate_impact_by_dp_pathogen(input_data = read_csv(file.path("tables", "daly_attributable_psa.csv")),
                                  input_rep  = 1:61)

# create graph of vaccine impact by infectious syndrome and pathogen

burden_averted_by_dp_pat(data_input  = deaths_associated_dp_pathogen, 
                         start_input = 3.6,
                         image_png   = file.path("figures", 
                                                 "Figure_deaths_associated_dp_pat.png"),
                         image_eps   = file.path("figures", 
                                                 "Figure_deaths_associated_dp_pat.eps"))

burden_averted_by_dp_pat(data_input  = deaths_attributable_dp_pathogen,
                         start_input = 3.6,
                         image_png   = file.path("figures", 
                                                 "Figure_deaths_attributable_dp_pat.png"),
                         image_eps   = file.path("figures", 
                                                 "Figure_deaths_attributable_dp_pat.eps"))

burden_averted_by_dp_pat(data_input  = daly_associated_dp_pathogen, 
                         start_input = 4.1,
                         image_png   = file.path("figures", 
                                                 "Figure_daly_associated_dp_pat.png"),
                         image_eps   = file.path("figures", 
                                                 "Figure_daly_associated_dp_pat.eps"))

burden_averted_by_dp_pat(data_input  = daly_attributable_dp_pathogen,
                         start_input = 3.9,
                         image_png   = file.path("figures", 
                                                "Figure_daly_attributable_dp_pat.png"),
                         image_eps   = file.path("figures", 
                                                "Figure_daly_attributable_dp_pat.eps"))

# ------------------------------------------------------------------------------



# ------------------------------------------------------------------------------
# [Table in appendix] vccine avertable burdens by region, infectious syndrome, and pathogen

avertable_death_associated_combination   <- 
  aggregate_impact_by_region_dp_pathogen(input_data = read_csv(file.path("tables", "deaths_associated_psa.csv")),
                                         input_rep  = 1:372,
                                         file_name = file.path ("tables", "Table_avertable_death_associated_combination.csv"))

avertable_death_attributable_combination <- 
  aggregate_impact_by_region_dp_pathogen(input_data = read_csv(file.path("tables", "deaths_attributable_psa.csv")),
                                         input_rep  = 1:372,
                                         file_name = file.path ("tables", "Table_avertable_death_attributable_combination.csv"))

avertable_daly_associated_combination     <- 
  aggregate_impact_by_region_dp_pathogen(input_data = read_csv(file.path("tables", "daly_associated_psa.csv")),
                                         input_rep  = 1:357,
                                         file_name = file.path ("tables", "Table_avertable_daly_associated_combination.csv"))

avertable_daly_attributable_combination   <- 
  aggregate_impact_by_region_dp_pathogen(input_data = read_csv(file.path("tables", "daly_attributable_psa.csv")),
                                         input_rep  = 1:322,
                                         file_name = file.path ("tables", "Table_avertable_daly_attributable_combination.csv"))

# ------------------------------------------------------------------------------
# [Table in appendix] vaccine avertable health burdens associated with and attributable
# to AMR by WHO region, pathogen, disease presentation, and age group

# vaccine avertable AMR burden estimates using the burden estimate median value

# create vaccine avertable burden file -- Baseline Scenario
AMR_death_data_updated <- update_death_burden(
  input_susceptible            = read_csv(file.path("tables", "susceptible_burden.csv")),
  input_associated             = read_csv(file.path("tables", "associated_burden.csv")),
  input_attributable           = read_csv(file.path("tables", "attributable_burden.csv")),
  burden_dt                    = death_burden_dt,
  AMR_burden_data_updated_file = file.path ("tables", "AMR_death_data_baseline.csv"),
  input_scenario               = "conservative")

AMR_daly_data_updated  <- update_death_burden(
  input_susceptible            = read_csv(file.path("tables", "daly_susceptible_burden.csv")),
  input_associated             = read_csv(file.path("tables", "daly_associated_burden.csv")),
  input_attributable           = read_csv(file.path("tables", "daly_attributable_burden.csv")),
  burden_dt                    = daly_burden_dt,
  AMR_burden_data_updated_file = file.path ("tables", "AMR_daly_data_baseline.csv"),
  input_scenario               = "conservative")

# create vaccine avertable burden file -- High-potential Scenario
AMR_death_data_updated  <- update_death_burden(
  input_susceptible            = read_csv(file.path("tables", "susceptible_burden.csv")),
  input_associated             = read_csv(file.path("tables", "associated_burden.csv")),
  input_attributable           = read_csv(file.path("tables", "attributable_burden.csv")),
  burden_dt                    = death_burden_dt,
  AMR_burden_data_updated_file = file.path ("tables", "AMR_death_data_high_potential.csv"),
  input_scenario               = "optimistic")

AMR_daly_data_updated <- update_death_burden(
  input_susceptible            = read_csv(file.path("tables", "daly_susceptible_burden.csv")),
  input_associated             = read_csv(file.path("tables", "daly_associated_burden.csv")),
  input_attributable           = read_csv(file.path("tables", "daly_attributable_burden.csv")),
  burden_dt                    = daly_burden_dt,
  AMR_burden_data_updated_file = file.path ("tables", "AMR_daly_data_high_potential.csv"),
  input_scenario               = "optimistic")
# ------------------------------------------------------------------------------



# ------------------------------------------------------------------------------
# return to source directory
setwd (source_wd)

# end time
end_time <- Sys.time ()
print (paste0 ("end time = ", end_time))
print (Sys.time () - start_time)
# ------------------------------------------------------------------------------
