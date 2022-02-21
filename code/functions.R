# ------------------------------------------------------------------------------
# functions.R
#
# functions for analysis to estimate vaccine averted AMR health burden
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# 2019 vaccine coverage for existing vaccine: Hib vaccine, PCV
existing_vaccine_coverage <- function(hib_coverage_file, pcv_coverage_file) {
  
  
  # "WPP2019_INT_F03_1_POPULATION_BY_AGE_ANNUAL_BOTH_SEXES.xlsx" data cleaning
  population_file <- read_excel("data/WPP2019 population by age.xlsx", 
                                col_names = FALSE)
  
  names(population_file) <- lapply(population_file[13, ], as.character)
  population_file <- population_file %>% filter(Type=="Country/Area")
  
  # Removing duplicates (year of 2020)
  population_file <- population_file[population_file$`Reference date (as of 1 July)` == "2019", 
                                     c("Region, subregion, country or area *", 
                                       "Reference date (as of 1 July)", "0")]
  
  # ----------------------------------------------------------------------------
  # "WPP2019_F01_LOCATIONS.XLSX" data cleaning
  WPP_iso3 <- read_excel("data/WPP2019 locations.XLSX", 
                         col_names = FALSE)
  names(WPP_iso3) <- lapply(WPP_iso3[13, ], as.character)
  WPP_iso3 <- WPP_iso3[-c(1:13),c(2,5)]
  
  # ----------------------------------------------------------------------------
  # adding ISO3 to WPP2019_POPULATION data frame
  population_iso3 <- left_join(population_file, WPP_iso3, 
                               by=c("Region, subregion, country or area *"
                                    = "Region, subregion, country or area*"))
  
  # changing population to the actual value
  population_iso3$`0`<- as.numeric(population_iso3$`0`) * 1000
  
  population_iso3 <- rename(population_iso3, "Population" = "0")
  
  # ----------------------------------------------------------------------------
  # Import WUENIC vaccine coverage data for HIB vaccine and PCV
  WUENIC <- read_excel("data/WUENIC.xlsx", sheet = "WUENIC_input_to_PDF")
  WUENIC <- WUENIC[, c("ISOCountryCode", "Year", "Vaccine", "WUENIC")]
  
  WUENIC <- WUENIC %>% filter(Year == "2019")
  
  WUENIC$"ISOCountryCode" = toupper(WUENIC$"ISOCountryCode")
  
  WUENIC_hib <- WUENIC %>% filter(Year == "2019" & Vaccine == "hib3")
  WUENIC_pcv <- WUENIC %>% filter(Year == "2019" & Vaccine == "pcv3")
  
  # ----------------------------------------------------------------------------
  # Import WHO region classification data
  country_region <- read_csv("data/country_income_region_classification.csv")
  country_region <- country_region[, c("iso3_code", "WHO_region")]
  
  # ----------------------------------------------------------------------------
  # Combine population data with current vaccine coverage data
  WUENIC_WPP_hib <- left_join(WUENIC_hib, population_iso3, by=c("ISOCountryCode" = "ISO3 Alpha-code"))
  WUENIC_WPP_pcv <- left_join(WUENIC_pcv, population_iso3, by=c("ISOCountryCode" = "ISO3 Alpha-code"))
  
  # Add region
  WUENIC_WPP_hib <- left_join(WUENIC_WPP_hib, country_region, by=c("ISOCountryCode" = "iso3_code"))
  WUENIC_WPP_pcv <- left_join(WUENIC_WPP_pcv, country_region, by=c("ISOCountryCode" = "iso3_code"))
  
  # Remove NA
  # population data is missiong for some countries / Palestine is not a member state of WHO
  WUENIC_WPP_hib <- na.omit(WUENIC_WPP_hib)
  WUENIC_WPP_pcv <- na.omit(WUENIC_WPP_pcv)
  
  # ----------------------------------------------------------------------------
  # Estimate regional coverage - Hib vaccine
  WUENIC_WPP_hib$vaccinated_pop <- WUENIC_WPP_hib$WUENIC * 1/100 * WUENIC_WPP_hib$Population 
  
  hib_vaccinated <- WUENIC_WPP_hib %>% 
    group_by(WHO_region) %>%
    summarise(vaccinated_pop = sum(vaccinated_pop))
  
  hib_total_pop <- WUENIC_WPP_hib %>% 
    group_by(WHO_region) %>%
    summarise(total_pop = sum(Population))
  
  hib_vaccinated$total_pop <- hib_total_pop$total_pop
  
  hib_vaccinated$vaccine_coverage <- hib_vaccinated$vaccinated_pop / hib_vaccinated$total_pop
  
  # Estimate regional coverage - PCV
  WUENIC_WPP_pcv$vaccinated_pop <- WUENIC_WPP_pcv$WUENIC * 1/100 * WUENIC_WPP_pcv$Population 
  
  pcv_vaccinated <- WUENIC_WPP_pcv %>% 
    group_by(WHO_region) %>%
    summarise(vaccinated_pop = sum(vaccinated_pop))
  
  pcv_total_pop <- WUENIC_WPP_pcv %>% 
    group_by(WHO_region) %>%
    summarise(total_pop = sum(Population))
  
  pcv_vaccinated$total_pop <- pcv_total_pop$total_pop
  
  pcv_vaccinated$vaccine_coverage <- pcv_vaccinated$vaccinated_pop / pcv_vaccinated$total_pop
  
  # save regional vaccine coverage for existing vaccine
  fwrite (x    = hib_vaccinated, 
          file = hib_coverage_file)
  
  fwrite (x    = pcv_vaccinated, 
          file = pcv_coverage_file)
  
  return ()
  
} # end of function -- existing_vaccine_coverage
# ------------------------------------------------------------------------------



# ------------------------------------------------------------------------------
# create and clean up the data frame of AMR burden

create_burden_table <- function(AMR_burden,
                                burden_file){
  
  AMR_burden <- AMR_burden[-c(1:2), -c(5:7)]
  
  names(AMR_burden) <- c("WHO_region", "Disease_presentation","Age_group","Pathogen", 
                         "Associated_resistant_mean", "Associated_resistant_lower",
                         "Associated_resistant_upper", "Attributable_resistance_mean",
                         "Attributable_resistance_lower", "Attributable_resistance_upper")
  
  AMR_burden[,5:10] <- lapply(AMR_burden[,5:10], function(x) gsub(",", "", x))
  
  AMR_burden[,5:10] <- lapply(AMR_burden[,5:10], as.numeric)
  
  AMR_burden$Age_group <- factor(AMR_burden$Age_group, 
                                 levels=unique(AMR_burden$Age_group), order=T)
  
  levels(AMR_burden$Age_group) <- c("EN", "LN", "PN",  "1 to 4", "5 to 9","10 to 14",
                                    "15 to 19", "20 to 24"  ,"25 to 29",  "30 to 34",
                                    "35 to 39", "40 to 44", "45 to 49",  "50 to 54" ,      
                                    "55 to 59", "60 to 64", "65 to 69", "70 to 74",
                                    "75 to 79", "80 to 84", "85 to 89", "90 to 94", "95 plus")
  
  AMR_burden <- AMR_burden %>%
    filter(AMR_burden$WHO_region != "unclassified")
  
# ------------------------------------------------------------------------------  
  # Estimate pre-vaccine burden for existing vaccines: HIB vaccine & PCV
  
  # Import current vaccine coverage estimates by region  
  hib_coverage <- read_csv("tables/hib coverage.csv")
  
  pcv_coverage <- read_csv("tables/pcv coverage.csv")
  
  hib_coverage <- data.table (WHO_region = c("Africa", "Americas", "Eastern Mediterranean", 
                                             "Europe", "South-East Asia", "Western Pacific"), 
                              hib_vaccine_coverage = hib_coverage$vaccine_coverage)
  
  pcv_coverage <- data.table (WHO_region = c("Africa", "Americas", "Eastern Mediterranean", 
                                             "Europe", "South-East Asia", "Western Pacific"), 
                              pcv_vaccine_coverage = pcv_coverage$vaccine_coverage)
  
  # Substitute the current burden estimates to pre-vaccine burden estimates for existing vaccines
  AMR_burden <- left_join(AMR_burden, hib_coverage, by=c("WHO_region" = "WHO_region"))
  
  AMR_burden <- left_join(AMR_burden, pcv_coverage, by=c("WHO_region" = "WHO_region"))
  
  AMR_burden <- data.table(AMR_burden)
  
  # apply prevaccine burden to HIB burden
  AMR_burden[Pathogen  == "Haemophilus influenzae" & Age_group == "PN",
             Associated_resistant_mean := Associated_resistant_mean /
               (3/48 + 4/48 * (1 - hib_vaccine_coverage *0.59) + 
                  4/48 * (1 - hib_vaccine_coverage * 0.92) + 
                  37/48 * (1 - hib_vaccine_coverage * 0.93))]
  
  AMR_burden[Pathogen  == "Haemophilus influenzae" & Age_group == "PN",
             Associated_resistant_lower := Associated_resistant_lower /
               (3/48 + 4/48 * (1 - hib_vaccine_coverage *0.59) + 
                  4/48 * (1 - hib_vaccine_coverage * 0.92) + 
                  37/48 * (1 - hib_vaccine_coverage * 0.93))]
  
  AMR_burden[Pathogen  == "Haemophilus influenzae" & Age_group == "PN",
             Associated_resistant_upper := Associated_resistant_upper /
               (3/48 + 4/48 * (1 - hib_vaccine_coverage *0.59) + 
                  4/48 * (1 - hib_vaccine_coverage * 0.92) + 
                  37/48 * (1 - hib_vaccine_coverage * 0.93))]
  
  AMR_burden[Pathogen  == "Haemophilus influenzae" & Age_group == "PN",
             Attributable_resistance_mean := Attributable_resistance_mean /
               (3/48 + 4/48 * (1 - hib_vaccine_coverage *0.59) + 
                  4/48 * (1 - hib_vaccine_coverage * 0.92) + 
                  37/48 * (1 - hib_vaccine_coverage * 0.93))]
  
  AMR_burden[Pathogen  == "Haemophilus influenzae" & Age_group == "PN",
             Attributable_resistance_lower := Attributable_resistance_lower /
               (3/48 + 4/48 * (1 - hib_vaccine_coverage *0.59) + 
                  4/48 * (1 - hib_vaccine_coverage * 0.92) + 
                  37/48 * (1 - hib_vaccine_coverage * 0.93))]
  
  AMR_burden[Pathogen  == "Haemophilus influenzae" & Age_group == "PN",
             Attributable_resistance_upper := Attributable_resistance_upper /
               (3/48 + 4/48 * (1 - hib_vaccine_coverage *0.59) + 
                  4/48 * (1 - hib_vaccine_coverage * 0.92) + 
                  37/48 * (1 - hib_vaccine_coverage * 0.93))]
  
  
  AMR_burden[Pathogen  == "Haemophilus influenzae" & Age_group == "1 to 4",
             Associated_resistant_mean := Associated_resistant_mean /
               (1 - hib_vaccine_coverage * 0.59)]
  
  AMR_burden[Pathogen  == "Haemophilus influenzae" & Age_group == "1 to 4",
             Associated_resistant_lower := Associated_resistant_lower /
               (1 - hib_vaccine_coverage * 0.59)]
  
  AMR_burden[Pathogen  == "Haemophilus influenzae" & Age_group == "1 to 4",
             Associated_resistant_upper := Associated_resistant_upper /
               (1 - hib_vaccine_coverage * 0.59)]
  
  AMR_burden[Pathogen  == "Haemophilus influenzae" & Age_group == "1 to 4",
             Attributable_resistance_mean := Attributable_resistance_mean /
               (1 - hib_vaccine_coverage * 0.59)]
  
  AMR_burden[Pathogen  == "Haemophilus influenzae" & Age_group == "1 to 4",
             Attributable_resistance_lower := Attributable_resistance_lower /
               (1 - hib_vaccine_coverage * 0.59)]
  
  AMR_burden[Pathogen  == "Haemophilus influenzae" & Age_group == "1 to 4",
             Attributable_resistance_upper := Attributable_resistance_upper /
               (1 - hib_vaccine_coverage * 0.59)]
  
  
  AMR_burden[Pathogen  == "Streptococcus pneumoniae" & Age_group == "PN",
             Associated_resistant_mean := Associated_resistant_mean /
               (3/48 + 4/48 * (1 - hib_vaccine_coverage *0.29) + 
                  4/48 * (1 - pcv_vaccine_coverage * 0.58) + 
                  37/48 * (1 - pcv_vaccine_coverage * 0.58))]
  
  AMR_burden[Pathogen  == "Streptococcus pneumoniae" & Age_group == "PN",
             Associated_resistant_lower := Associated_resistant_lower /
               (3/48 + 4/48 * (1 - hib_vaccine_coverage *0.29) + 
                  4/48 * (1 - pcv_vaccine_coverage * 0.58) + 
                  37/48 * (1 - pcv_vaccine_coverage * 0.58))]
  
  AMR_burden[Pathogen  == "Streptococcus pneumoniae" & Age_group == "PN",
             Associated_resistant_upper := Associated_resistant_upper /
               (3/48 + 4/48 * (1 - hib_vaccine_coverage *0.29) + 
                  4/48 * (1 - pcv_vaccine_coverage * 0.58) + 
                  37/48 * (1 - pcv_vaccine_coverage * 0.58))]
  
  AMR_burden[Pathogen  == "Streptococcus pneumoniae" & Age_group == "PN",
             Attributable_resistance_mean := Attributable_resistance_mean /
               (3/48 + 4/48 * (1 - hib_vaccine_coverage *0.29) + 
                  4/48 * (1 - pcv_vaccine_coverage * 0.58) + 
                  37/48 * (1 - pcv_vaccine_coverage * 0.58))]
  
  AMR_burden[Pathogen  == "Streptococcus pneumoniae" & Age_group == "PN",
             Attributable_resistance_lower := Attributable_resistance_lower /
               (3/48 + 4/48 * (1 - hib_vaccine_coverage *0.29) + 
                  4/48 * (1 - pcv_vaccine_coverage * 0.58) + 
                  37/48 * (1 - pcv_vaccine_coverage * 0.58))]
  
  AMR_burden[Pathogen  == "Streptococcus pneumoniae" & Age_group == "PN",
             Attributable_resistance_upper := Attributable_resistance_upper /
               (3/48 + 4/48 * (1 - hib_vaccine_coverage *0.29) + 
                  4/48 * (1 - pcv_vaccine_coverage * 0.58) + 
                  37/48 * (1 - pcv_vaccine_coverage * 0.58))]
  
  
  AMR_burden[Pathogen  == "Streptococcus pneumoniae" & Age_group == "1 to 4",
             Associated_resistant_mean := Associated_resistant_mean /
               (1 - pcv_vaccine_coverage * 0.58)]
  
  AMR_burden[Pathogen  == "Streptococcus pneumoniae" & Age_group == "1 to 4",
             Associated_resistant_lower := Associated_resistant_lower /
               (1 - pcv_vaccine_coverage * 0.58)]
  
  AMR_burden[Pathogen  == "Streptococcus pneumoniae" & Age_group == "1 to 4",
             Associated_resistant_upper := Associated_resistant_upper /
               (1 - pcv_vaccine_coverage * 0.58)]
  
  AMR_burden[Pathogen  == "Streptococcus pneumoniae" & Age_group == "1 to 4",
             Attributable_resistance_mean := Attributable_resistance_mean /
               (1 - pcv_vaccine_coverage * 0.58)]
  
  AMR_burden[Pathogen  == "Streptococcus pneumoniae" & Age_group == "1 to 4",
             Attributable_resistance_lower := Attributable_resistance_lower /
               (1 - pcv_vaccine_coverage * 0.58)]
  
  AMR_burden[Pathogen  == "Streptococcus pneumoniae" & Age_group == "1 to 4",
             Attributable_resistance_upper := Attributable_resistance_upper /
               (1 - pcv_vaccine_coverage * 0.58)]
  
  fwrite (x    = AMR_burden,
          file = burden_file)
  
  return(AMR_burden)
} # end of function -- create_death_burden_table

# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# create and clean up the data frame of WHO vaccine profile

create_vaccine_profile_table <- function(vaccine_profile,
                                         vaccine_profile_file){
  
  vaccine_profile <- vaccine_profile[, c(2, 3, 4, 5, 9, 10, 11, 12)]
  
  vaccine_profile <- vaccine_profile %>%
    rename("Efficacy"     = "Efficacy (%)",
           "Coverage"     = "Coverage in target group",
           "Duration"     = "Duration of protection", 
           "DP"           = "Disease presentation",
           "VC"           = "Vaccination (conservative scenario - 1 time vaccination)",
           "VO"           = "Vaccination (optimistic scenario - plus repeat/booster vaccinations)",
           "MainAnalysis" = "Main analysis")
  
  vaccine_profile$Efficacy <- vaccine_profile$Efficacy * 1/100
  vaccine_profile$Coverage <- vaccine_profile$Coverage * 1/100
  
  fwrite (x    = vaccine_profile, 
          file = vaccine_profile_file)
  
  vaccine_profile <- vaccine_profile %>% filter(MainAnalysis == "Yes")
  
  return(vaccine_profile)
} # end of function -- create_vaccine_profile_table

# ------------------------------------------------------------------------------
# create combined table: disease burden + vaccine profile

create_combined_table <- function(death_burden_dt, 
                                  vaccine_profile_dt,
                                  attributable_burden_file,
                                  associated_burden_file){
  
# create combined table
  combined_table <- left_join(death_burden_dt, vaccine_profile_dt, 
                              by=c("Pathogen" = "Pathogen"))
  combined_table <- data.table(combined_table)
  
# Streptococcus pneumoniae vaccine specification
  combined_table[Pathogen             == "Streptococcus pneumoniae" &
                 Disease_presentation == "LRI and thorax infections",
                 Efficacy             := 0.5]
  
# separate burden attributable to AMR and associated with AMR  
  
  attributable_burden <- combined_table[, -c("Associated_resistant_mean",
                                              "Associated_resistant_lower",
                                              "Associated_resistant_upper")]
  
  attributable_burden <- attributable_burden %>% 
    rename("burden_lower_value" = "Attributable_resistance_lower",
           "burden_mean_value"  = "Attributable_resistance_mean",
           "burden_upper_value" = "Attributable_resistance_upper")
  
   fwrite (x    = attributable_burden,
           file = attributable_burden_file)
  
  associated_burden <- combined_table[, -c("Attributable_resistance_mean",
                                           "Attributable_resistance_lower",
                                           "Attributable_resistance_upper")]
  
  associated_burden <- associated_burden %>%
    rename("burden_lower_value" = "Associated_resistant_lower",
           "burden_mean_value"  = "Associated_resistant_mean",
           "burden_upper_value" = "Associated_resistant_upper")
  
  fwrite (x    = associated_burden,
          file = associated_burden_file)
  
  return(combined_table)
  
} # end of function -- create_combined_table

# ------------------------------------------------------------------------------
# create graph of death trend by pathogen across all age groups

create_death_by_pathogen_graph <- function(pathogen){
  dt <- death_burden_dt[death_burden_dt$Pathogen == pathogen, ]
  
  dt <- left_join(dt, vaccine_profile_dt, by=c("Pathogen" = "Pathogen"))
  
  dt <- dt %>% 
    group_by(Pathogen, Efficacy, Coverage, Age_group) %>%
    summarise(sum_Attributable_resistance_mean=sum(Attributable_resistance_mean), .groups = 'drop') 
  
  dt <- dt %>% 
    mutate(v_Attributable_resistance_mean = sum_Attributable_resistance_mean * Efficacy * Coverage,
           va_Attributable_resistance_mean = v_Attributable_resistance_mean-sum_Attributable_resistance_mean)
  
  ggplot(dt, aes(x=Age_group)) +
    geom_line(aes(y=sum_Attributable_resistance_mean, group=1, colour="non-vaccinated")) +
    ylab("Number of Vaccine Attributable Deaths") + 
    xlab("Age") +
    ggtitle(paste(pathogen,"(global)")) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5)) +
    expand_limits(y=0)
  
  # save plot
  ggsave (filename = paste(pathogen,"_global.png"),
          path = "figures",
          width = 15, 
          height = 6, 
          dpi = 600)
} # end of function -- create_death_by_pathogen_graph

# ------------------------------------------------------------------------------



# ------------------------------------------------------------------------------
# log-normal distribution
 uncertainty_analysis_baseline <- function(psa, 
                                           tolerance,
                                           data){
  
  burden_dt <- data.table(data)
  
  # minor changes to fit lognormal distribution
  
  burden_dt <- burden_dt[!(burden_dt$burden_mean_value == "0" & 
                           burden_dt$burden_lower_value == "0" &  
                           burden_dt$burden_upper_value == "0"),]
  
  # add a small 0.25% value to higher bounds where higher bound equals mid value
  burden_dt [burden_mean_value == burden_upper_value, 
             burden_upper_value := burden_upper_value * 1.0025]
  
  # + 0.5 for the value which is 0 (attributable to resistance)
  
  burden_dt$burden_lower_value <- burden_dt$burden_lower_value + 0.5
  
  burden_dt$burden_mean_value  <- burden_dt$burden_mean_value + 0.5
  
  burden_dt$burden_upper_value <- burden_dt$burden_upper_value + 0.5
  
  # ----------------------------------------------------------------------------
  # estimate mean log & sd log (attributable to resistance)
  burden_dt [, burden_sd_log := suppressMessages (suppressWarnings (get.lnorm.par (p           = c (0.025, 0.5, 0.975),
                                                                                   q           = c (burden_lower_value, burden_mean_value, burden_upper_value),
                                                                                   show.output = FALSE,
                                                                                   plot        = FALSE,
                                                                                   tol         = tolerance) ["sdlog"]) ),
             by = .(WHO_region,Disease_presentation,Age_group,Pathogen)]
  
  burden_dt [, burden_mean_log := suppressMessages (suppressWarnings (get.lnorm.par (p           = c (0.025, 0.5, 0.975),
                                                                                     q           = c (burden_lower_value, burden_mean_value, burden_upper_value),
                                                                                     show.output = FALSE,
                                                                                     plot        = FALSE,
                                                                                     tol         = tolerance) ["meanlog"]) ),
             by = .(WHO_region,Disease_presentation,Age_group,Pathogen)]
  
  
  # ----------------------------------------------------------------------------
  # create the table for psa
  
  burden_dt [, run_id := 0]
  
  vaccine_impact_psa <- data.table(burden_dt [0, ])
  
  for (i in 1:psa) {
    # copy data table
    dt <- copy (burden_dt)
    # set run id
    dt [, run_id := i]
    # add combined different scenarios
    vaccine_impact_psa <- rbindlist (list (vaccine_impact_psa, dt),
                                     use.names = TRUE)
  }
# ------------------------------------------------------------------------------
  # rlnorm will generate 'psa' number of values (attributable to resistance)
  
  vaccine_impact_psa [, burden_psa :=
                        rlnorm ( n       = psa,
                                 meanlog = mean (burden_mean_log),
                                 sdlog   = mean (burden_sd_log)),
                      by = .(WHO_region,Disease_presentation,Age_group,Pathogen)]
  
  
  return(vaccine_impact_psa)
  
} # end of function -- uncertainty_analysis_baseline

# ------------------------------------------------------------------------------
# applying vaccine impact on vaccine target population

# applying vaccine target age group
  estimate_vaccine_impact <- function(i, data, burden_psa = burden_psa,
                                      scenario = "conservative"){
  
  if(is.numeric(i) == TRUE) {
    vaccine_age <- data %>%
      filter(run_id == i)
  } else {
    vaccine_age <- data
  }
  
  if((scenario == "optimistic") == TRUE){
  vaccine_age <- vaccine_age %>%
    mutate(va_age =

             ifelse(((VC=="4 week (effective at 6 week)" & VO=="under 5") |
                     (VC=="4 week (effective at 6 week), 60 year" &
                      VO == "under 5, 60 years and above")) & 
                     Age_group == "PN", 
                     burden_psa * Efficacy * Coverage * 47/48,
             ifelse(((VC=="4 week (effective at 6 week)" & VO=="under 5") |
                     (VC=="4 week (effective at 6 week), 60 year" &
                      VO == "under 5, 60 years and above")) &
                     Age_group == "1 to 4", 
                     burden_psa * Efficacy * Coverage,
             ifelse(((VC=="4 week (effective at 6 week)" & VO=="under 5") |
                     (VC=="4 week (effective at 6 week), 60 year" &
                      VO == "under 5, 60 years and above")) & 
                     Duration=="5 years" & Age_group == "5 to 9", 
                     burden_psa * Efficacy * Coverage * 1/5 *5/52,
                                                
                                                              
             ifelse(VC=="10 week (effective at 12 week)" & VO=="under 5" &
                    Age_group == "PN",
                    burden_psa * Efficacy * Coverage * 41/48,
             ifelse(VC=="10 week (effective at 12 week)" & VO=="under 5" & 
                    Age_group == "1 to 4", 
                    burden_psa * Efficacy * Coverage,
             ifelse(VC=="10 week (effective at 12 week)" & VO=="under 5" & 
                    Age_group == "5 to 9", 
                    burden_psa * Efficacy * Coverage * 1/5 * 11/52,

             ifelse(VC=="9 month" & VO=="below 35" & 
                    Age_group == "PN",
                    burden_psa * Efficacy * Coverage * 4/11,
             ifelse(VC=="9 month" & VO=="below 35" & 
                   (Age_group == "1 to 4"   | Age_group == "5 to 9"|
                    Age_group == "10 to 14" | Age_group == "15 to 19"|
                    Age_group == "20 to 24" | Age_group == "25 to 29"|
                    Age_group == "30 to 34"), 
                    burden_psa * Efficacy * Coverage,
                    
             ifelse((VO=="above 10") & 
                   (Age_group == "10 to 14"| Age_group == "14 to 19"|
                    Age_group == "20 to 24"| Age_group == "25 to 29"|
                    Age_group == "30 to 34"| Age_group == "35 to 39"| 
                    Age_group == "40 to 44"| Age_group == "45 to 49"|
                    Age_group == "50 to 54"| Age_group == "55 to 59"),
                    burden_psa * Efficacy * Coverage,                         
                                                                                                       
                                                                                  
             ifelse((VO=="60 years and above"|VO=="under 5, 60 years and above"|
                     VO=="above 10") & 
                   (Age_group == "60 to 64"| Age_group == "65 to 69"|
                    Age_group == "70 to 74"| Age_group == "75 to 79"| 
                    Age_group == "80 to 84"| Age_group == "85 to 89"|
                    Age_group == "90 to 94"| Age_group == "95 plus"),
                    burden_psa * Efficacy * Coverage,                    
                                                                                                                                    
                     0)))))))))))
    
  } else {
  vaccine_age <- vaccine_age %>%
    mutate(va_age =
            ifelse((VC=="4 week (effective at 6 week)" | 
                   VC == "4 week (effective at 6 week), 60 year") &
                   Duration=="2 years" & Age_group == "PN", 
                   burden_psa * Efficacy * Coverage * 47/48,
            ifelse((VC=="4 week (effective at 6 week)" | 
                   VC == "4 week (effective at 6 week), 60 year") &
                   Duration=="2 years" & Age_group == "1 to 4", 
                burden_psa * Efficacy * Coverage * (1/4 + 1/4*5/52),
             
             
            ifelse((VC=="4 week (effective at 6 week)" | 
                   VC=="4 week (effective at 6 week), 60 year") & 
                   Duration=="5 years" & Age_group == "PN", 
                   burden_psa * Efficacy * Coverage * 47/48,
            ifelse((VC=="4 week (effective at 6 week)"| 
                    VC=="4 week (effective at 6 week), 60 year") & 
                    Duration=="5 years" & Age_group == "1 to 4", 
                   burden_psa * Efficacy * Coverage,
            ifelse((VC=="4 week (effective at 6 week)" | 
                   VC=="4 week (effective at 6 week), 60 year") & 
                   Duration=="5 years" & Age_group == "5 to 9", 
                   burden_psa * Efficacy * Coverage * 1/5 *5/52,
                                        
                   
            ifelse(VC=="10 week (effective at 12 week)" & Duration=="2 years" 
                   & Age_group == "PN", 
                   burden_psa * Efficacy * Coverage * 41/48,
            ifelse(VC=="10 week (effective at 12 week)" & Duration=="2 years" 
                   & Age_group == "1 to 4", 
                   burden_psa * Efficacy * Coverage * (1/4 + 1/4 * 11/52),
                                                    
                                                                                    
            ifelse(VC=="10 week (effective at 12 week)" & Duration=="5 years"
                   & Age_group == "PN",
                   burden_psa * Efficacy * Coverage * 41/48,
            ifelse(VC=="10 week (effective at 12 week)" & Duration=="5 years" 
                   & Age_group == "1 to 4", 
                   burden_psa * Efficacy * Coverage,
            ifelse(VC=="10 week (effective at 12 week)" & Duration=="5 years" 
                   & Age_group == "5 to 9", 
                   burden_psa * Efficacy * Coverage * 1/5 * 11/52,
                   
                   
            ifelse(VC=="9 month" & Duration=="20 years"
                   & Age_group == "PN",
                   burden_psa * Efficacy * Coverage * 4/11,
            ifelse(VC=="9 month" & Duration=="20 years" 
                   & (Age_group == "1 to 4"   | Age_group == "5 to 9"|
                      Age_group == "10 to 14" | Age_group == "15 to 19"), 
                   burden_psa * Efficacy * Coverage,
            ifelse(VC=="9 month" & Duration=="20 years" 
                   & Age_group == "20 to 24", 
                   burden_psa * Efficacy * Coverage * 1/5 * 8/12,
                                                             
                                                                                                                                                                    
            ifelse(VC=="4 week (effective at 6 week), 60 year" & 
                   Duration=="2 years" & Age_group == "60 to 64", 
                   burden_psa * Efficacy * Coverage * 2/5,                    
                                                                                                                                                                                           
            ifelse(VC=="60 year" & Duration=="1 year" & Age_group == "60 to 64", 
                   burden_psa * Efficacy * Coverage * 1/5,
                   
            ifelse((VC=="60 year" | VC== "4 week (effective at 6 week), 60 year") & 
                   Duration=="5 years" & Age_group == "60 to 64", 
                   burden_psa * Efficacy * Coverage,
                   
            ifelse(VC=="60 year" & Duration=="10 years" & 
                   (Age_group == "60 to 64"| Age_group == "65 to 69"), 
                  burden_psa * Efficacy * Coverage,
                                 
                   0))))))))))))))))))
  }
  
  # applying vaccine target disease presentation
  vaccine_impact <- data.table(vaccine_age)
  
  vaccine_impact[, va_health_burden := 0]
  
  vaccine_impact[DP  == "All" |
                   
                   (DP  == "BSI" & Disease_presentation == "BSI") |
                   
                   (DP  == "BSI, LRI and thorax infections" & 
                      (Disease_presentation == "BSI" | 
                       Disease_presentation == "LRI and thorax infections")) |
                   
                   (DP  == "BSI, CNS infections, Cardiac infections, LRI" & 
                      (Disease_presentation == "BSI" | 
                       Disease_presentation == "CNS infections" | 
                       Disease_presentation == "Cardiac infections" |
                       Disease_presentation == "LRI and thorax infections")) |
                   
                   (DP  == "Diarrhoea" & Disease_presentation == "Diarrhoea") |  
                   
                   (DP  == "LRI" & Disease_presentation == "LRI and thorax infections") |
                   
                   (DP  == "UTI" & Disease_presentation == "UTI"),
                 
                 va_health_burden := va_age]
  
  return(vaccine_impact)

 } # end of function -- estimate_vaccine_impact

# ------------------------------------------------------------------------------

  
    
#-------------------------------------------------------------------------------
# [table in appendix] vaccine avertable health burdens associated with and attributable
# to AMR by WHO region, pathogen, disease presentation, and age group
  
  update_death_burden <- function(combined_dt,
                                  death_burden_dt,
                                  AMR_burden_data_updated_file){
    
  associated_mean     <- estimate_vaccine_impact(i = "ANY", data = combined_dt, 
                                                 burden_psa = combined_dt$Associated_resistant_mean)[,va_health_burden]
    
    
  associated_lower    <- estimate_vaccine_impact(i = "ANY", data = combined_dt, 
                                                 burden_psa = combined_dt$Associated_resistant_lower)[,va_health_burden]
    
  associated_upper    <- estimate_vaccine_impact(i = "ANY", data = combined_dt,
                                                 burden_psa = combined_dt$Associated_resistant_upper)[,va_health_burden]
    
    
  attributable_mean   <- estimate_vaccine_impact(i = "ANY", data = combined_dt, 
                                                 burden_psa = combined_dt$Attributable_resistance_mean)[,va_health_burden]
    
    
  attributable_lower  <- estimate_vaccine_impact(i = "ANY", data = combined_dt, 
                                                 burden_psa = combined_dt$Attributable_resistance_lower)[,va_health_burden]
    
    
  attributable_upper  <- estimate_vaccine_impact(i = "ANY", data = combined_dt, 
                                                 burden_psa = combined_dt$Attributable_resistance_upper)[,va_health_burden]
    
  x <- data.frame("vaccine_avertable_deaths_associated_with_resistance_(mean)"   =  associated_mean,
                  "vaccine_avertable_deaths_associated_with_resistance_(lower)"  =  associated_lower,
                  "vaccine_avertable_deaths_associated_with_resistance_(upper)"  =  associated_upper,
                  "vaccine_avertable_deaths_attributable_to_resistance_(mean)"   =  attributable_mean,
                  "vaccine_avertable_deaths_attributable_to_resistance_(lower)"  =  attributable_lower,
                  "vaccine_avertable_deaths_attributable_to_resistance_(upper)"  =  attributable_upper)
    
  AMR_burden_data_updated <- cbind(death_burden_dt, x)
    
  AMR_burden_data_updated <- AMR_burden_data_updated %>%
    rename("deaths associated with resistance (mean)"   =  "Associated_resistant_mean",
           "deaths associated with resistance (lower)"  =  "Associated_resistant_lower",
           "deaths associated with resistance (upper)"  =  "Associated_resistant_upper",
           "deaths attributable to resistance (mean)"   =  "Attributable_resistance_mean",
           "deaths attributable to resistance (lower)"  =  "Attributable_resistance_lower",
           "deaths attributable to resistance (upper)"  =  "Attributable_resistance_upper")
    
    # save as xlsx
  fwrite (x    = AMR_burden_data_updated,
          file = AMR_burden_data_updated_file)
    
 return(AMR_burden_data_updated)}
#-------------------------------------------------------------------------------
# [table 2] Deaths and DALYs associated with and attributable to bacterial antimicrobial resistance
# globally and by WHO_region, 2019

  aggregate_impact_by_region <- function(input_data, input_scenario = "conservative"){
    
    impact_by_region <- data.table(WHO_region     = character(),
                                   averted_burden = numeric(),
                                   run_id         = numeric())
    
    for(i in 1:run){
      dt <- estimate_vaccine_impact(i, data= input_data, scenario = input_scenario)
      dt <- dt %>%
        group_by(WHO_region, run_id) %>%
        summarise(averted_burden=sum(va_health_burden), .groups = 'drop')
      impact_by_region <- rbindlist (list (impact_by_region, dt),
                                     use.names = TRUE)} 
    
# -------------------------------------------------------------------------
    WHOregion <- unique(death_burden_dt$WHO_region)
    
    burden_averted_regional <- data.table("50%"=numeric(), "2.5%"=numeric(), "97.5%"=numeric())
    
    impact_by_region_dt <- impact_by_region
    
    for(i in WHOregion){
      dt <- impact_by_region_dt %>% filter(WHO_region == i)
      dt <- quantile(x = dt$averted_burden, probs = c (0.5, 0.025, 0.975))
      dt <- data.table(t(dt))
      burden_averted_regional <- rbindlist (list (burden_averted_regional, dt),
                                            use.names = FALSE)
    }
    
    burden_averted_global <- impact_by_region_dt %>%
      group_by(run_id) %>%
      summarise(averted_burden=sum(averted_burden), .groups = 'drop')
    
    burden_averted_global <- quantile(x = burden_averted_global$averted_burden, probs = c (0.5, 0.025, 0.975))
    burden_averted_global <- data.table(t(burden_averted_global))
    
    Burden_Averted <- rbind(burden_averted_global, burden_averted_regional)
    
    Burden_Averted <- cbind(data.table(Counts=c("Global", WHOregion)), Burden_Averted)
    
    Burden_Averted <- Burden_Averted[,c("Counts", "2.5%", "50%", "97.5%")]
    
    return(Burden_Averted)
    
  } # end of function -- aggregate_impact_by_region

# -------------------------------------------------------------------------
# table 2 summary
  estimate_bruden_averted_by_region <- function(input_data){
    
    burden_averted <- data.table(input_data)
    
    burden_averted [,2:7] <- lapply(burden_averted[,2:7], function(x) comma(x,  format = "d"))
    
    burden_averted [,"Associated with resistance" := 
                      paste(burden_averted$"50%.x","(",burden_averted$"2.5%.x","-",
                            burden_averted$"97.5%.x",")")]
    
    burden_averted [,"Attributable to resistance" := 
                      paste(burden_averted$"50%.y","(",burden_averted$"2.5%.y","-",
                            burden_averted$"97.5%.y",")")]
    
    burden_averted <- burden_averted [, c("Counts", "Associated with resistance", 
                                          "Attributable to resistance")]
    
    burden_averted <- burden_averted [c(2,3,4,5,6,7,1),]
    
    return(burden_averted)
  }
  
# ------------------------------------------------------------------------------
# Figure 1 -- Create vaccine averted burden by region

create_burden_averted_by_region_graph  <- function(Attributable_burden_averted,
                                                   Associated_burden_averted,
                                                   ylim_max,
                                                   ylabel){

Associated_burden_averted$Resistance   <- "Associated with resistance" 

Attributable_burden_averted$Resistance <- "Attributable to resistance" 
  
burden_averted_by_region <- rbind(Associated_burden_averted, Attributable_burden_averted)

burden_averted_by_region <- burden_averted_by_region %>% rename("lower_value"  = "2.5%",
                                                                "median_value" = "50%",
                                                                "upper_value"  = "97.5%")

burden_averted_by_region  <- burden_averted_by_region %>% filter(Counts != "Global")
                   
ggplot(burden_averted_by_region, aes(x = reorder(Counts, -median_value), y=median_value, fill=Resistance)) +
  geom_bar(stat = "identity", position="dodge") +
  scale_fill_manual(values = c("#D4E3FF","#054C70")) +
  labs(x = "WHO region", y = paste(ylabel)) +
  ylim(0, ylim_max) +
  geom_errorbar(aes(ymin=lower_value, ymax=upper_value), width=0.25,
                size=0.5, position=position_dodge(0.9)) +
  theme_classic() +
  theme(legend.position = c(0.8, 0.9))

} # end of function -- create_burden_averted_by_region_graph

# ------------------------------------------------------------------------------
# Figure 2 Create vaccine averted burden by disease presentation

aggregate_impact_by_dp <- function(data_input, input_scenario = "conservative"){
  
  # create table for avertable burden to AMR by disease presentation
  impact_by_dp <- data.table(Disease_presentation = character(), 
                             averted_burden       = numeric(), 
                             run_id               = numeric())
  
  for(i in 1:run){
    dt <- estimate_vaccine_impact(i, data = data_input, scenario = input_scenario)
    dt <- dt %>%
      group_by(Disease_presentation, run_id) %>%
      summarise(averted_burden = sum(va_health_burden), .groups = 'drop')
    impact_by_dp <- rbindlist (list (impact_by_dp, dt),
                               use.names = TRUE) 
  }
  
  # -------------------------------------------------------------------------
  
  DiseasePresentation <- unique(death_burden_dt$Disease_presentation)
  
  death_averted_dp    <- data.table("50%"=numeric(), "2.5%"=numeric(), "97.5%"=numeric())
  
  impact_by_dp_dt     <- impact_by_dp
  
  for(i in DiseasePresentation){
    dt <- impact_by_dp_dt %>% filter(Disease_presentation == i)
    dt <- quantile(x = dt$averted_burden, probs = c (0.5, 0.025, 0.975))
    dt <- data.table(t(dt))
    death_averted_dp <- rbindlist (list (death_averted_dp, dt),
                                   use.names = FALSE)
  }
  
  Death_Aveted <- cbind(data.table(Counts=DiseasePresentation, death_averted_dp))
  
  Death_Aveted <- Death_Aveted[, c("Counts", "2.5%", "50%", "97.5%")]
  
  return(Death_Aveted)} # end of function -- aggregate_impact_by_dp


# create Figure 2

create_burden_averted_by_dp_graph <- function(Attributable_burden_averted,
                                              Associated_burden_averted,
                                              ylim_max,
                                              ylabel){
  
  Associated_burden_averted$Resistance   <- "Associated with resistance" 
  
  Attributable_burden_averted$Resistance <- "Attributable to resistance" 
  
  
    burden_averted_by_dp <- rbind(Associated_burden_averted, Attributable_burden_averted)
  
    burden_averted_by_dp<-  burden_averted_by_dp %>% rename("lower_value"  = "2.5%",
                                                            "median_value" = "50%",
                                                            "upper_value"  = "97.5%")

    burden_averted_by_dp$Counts <- gsub("Bone and joint infections", "Bone+",
                                        burden_averted_by_dp$Counts)
    
    burden_averted_by_dp$Counts <- gsub("Cardiac infections", "Cardiac",
                                        burden_averted_by_dp$Counts)
    
    burden_averted_by_dp$Counts <- gsub("CNS infections", "CNS",
                                        burden_averted_by_dp$Counts)
    
    burden_averted_by_dp$Counts <- gsub("Intra-abdominal infections", "Intra-abdominal",
                                        burden_averted_by_dp$Counts)
    
    burden_averted_by_dp$Counts <- gsub("TB", "Tuberculosis",
                                        burden_averted_by_dp$Counts)
    
    burden_averted_by_dp$Counts <- gsub("LRI and thorax infections", "LRI+",
                                        burden_averted_by_dp$Counts)
    
    burden_averted_by_dp$Counts <- gsub("Bacterial skin infections", "Skin",
                                        burden_averted_by_dp$Counts)
    
    burden_averted_by_dp$Counts <- gsub("Typhoid, paratyphoid, and iNTS", "TF-PF-iNTS",
                                        burden_averted_by_dp$Counts)
      
    burden_averted_by_dp_graph <- ggplot(burden_averted_by_dp, 
                                       aes(x = reorder(Counts, -median_value), 
                                           y = median_value, fill = Resistance)) +

  geom_bar(stat = "identity", position="dodge") +
    scale_fill_manual(values = c("#D4E3FF","#054C70")) +
    labs(x = "Infectious syndrome", y = paste(ylabel)) + 
    ylim(0,ylim_max) +
    geom_errorbar(aes(ymin=lower_value, ymax=upper_value), width=0.25,
                  size=0.5, position=position_dodge(0.9)) +
    theme_classic() +
    theme(legend.position = c(0.9, 0.9))
  
  return(burden_averted_by_dp_graph)
  } # end of function -- create_burden_averted_by_dp_graph

# ------------------------------------------------------------------------------
# Figure 3 Create vaccine averted burden by pathogen

aggregate_impact_by_pathogen <- function(data_input, input_scenario = "conservative"){
  
  impact_by_pathogen <- data.table(Pathogen             = character(), 
                                   averted_burden       = numeric(), 
                                   run_id               = numeric())
  
  for(i in 1:run){
    dt <- estimate_vaccine_impact(i, data = data_input, scenario = input_scenario)
    dt <- dt %>%
      group_by(Pathogen, run_id) %>%
      summarise(averted_burden = sum(va_health_burden), .groups = 'drop')
    impact_by_pathogen <- rbindlist (list (impact_by_pathogen, dt),
                                     use.names = TRUE)
  }
  # -------------------------------------------------------------------------
  burden_averted_pathogen <- data.table("50%"=numeric(), "2.5%"=numeric(), "97.5%"=numeric())
  
  impact_by_pathogen_dt <- impact_by_pathogen
  
  for(i in pathogenlist){
    dt <- impact_by_pathogen_dt %>% filter(Pathogen == i)
    dt <- quantile(x = dt$averted_burden, probs = c (0.5, 0.025, 0.975))
    dt <- data.table(t(dt))
    burden_averted_pathogen <- rbindlist (list (burden_averted_pathogen, dt),
                                          use.names = FALSE)
  }
  
  Burden_Averted <- cbind(data.table(Counts=pathogenlist, burden_averted_pathogen))
  
  Burden_Averted <- Burden_Averted[,c("Counts", "2.5%", "50%", "97.5%")]
  
  return(Burden_Averted)} # end of function -- aggregate_impact_by_pathogen

# create Figure 3

create_burden_averted_by_pathogen_graph <- function(Attributable_burden_averted,
                                                    Associated_burden_averted,
                                                    ylim_max,
                                                    ylabel){
  
  Associated_burden_averted$Resistance    <- "Associated with resistance" 
  
  Attributable_burden_averted$Resistance  <- "Attributable to resistance" 
  
  
  burden_averted_by_pathogen <- rbind(Associated_burden_averted, Attributable_burden_averted)
  
  burden_averted_by_pathogen<-  burden_averted_by_pathogen %>% rename("lower_value" = "2.5%",
                                                                      "median_value" = "50%",
                                                                      "upper_value" = "97.5%")
  
  burden_averted_by_pathogen_graph  <- ggplot(burden_averted_by_pathogen, 
                                              aes(x = reorder(Counts, -median_value), 
                                                  y=median_value, fill=Resistance)) +
    geom_bar(stat = "identity", position="dodge") +
    scale_fill_manual(values = c("#D4E3FF","#054C70")) +
    labs(x = "Pathogen", y = paste(ylabel)) + 
    ylim(0,ylim_max) +
    geom_errorbar(aes(ymin = lower_value, ymax = upper_value), width=0.25,
                  size=0.5, position=position_dodge(0.9)) +
    theme_classic(base_size=8) +
    theme(legend.position = c(0.9, 0.9))
  
  return(burden_averted_by_pathogen_graph)
} # end of function -- create_burden_averted_by_pathogen_graph
# ------------------------------------------------------------------------------

