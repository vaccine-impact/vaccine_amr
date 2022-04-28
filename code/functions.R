# ------------------------------------------------------------------------------
# functions.R
#
# functions for the main analysis
# to estimate the vaccine avertable bacterial antimicrobial resistance burden 
# based on profiles of existing and future vaccines 
# at the global and regional levels
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# estimate vaccine coverage of existing vaccine -- HIB vaccine, PCV

existing_vaccine_coverage <- function(year,
                                      hib_coverage_file,
                                      pcv_coverage_file){
  
  # import WPP population data
  population_file <- read_excel("data/WPP_population_by_age.xlsx",
                                col_names = FALSE)
  
  names(population_file) <- lapply(population_file[13, ], as.character)
  population_file <- population_file %>% filter(Type=="Country/Area")
  
  # extracting data of the year
  population_file <- population_file[population_file$`Reference date (as of 1 July)` == year, 
                                     c("Region, subregion, country or area *", 
                                       "Reference date (as of 1 July)", "0")]
  
 # ----------------------------------------------------------------------------
  # import ISO3 information data
  WPP_iso3 <- read_excel("data/WPP_locations.XLSX", 
                         col_names = FALSE)
  
  names(WPP_iso3) <- lapply(WPP_iso3[13, ], as.character)
  
  WPP_iso3 <- WPP_iso3[-c(1:13),c(2,5)]
  
  # ----------------------------------------------------------------------------
  # adding ISO3 to WPP2019_POPULATION data frame
  population_iso3 <- left_join(population_file, WPP_iso3, 
                               by=c("Region, subregion, country or area *"
                                   = "Region, subregion, country or area*"))
  
  # extract the population at age of 0 and transform it to the actual value
  population_iso3$`0`<- as.numeric(population_iso3$`0`) * 1000
  
  population_iso3 <- rename(population_iso3, "Population" = "0")
  
  # ----------------------------------------------------------------------------
  # import WUENIC vaccine coverage data for HIB vaccine and PCV
  WUENIC <- read_excel("data/WUENIC.xlsx", sheet = "WUENIC_input_to_PDF")
  WUENIC <- WUENIC[, c("ISOCountryCode", "Year", "Vaccine", "WUENIC")]
  
  WUENIC$"ISOCountryCode" = toupper(WUENIC$"ISOCountryCode")
  
  WUENIC_hib <- WUENIC %>% filter(Year == year & Vaccine == "hib3")
  WUENIC_pcv <- WUENIC %>% filter(Year == year & Vaccine == "pcv3")
  
  # ----------------------------------------------------------------------------
  # import WHO region classification data
  country_region <- read_csv("data/country_income_region_classification.csv")
  country_region <- country_region[, c("iso3_code", "WHO_region")]
  
  # ----------------------------------------------------------------------------
  # combine population data with current vaccine coverage data
  WUENIC_WPP_hib <- left_join(WUENIC_hib, population_iso3, 
                              by=c("ISOCountryCode" = "ISO3 Alpha-code"))
  WUENIC_WPP_pcv <- left_join(WUENIC_pcv, population_iso3, 
                              by=c("ISOCountryCode" = "ISO3 Alpha-code"))
  
  # add region
  WUENIC_WPP_hib <- left_join(WUENIC_WPP_hib, country_region, 
                              by=c("ISOCountryCode" = "iso3_code"))
  WUENIC_WPP_pcv <- left_join(WUENIC_WPP_pcv, country_region, 
                              by=c("ISOCountryCode" = "iso3_code"))
  
  # remove NA
  # -- population data is missing for some countries
  # -- Palestine is not a member state of WHO
  WUENIC_WPP_hib <- na.omit(WUENIC_WPP_hib)
  WUENIC_WPP_pcv <- na.omit(WUENIC_WPP_pcv)
  
  # ----------------------------------------------------------------------------
  # estimate regional coverage -- HIB vaccine
  WUENIC_WPP_hib$vaccinated_pop <- WUENIC_WPP_hib$WUENIC * 1/100 * WUENIC_WPP_hib$Population 
  
  hib_vaccinated <- WUENIC_WPP_hib %>% 
    group_by(WHO_region) %>%
    summarise(vaccinated_pop = sum(vaccinated_pop))
  
  hib_total_pop <- WUENIC_WPP_hib %>% 
    group_by(WHO_region) %>%
    summarise(total_pop = sum(Population))
  
  hib_vaccinated$total_pop <- hib_total_pop$total_pop
  
  hib_vaccinated$vaccine_coverage <- hib_vaccinated$vaccinated_pop / hib_vaccinated$total_pop
  
  # estimate regional coverage -- PCV
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
# create and clean up the AMR burden data

create_burden_table <- function(AMR_burden,
                                burden_file){
  
  AMR_burden <- AMR_burden[-c(1:2),]
  
  names(AMR_burden) <- c("WHO_region", "Disease_presentation","Age_group","Pathogen",
                         "Susceptible_mean", "Susceptible_lower", "Susceptible_upper",
                         "Associated_resistant_mean", "Associated_resistant_lower",
                         "Associated_resistant_upper", "Attributable_resistance_mean",
                         "Attributable_resistance_lower", "Attributable_resistance_upper")
  
  AMR_burden[,5:13] <- lapply(AMR_burden[,5:13], function(x) gsub(",", "", x))
  
  AMR_burden[,5:13] <- lapply(AMR_burden[,5:13], as.numeric)
  
  AMR_burden$Age_group <- factor(AMR_burden$Age_group, 
                                 levels=unique(AMR_burden$Age_group), order=T)
  
  levels(AMR_burden$Age_group) <- c("EN", "LN", "PN",  "1 to 4", "5 to 9","10 to 14",
                                    "15 to 19", "20 to 24", "25 to 29", "30 to 34",
                                    "35 to 39", "40 to 44", "45 to 49", "50 to 54" ,      
                                    "55 to 59", "60 to 64", "65 to 69", "70 to 74",
                                    "75 to 79", "80 to 84", "85 to 89", "90 to 94", 
                                    "95 plus")
  
  AMR_burden <- AMR_burden %>%
    filter(AMR_burden$WHO_region != "unclassified")
  
  # save the file
  fwrite (x    = AMR_burden,
          file = burden_file)
  
  return(AMR_burden)
} # end of function -- create_burden_table
# ------------------------------------------------------------------------------  



# ------------------------------------------------------------------------------
# create and clean up the data frame of WHO vaccine profile

create_vaccine_profile_table <- function(vaccine_profile,
                                         vaccine_profile_file){
  
  vaccine_profile <- vaccine_profile[, c("Vaccine - Pathogen", 
                                         "Pathogen", 
                                         "Efficacy (%)",
                                         "Coverage in target group",
                                         "Duration of protection",
                                         "Disease presentation",
                                         "Age of vaccination (Baseline)",
                                         "Age of vaccination (High-potential)",
                                         "Main analysis",
                                         "Elderly age group with highest burden")]
  
  vaccine_profile <- vaccine_profile %>%
    rename("Vaccine_pathogen" = "Vaccine - Pathogen",
           "Efficacy"     = "Efficacy (%)",
           "Coverage"     = "Coverage in target group",
           "Duration"     = "Duration of protection",
           "DP"           = "Disease presentation",
           "VC"           = "Age of vaccination (Baseline)",
           "VO"           = "Age of vaccination (High-potential)",
           "MainAnalysis" = "Main analysis",
           "HBG"          = "Elderly age group with highest burden")
  
  vaccine_profile$Efficacy <- vaccine_profile$Efficacy * 1/100
  vaccine_profile$Coverage <- vaccine_profile$Coverage * 1/100
  
  fwrite (x    = vaccine_profile, 
          file = vaccine_profile_file)
  
  return(vaccine_profile)
} # end of function -- create_vaccine_profile_table

# ------------------------------------------------------------------------------
# create combined table: disease burden + vaccine profile

create_combined_table <- function(death_burden_dt,
                                  vaccine_profile_dt,
                                  susceptible_burden_file,
                                  attributable_burden_file,
                                  associated_burden_file){
  
  # filter vaccines for inclusion in the main anlaysis
  vaccine_profile <- vaccine_profile_dt %>% 
    filter(MainAnalysis == "Yes")
  
  vaccine_profile <- data.table(vaccine_profile[-c(3:4,8), -1])

  vaccine_profile[Pathogen == "Escherichia coli", DP := "BSI, Diarrhoea, UTI"]

  # create combined table
  combined_table <- data.table(left_join(death_burden_dt, vaccine_profile,
                                         by=c("Pathogen" = "Pathogen")))

  # Streptococcus pneumoniae vaccine specification
  combined_table[Pathogen             == "Streptococcus pneumoniae" &
                 Disease_presentation == "LRI and thorax infections",
                 Efficacy             := 0.25]
  
  # Escherichia coli vaccine specification
  combined_table[Pathogen             == "Escherichia coli" &
                 Disease_presentation == "Diarrhoea",
                 Efficacy             := 0.6]
  
  combined_table[Pathogen             == "Escherichia coli" &
                 Disease_presentation == "Diarrhoea",
                 VC                   := "6 months"]
  
  combined_table[Pathogen             == "Escherichia coli" &
                 Disease_presentation == "Diarrhoea",
                 VO                   := "6 months"]

  
  # Using maternal and infant vaccines
  combined_table[Pathogen == "Klebsiella pneumoniae" & Disease_presentation == "BSI" &
                   (Age_group=="EN" | Age_group=="LN" | Age_group=="PN"), 
                 VC := "0 weeks ~"]
  
  # separate burden type
  susceptible_burden <- combined_table[, -c("Associated_resistant_mean",
                                            "Associated_resistant_lower",
                                            "Associated_resistant_upper",
                                            "Attributable_resistance_lower",
                                            "Attributable_resistance_mean",
                                            "Attributable_resistance_upper")]
  
  susceptible_burden  <- susceptible_burden  %>% 
    rename("burden_lower_value" = "Susceptible_lower",
           "burden_mean_value"  = "Susceptible_mean",
           "burden_upper_value" = "Susceptible_upper")
  
  associated_burden <- combined_table[, -c("Susceptible_lower",
                                           "Susceptible_mean",
                                           "Susceptible_upper",
                                           "Attributable_resistance_mean",
                                           "Attributable_resistance_lower",
                                           "Attributable_resistance_upper")]
  
  associated_burden <- associated_burden %>%
    rename("burden_lower_value" = "Associated_resistant_lower",
           "burden_mean_value"  = "Associated_resistant_mean",
           "burden_upper_value" = "Associated_resistant_upper")
  
  attributable_burden <- combined_table[, -c("Susceptible_lower",
                                             "Susceptible_mean",
                                             "Susceptible_upper",
                                             "Associated_resistant_mean",
                                             "Associated_resistant_lower",
                                             "Associated_resistant_upper")]
  
  attributable_burden <- attributable_burden %>% 
    rename("burden_lower_value" = "Attributable_resistance_lower",
           "burden_mean_value"  = "Attributable_resistance_mean",
           "burden_upper_value" = "Attributable_resistance_upper")
  
  fwrite (x    = susceptible_burden,
          file = susceptible_burden_file)
  
  fwrite (x    = associated_burden,
          file = associated_burden_file)
  
  fwrite (x    = attributable_burden,
          file = attributable_burden_file)
  
  return(combined_table)
  
} # end of function -- create_combined_table
# ------------------------------------------------------------------------------
# estimate pre-vaccine burden for existing vaccines -- HIB vaccine & PCV

estimate_prevaccination_burden <- function(burden_input,
                                           burden_file){
  
  # import current vaccine coverage estimates by region  
  hib_coverage_2019 <- read_csv("tables/hib coverage 2019.csv")
  hib_coverage_2018 <- read_csv("tables/hib coverage 2018.csv")
  
  pcv_coverage_2019 <- read_csv("tables/pcv coverage 2019.csv")
  pcv_coverage_2018 <- read_csv("tables/pcv coverage 2018.csv")
  
  hib_coverage <- data.table (WHO_region = c("Africa", "Americas", "Eastern Mediterranean", 
                                             "Europe", "South-East Asia", "Western Pacific"), 
                              hib_vaccine_coverage_2019 = hib_coverage_2019$vaccine_coverage,
                              hib_vaccine_coverage_2018 = hib_coverage_2018$vaccine_coverage)
  
  pcv_coverage <- data.table (WHO_region = c("Africa", "Americas", "Eastern Mediterranean", 
                                             "Europe", "South-East Asia", "Western Pacific"), 
                              pcv_vaccine_coverage_2019 = pcv_coverage_2019$vaccine_coverage,
                              pcv_vaccine_coverage_2018 = pcv_coverage_2018$vaccine_coverage)
  
  # substitute the current burden estimates to pre-vaccine burden estimates for existing vaccines
  AMR_burden <- burden_input
  
  AMR_burden <- left_join(AMR_burden, hib_coverage, by=c("WHO_region" = "WHO_region"))
  
  AMR_burden <- left_join(AMR_burden, pcv_coverage, by=c("WHO_region" = "WHO_region"))
  
  AMR_burden <- data.table(AMR_burden)
  
  # estimate HIB pre-vaccination burden
  AMR_burden[Pathogen  == "Haemophilus influenzae" & Age_group == "PN",
             burden_mean_value := burden_mean_value /
               (3/48 + 4/48 * (1 - hib_vaccine_coverage_2019 * 0.59 * 0.95) +
                  4/48 * (1 - hib_vaccine_coverage_2019 * 0.92 * 0.95) + 
                  37/48 * (1 - hib_vaccine_coverage_2019 * 0.93 * 0.95))]
  
  AMR_burden[Pathogen  == "Haemophilus influenzae" & Age_group == "PN",
             burden_lower_value := burden_lower_value /
               (3/48 + 4/48 * (1 - hib_vaccine_coverage_2019 * 0.59 * 0.95) +
                  4/48 * (1 - hib_vaccine_coverage_2019 * 0.92 * 0.95) + 
                  37/48 * (1 - hib_vaccine_coverage_2019 * 0.93 * 0.95))]
  
  AMR_burden[Pathogen  == "Haemophilus influenzae" & Age_group == "PN",
             burden_upper_value := burden_upper_value /
               (3/48 + 4/48 * (1 - hib_vaccine_coverage_2019 * 0.59 * 0.95) + 
                  4/48 * (1 - hib_vaccine_coverage_2019 * 0.92 * 0.95) + 
                  37/48 * (1 - hib_vaccine_coverage_2019 * 0.93 * 0.95))]
  
  AMR_burden[Pathogen  == "Haemophilus influenzae" & Age_group == "1 to 4",
             burden_mean_value := burden_mean_value /
               (1 - hib_vaccine_coverage_2018 * 0.93 * 0.95)]
  
  AMR_burden[Pathogen  == "Haemophilus influenzae" & Age_group == "1 to 4",
             burden_lower_value := burden_lower_value /
               (1 - hib_vaccine_coverage_2018 * 0.93 * 0.95)]
  
  AMR_burden[Pathogen  == "Haemophilus influenzae" & Age_group == "1 to 4",
             burden_upper_value := burden_upper_value /
               (1 - hib_vaccine_coverage_2018 * 0.93 * 0.95)]
  
  # estimate streptococcus pneumoniae pre-vaccination burden
  AMR_burden[Pathogen  == "Streptococcus pneumoniae" & Age_group == "PN",
             burden_mean_value := burden_mean_value /
               (3/48 + 4/48 * (1 - hib_vaccine_coverage_2019 * 0.29) + 
                  4/48 * (1 - pcv_vaccine_coverage_2019 * 0.58) + 
                  37/48 * (1 - pcv_vaccine_coverage_2019 * 0.58))]
  
  AMR_burden[Pathogen  == "Streptococcus pneumoniae" & Age_group == "PN",
             burden_lower_value := burden_lower_value /
               (3/48 + 4/48 * (1 - hib_vaccine_coverage_2019 * 0.29) + 
                  4/48 * (1 - pcv_vaccine_coverage_2019 * 0.58) + 
                  37/48 * (1 - pcv_vaccine_coverage_2019 * 0.58))]
  
  AMR_burden[Pathogen  == "Streptococcus pneumoniae" & Age_group == "PN",
             burden_upper_value := burden_upper_value /
               (3/48 + 4/48 * (1 - hib_vaccine_coverage_2019 * 0.29) + 
                  4/48 * (1 - pcv_vaccine_coverage_2019 * 0.58) + 
                  37/48 * (1 - pcv_vaccine_coverage_2019 * 0.58))]
  
  AMR_burden[Pathogen  == "Streptococcus pneumoniae" & Age_group == "1 to 4",
             burden_mean_value := burden_mean_value /
               (1 - pcv_vaccine_coverage_2018 * 0.58)]
  
  AMR_burden[Pathogen  == "Streptococcus pneumoniae" & Age_group == "1 to 4",
             burden_lower_value := burden_lower_value /
               (1 - pcv_vaccine_coverage_2018 * 0.58)]
  
  AMR_burden[Pathogen  == "Streptococcus pneumoniae" & Age_group == "1 to 4",
             burden_upper_value := burden_upper_value /
               (1 - pcv_vaccine_coverage_2018 * 0.58)]
  
  # save the file
  fwrite (x    = AMR_burden,
          file = burden_file)
  
} # end of function -- estimate_prevaccination_burden
# ------------------------------------------------------------------------------



# ------------------------------------------------------------------------------
# create burden trend graph by pathogen across all age groups

create_burden_by_pathogen_graph <- function(pathogen,
                                            input_data,
                                            ylabel,
                                            burden_type){
  burden_data <- input_data[input_data$Pathogen == pathogen, ]
  
  burden_data <- burden_data %>% 
    group_by(Pathogen, WHO_region, Age_group) %>%
    summarise(Associated_resistant_mean=sum(Associated_resistant_mean), .groups = 'drop') 
  
  ggplot(burden_data, aes(x=Age_group)) +
    geom_line(aes(y=Associated_resistant_mean, group=WHO_region, colour=WHO_region)) +
    ylab(paste(ylabel)) + 
    xlab("Age Group") +
    ggtitle(paste(pathogen)) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5))
  
  # save plot
  ggsave (filename = paste(pathogen,"_",burden_type,".png"),
          path = "figures",
          width = 15,
          height = 6,
          dpi = 600)
} # end of function -- create_burden_by_pathogen_graph

# ------------------------------------------------------------------------------



# ------------------------------------------------------------------------------
# conduct uncertainty analysis

 uncertainty_analysis_baseline <- function(psa, 
                                           tolerance,
                                           data,
                                           mode){
  
  burden_dt <- data.table(data)

  # -------------------------------------------------------------------------  
  # minor changes to fit lognormal distribution
    burden_dt <- burden_dt[!(burden_dt$burden_mean_value == "0" & 
                             burden_dt$burden_lower_value == "0" &  
                             burden_dt$burden_upper_value == "0"),]
    
    if ((mode == "death") == TRUE) {
    # add a small value to zero to avoid log(0)
    burden_dt [burden_lower_value == 0, 
               burden_lower_value := burden_lower_value + 1e-6]
    
    burden_dt [burden_mean_value == 0,
               burden_mean_value:= burden_mean_value + 1e-6]
    
    # adjusting value that generate NAs
    burden_dt [WHO_region == "Africa" &
                 Disease_presentation == "Intra-abdominal infections" & 
                 Age_group == "95 plus" & 
                 Pathogen == "Escherichia coli",
               burden_lower_value:= floor(burden_lower_value)]
    
    burden_dt [WHO_region == "Africa" &
                 Disease_presentation == "Intra-abdominal infections" & 
                 Age_group == "95 plus" & 
                 Pathogen == "Escherichia coli",
               burden_mean_value:= floor(burden_mean_value)]
    
    burden_dt [WHO_region == "Africa" &
                 Disease_presentation == "Intra-abdominal infections" & 
                 Age_group == "95 plus" & 
                 Pathogen == "Escherichia coli",
               burden_upper_value:= floor(burden_upper_value)]
    
    burden_dt [WHO_region == "Western Pacific" &
                 Disease_presentation == "Intra-abdominal infections" & 
                 Age_group == "65 to 69" & 
                 Pathogen == "Staphylococcus aureus",
               burden_lower_value:= floor(burden_lower_value)]
    
    burden_dt [WHO_region == "Western Pacific" &
                 Disease_presentation == "Intra-abdominal infections" & 
                 Age_group == "65 to 69" & 
                 Pathogen == "Staphylococcus aureus",
               burden_mean_value:= floor(burden_mean_value)]
    
    burden_dt [WHO_region == "Western Pacific" &
                 Disease_presentation == "Intra-abdominal infections" & 
                 Age_group == "65 to 69" & 
                 Pathogen == "Staphylococcus aureus",
               burden_upper_value:= floor(burden_upper_value)]
    
  } else {
    # add a small value to zero values to avoid log(0)
    burden_dt [burden_lower_value == 0, 
               burden_lower_value := burden_lower_value + 0.05]
    
    burden_dt [burden_mean_value == 0,
               burden_mean_value := burden_mean_value + 0.5]
    
    
    # adjusting value that generate NAs
    burden_dt [burden_lower_value == burden_mean_value, 
               burden_lower_value := burden_lower_value * 0.939]
    
    burden_dt [burden_mean_value == burden_upper_value, 
               burden_mean_value := burden_mean_value * 0.9393]
    
    burden_dt [(WHO_region == "South-East Asia" &
                  Disease_presentation == "BSI" & 
                  Age_group == "75 to 79" & 
                  Pathogen == "Klebsiella pneumoniae")|
                 (WHO_region == "South-East Asia" &
                    Disease_presentation == "Intra-abdominal infections" & 
                    Age_group == "45 to 49" & 
                    Pathogen == "Staphylococcus aureus")|
                 (WHO_region == "South-East Asia" &
                    Disease_presentation == "BSI" & 
                    Age_group == "PN" & 
                    Pathogen == "Non-typhoidal Salmonella") |
                 (WHO_region == "South-East Asia" &
                    Disease_presentation == "BSI" & 
                    Age_group == "40 to 44" & 
                    Pathogen == "Acinetobacter baumannii") |
                 (WHO_region == "Western Pacific" &
                    Disease_presentation == "BSI" & 
                    Age_group == "60 to 64" & 
                    Pathogen == "Staphylococcus aureus") |
                 (WHO_region == "Africa" &
                    Disease_presentation == "LRI and thorax infections" & 
                    Age_group == "1 to 4" & 
                    Pathogen == "Haemophilus influenzae")
               
               , burden_mean_value:= burden_mean_value * 0.94]
    }  
  # ----------------------------------------------------------------------------
  # log-normal distribution
  burden_dt [, burden_sd_log := suppressMessages (suppressWarnings (get.lnorm.par (
                                p           = c (0.025, 0.5, 0.975),
                                q           = c (burden_lower_value, burden_mean_value, burden_upper_value),
                                show.output = FALSE,
                                plot        = FALSE,
                                tol         = tolerance) ["sdlog"]) ),
             
             by = .(WHO_region,Disease_presentation,Age_group,Pathogen)]
  
  burden_dt [, burden_mean_log := suppressMessages (suppressWarnings (get.lnorm.par (
                                p           = c (0.025, 0.5, 0.975),
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
  # rlnorm will generate values based on lognormal distribution
  
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
  estimate_vaccine_impact <- function(data,
                                      scenario = "conservative"){
  
  va <- data.table(data)
  
  if((scenario == "optimistic") == TRUE){
    va[, va_age := 0]
    
    va[VO == "0 weeks (maternal)" & Duration == "6 months" & 
         (Age_group == "EN" | Age_group == "LN"),
       va_age := burden_psa * Efficacy * Coverage]
    va[VO == "0 weeks (maternal)" & Duration == "6 months" & 
         (Age_group=="PN"),
       va_age := burden_psa * Efficacy * Coverage * 5/11]
    
    va[(VO == "6 weeks" | VO == "6 weeks & elderly age group") &
         Duration == "5 years" & Age_group == "PN",
       va_age := burden_psa * Efficacy * Coverage * 45/48]
    va[(VO == "6 weeks" | VO == "6 weeks & elderly age group") &
         Duration == "5 years" & Age_group == "1 to 4",
       va_age := burden_psa * Efficacy * Coverage]
    va[(VO == "6 weeks" | VO == "6 weeks & elderly age group") &
         Duration == "5 years" & Age_group == "5 to 9",
       va_age := burden_psa * Efficacy * Coverage * 1/5 * 7/52]
    
    va[VO == "6 weeks & 9 months" & Duration == "5 years" & Age_group == "PN",
       va_age := burden_psa * Efficacy * Coverage * 45/48]      
    va[VO == "6 weeks & 9 months" & Duration == "5 years" & Age_group == "1 to 4",
       va_age := burden_psa * Efficacy * Coverage]       
    va[VO == "6 weeks & 9 months" & Duration == "5 years" & Age_group == "5 to 9",
       va_age := burden_psa * Efficacy * Coverage * 1/5 * 8/12]
    
    va[VO == "6 months" & Duration == "5 years" & Age_group == "PN",
       va_age := burden_psa * Efficacy * Coverage * 7/11]       
    va[VO == "6 months" & Duration == "5 years" & Age_group == "1 to 4",
       va_age := burden_psa * Efficacy * Coverage]       
    va[VO == "6 months" & Duration == "5 years" & Age_group == "5 to 9",
       va_age := burden_psa * Efficacy * Coverage * 1/5 * 5/12]       
    
    va[VO == "9 months" & Duration == "5 years" & Age_group == "PN",
       va_age := burden_psa * Efficacy * Coverage * 4/11]      
    va[VO == "9 months" & Duration == "5 years" & Age_group == "1 to 4",
       va_age := burden_psa * Efficacy * Coverage]       
    va[VO == "9 months" & Duration == "5 years" & Age_group == "5 to 9",
       va_age := burden_psa * Efficacy * Coverage * 1/5 * 8/12]
    
    va[VO == "9 months" & Duration == "15 years" & Age_group == "PN",
       va_age := burden_psa * Efficacy * Coverage * 4/11] 
    va[VO == "9 months" & Duration == "15 years" &
         (Age_group == "1 to 4"| Age_group == "5 to 9"| Age_group == "10 to 14"),
       va_age := burden_psa * Efficacy * Coverage] 
    va[VO == "9 months" & Duration == "15 years" & Age_group == "15 to 19",
       va_age := burden_psa * Efficacy * Coverage * 1/5 * 8/12] 
    
    va[VO == "15 years" & Duration == "10 years" & 
         (Age_group == "15 to 19"| Age_group == "20 to 24"),
       va_age := burden_psa * Efficacy * Coverage] 
                                                                                                               
    va[VO == "6 weeks & elderly age group" &
         HBG == "65 years" & Duration == "5 years" & 
         Age_group == "65 to 69",
       va_age := burden_psa * Efficacy * Coverage]   
    va[VO == "6 weeks & elderly age group" &
         HBG == "65 years" & Duration == "10 years" & 
         (Age_group == "65 to 69"| Age_group == "70 to 74"),
       va_age := burden_psa * Efficacy * Coverage] 
                
    va[VO == "6 weeks & elderly age group" &
         HBG == "70 years" & Duration == "5 years" & 
         Age_group == "70 to 74",
       va_age := burden_psa * Efficacy * Coverage] 
                                                                                                                                    
    va[VO == "6 weeks & elderly age group" &
         HBG == "75 years" & Duration == "5 years" & 
         Age_group == "75 to 79",
       va_age := burden_psa * Efficacy * Coverage]
    
    va[VO == "10 years + boost every 10 years" &
         (Age_group != "EN"| Age_group != "LN" |
          Age_group != "PN"| Age_group != "1 to 4"|
          Age_group != "5 to 9"),
       va_age := burden_psa * Efficacy * Coverage]
    
    va[(VO=="All age groups" | VO == "0 weeks + boost every 10 years"), 
       va_age := burden_psa * Efficacy * Coverage]

    # using existing Hib vaccine efficacy
    va[Pathogen == "Haemophilus influenzae" & Efficacy == "0.93" &
         VO == "6, 10, 14 weeks" & Age_group == "PN",
       va_age := burden_psa * (4/48 * 0.59 * Coverage +
                               4/48 * 0.92 * Coverage + 
                              37/48 * 0.93 * Coverage) * 0.95]
    
    va[Pathogen == "Haemophilus influenzae" & Efficacy == "0.93" &
         VO == "6, 10, 14 weeks" & Age_group == "1 to 4",
       va_age := burden_psa * Efficacy * Coverage * 0.95]
    
    # using existing PCV efficacy
    va[Pathogen == "Streptococcus pneumoniae" &
         VO == "6, 10, 14 weeks & elderly age group" & Age_group == "PN",
       va_age := burden_psa * (4/48 * 0.29 * Coverage + 
                               4/48 * 0.58 * Coverage+ 
                              37/48 * 0.58 * Coverage)]
    
    va[Pathogen == "Streptococcus pneumoniae" & 
       Disease_presentation == "LRI and thorax infections" &
         VO == "6, 10, 14 weeks & elderly age group" & Age_group == "PN",
       va_age := burden_psa * (4/48 * 0.25 * Coverage + 
                               4/48 * 0.25 * Coverage + 
                              37/48 * 0.25 * Coverage)]
    
    va[Pathogen == "Streptococcus pneumoniae" &
         VO == "6, 10, 14 weeks & elderly age group" &
       (Age_group == "1 to 4" | Age_group == "75 to 79"),
       va_age := burden_psa * Efficacy * Coverage]
    
    
    # apply vaccine assumptions to the ETEC proportion of the disease
    va[Pathogen == "Escherichia coli" & Disease_presentation == "Diarrhoea",
       va_age := va_age * 0.4397]

  } else {
    va[, va_age := 0]
    
    va[VC == "0 weeks (maternal)" & Duration == "6 months" & 
         (Age_group == "EN" | Age_group == "LN"),
       va_age := burden_psa * Efficacy * Coverage]
    va[VC == "0 weeks (maternal)" & Duration == "6 months" & 
         (Age_group == "PN"),
       va_age := burden_psa * Efficacy * Coverage * 5/11]
    
    va[(VC == "6 weeks" | 
          VC == "6 weeks & elderly age group") &
         Duration == "5 years" & Age_group == "PN",
       va_age := burden_psa * Efficacy * Coverage * 45/48]
    va[(VC == "6 weeks" | 
          VC == "6 weeks & elderly age group") &
         Duration == "5 years" & Age_group == "1 to 4",
       va_age := burden_psa * Efficacy * Coverage]
    va[(VC == "6 weeks" | 
          VC == "6 weeks & elderly age group") &
         Duration=="5 years" & Age_group == "5 to 9",
       va_age := burden_psa * Efficacy * Coverage * 1/5 * 7/52]
    
    va[VC == "6 weeks & 9 months" & Duration == "5 years" & Age_group == "PN",
       va_age := burden_psa * Efficacy * Coverage * 45/48]      
    va[VC == "6 weeks & 9 months" & Duration == "5 years" & Age_group == "1 to 4",
       va_age := burden_psa * Efficacy * Coverage]       
    va[VC == "6 weeks & 9 months" & Duration == "5 years" & Age_group == "5 to 9",
       va_age := burden_psa * Efficacy * Coverage * 1/5 * 8/12]
    
    va[VC == "6 months" & Duration == "5 years" & Age_group == "PN",
       va_age := burden_psa * Efficacy * Coverage * 7/11]       
    va[VC == "6 months" & Duration == "5 years" & Age_group == "1 to 4",
       va_age := burden_psa * Efficacy * Coverage]       
    va[VC == "6 months" & Duration == "5 years" & Age_group == "5 to 9",
       va_age := burden_psa * Efficacy * Coverage * 1/5 * 5/12]       
    
    va[VC == "9 months" & Duration == "5 years" & Age_group == "PN",
       va_age := burden_psa * Efficacy * Coverage * 4/11]       
    va[VC == "9 months" & Duration == "5 years" & Age_group == "1 to 4",
       va_age := burden_psa * Efficacy * Coverage]       
    va[VC == "9 months" & Duration == "5 years" & Age_group == "5 to 9",
       va_age := burden_psa * Efficacy * Coverage * 1/5 * 8/12]    
    
    va[VC == "9 months" & Duration == "15 years" & Age_group == "PN",
       va_age := burden_psa * Efficacy * Coverage * 4/11] 
    va[VC == "9 months" & Duration == "15 years" &
         (Age_group == "1 to 4"| Age_group == "5 to 9"| Age_group == "10 to 14"),
       va_age := burden_psa * Efficacy * Coverage] 
    va[VC == "9 months" & Duration == "15 years" & Age_group == "15 to 19",
       va_age := burden_psa * Efficacy * Coverage * 1/5 * 8/12] 
    
    va[VC == "15 years" & Duration == "10 years" & 
         (Age_group == "15 to 19"| Age_group == "20 to 24"),
       va_age := burden_psa * Efficacy * Coverage] 
    
    va[VC == "6 weeks & elderly age group" &
         HBG == "65 years" & Duration == "5 years" & 
         Age_group == "65 to 69",
       va_age := burden_psa * Efficacy * Coverage]   
    va[VC == "6 weeks & elderly age group" &
         HBG == "65 years" & Duration == "10 years" & 
         (Age_group == "65 to 69"| Age_group == "70 to 74"),
       va_age := burden_psa * Efficacy * Coverage] 
    
    va[VC == "6 weeks & elderly age group" &
         HBG == "70 years" & Duration == "5 years" & 
         Age_group == "70 to 74",
       va_age := burden_psa * Efficacy * Coverage] 
    
    va[VC == "6 weeks & elderly age group" &
         HBG == "75 years" & Duration == "5 years" & 
         Age_group == "75 to 79",
       va_age := burden_psa * Efficacy * Coverage] 
    
    va[VC == "10 years + boost every 10 years" &
         (Age_group != "EN"| Age_group != "LN" |
          Age_group != "PN"| Age_group != "1 to 4"|
          Age_group != "5 to 9"), 
       va_age := burden_psa * Efficacy * Coverage]
    
    va[VC == "0 weeks + boost every 10 years", 
       va_age := burden_psa * Efficacy * Coverage]
    
    # using both vaccines for Klebsiella pneumoniae   
    va[Pathogen == "Klebsiella pneumoniae" & 
         Disease_presentation == "BSI" & VC == "0 weeks ~" & 
         (Age_group=="EN" | Age_group=="LN" | Age_group=="PN"),
       va_age := burden_psa * Efficacy * Coverage]
    
    # using existing Hib vaccine efficacy
    va[Pathogen == "Haemophilus influenzae" & Efficacy == "0.93" &
         VC == "6, 10, 14 weeks" & Age_group == "PN",
       va_age := burden_psa * (4/48 * 0.59 * Coverage +
                                4/48 * 0.92 * Coverage + 
                               37/48 * 0.93 * Coverage) * 0.95]
    
    va[Pathogen == "Haemophilus influenzae" & Efficacy == "0.93" &
        VC == "6, 10, 14 weeks" & Age_group == "1 to 4",
       va_age := burden_psa * 0.93 * Coverage * 0.95]
    
    # using existing PCV efficacy
    va[Pathogen == "Streptococcus pneumoniae" &
        VC == "6, 10, 14 weeks" & Age_group == "PN",
      va_age := burden_psa * (4/48 * 0.29 * Coverage + 
                              4/48 * 0.58 * Coverage + 
                             37/48 * 0.58 * Coverage)]
    
    va[Pathogen == "Streptococcus pneumoniae" & 
         Disease_presentation == "LRI and thorax infections" &
         VC == "6, 10, 14 weeks" & Age_group == "PN",
       va_age := burden_psa * (4/48 * 0.25 * Coverage + 
                               4/48 * 0.25 * Coverage + 
                              37/48 * 0.25 * Coverage)]
    
    va[Pathogen == "Streptococcus pneumoniae" &
         VC == "6, 10, 14 weeks" & Age_group == "1 to 4",
       va_age := burden_psa * Efficacy * Coverage]
    
    # apply vaccine assumptions to the ETEC proportion of the disease
    va[Pathogen == "Escherichia coli" & Disease_presentation == "Diarrhoea",
       va_age := va_age * 0.4397]
    
      }
    
  # applying vaccine target disease presentation
  vaccine_impact <- va
  
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

                   (DP  == "BSI, Diarrhoea, UTI" & 
                      (Disease_presentation == "BSI" | 
                       Disease_presentation == "Diarrhoea" | 
                       Disease_presentation == "UTI")) |
                                      
                   (DP  == "Diarrhoea" & Disease_presentation == "Diarrhoea") |  
                   
                   (DP  == "LRI" & Disease_presentation == "LRI and thorax infections") |
                   
                   (DP  == "UTI" & Disease_presentation == "UTI") |
                   
                   (DP  == "gonorrhoea and chlamydia" & 
                      Disease_presentation == "Gonorrhoea and chlamydia"),
                 
                 va_health_burden := va_age]
  
  return(vaccine_impact)

 } # end of function -- estimate_vaccine_impact

# ------------------------------------------------------------------------------

  
    
#-------------------------------------------------------------------------------
# [table in appendix] vaccine avertable health burdens associated with and 
# attributable to AMR by WHO region, pathogen, disease presentation, and age group
  
  update_death_burden <- function(input_associated,
                                  input_attributable,
                                  burden_dt,
                                  AMR_burden_data_updated_file,
                                  input_scenario){
    
  associated_mean <- data.table(input_associated)
  associated_mean[, burden_psa := burden_mean_value]
  associated_mean <- 
    estimate_vaccine_impact(data = associated_mean,
                            scenario = input_scenario)[,va_health_burden]    
    
  associated_lower <- data.table(input_associated)
  associated_lower[, burden_psa := burden_lower_value]
  associated_lower <-
    estimate_vaccine_impact(data = associated_lower,
                            scenario = input_scenario)[,va_health_burden]    
  
  associated_upper <- data.table(input_associated)
  associated_upper[, burden_psa := burden_upper_value]
  associated_upper <-
    estimate_vaccine_impact(data = associated_upper,
                            scenario = input_scenario)[,va_health_burden]    
  
  attributable_mean <- data.table(input_attributable)
  attributable_mean[, burden_psa := burden_mean_value]
  attributable_mean <-
    estimate_vaccine_impact(data = attributable_mean,
                            scenario = input_scenario)[,va_health_burden]    
  
  attributable_lower <- data.table(input_attributable)
  attributable_lower[, burden_psa := burden_lower_value]
  attributable_lower <-
    estimate_vaccine_impact(data = attributable_lower,
                            scenario = input_scenario)[,va_health_burden]    
  
  attributable_upper <- data.table(input_attributable)
  attributable_upper[, burden_psa := burden_upper_value]
  attributable_upper <-
    estimate_vaccine_impact(data = attributable_upper,
                            scenario = input_scenario)[,va_health_burden]    
  
  dt <- data.table("vaccine avertable-associated with resistance (mean)"   =  associated_mean,
                   "vaccine avertable-associated with resistance (lower)"  =  associated_lower,
                   "vaccine avertable-associated with resistance (upper)"  =  associated_upper,
                   "vaccine avertable-attributable to resistance (mean)"   =  attributable_mean,
                   "vaccine avertable-attributable to resistance (lower)"  =  attributable_lower,
                   "vaccine avertable-attributable to resistance (upper)"  =  attributable_upper)
    
  AMR_burden_data_updated <- cbind(burden_dt, dt)
    
  AMR_burden_data_updated <- AMR_burden_data_updated %>%
    rename("associated with resistance (mean)"   =  "Associated_resistant_mean",
           "associated with resistance (lower)"  =  "Associated_resistant_lower",
           "associated with resistance (upper)"  =  "Associated_resistant_upper",
           "attributable to resistance (mean)"   =  "Attributable_resistance_mean",
           "attributable to resistance (lower)"  =  "Attributable_resistance_lower",
           "attributable to resistance (upper)"  =  "Attributable_resistance_upper")
  
    AMR_burden_data_updated <- AMR_burden_data_updated %>% 
    arrange(Pathogen, WHO_region, Disease_presentation, Age_group)
    
    # save as xlsx
  fwrite (x    = AMR_burden_data_updated,
          file = AMR_burden_data_updated_file)
    
 return(AMR_burden_data_updated)}
#-------------------------------------------------------------------------------
  
  
  
#-------------------------------------------------------------------------------
  
# deaths and DALYs associated with and attributable to bacterial antimicrobial resistance
# globally and by WHO_region, 2019

  aggregate_impact_by_region <- function(input_data, 
                                         mode){
   
    if (mode == "incremental") {
      input_data$va_health_burden <- input_data$va_incre
      
    } else if (mode == "high_potential") {
      input_data$va_health_burden <- input_data$va_high
      
    } else {
      input_data$va_health_burden <- input_data$va_base
    }
    
    impact_by_region <- input_data %>%
        group_by(WHO_region, run_id) %>%
        summarise(averted_burden=sum(va_health_burden), .groups = 'drop')
  # -------------------------------------------------------------------------
    WHOregion <- unique(death_burden_dt$WHO_region)
    
    burden_averted_regional <- data.table("50%"=numeric(), "2.5%"=numeric(), "97.5%"=numeric())
    
    for(i in WHOregion){
      dt <- impact_by_region %>% filter(WHO_region == i)
      dt <- quantile(x = dt$averted_burden, probs = c (0.5, 0.025, 0.975))
      dt <- data.table(t(dt))
      burden_averted_regional <- rbindlist (list (burden_averted_regional, dt),
                                            use.names = FALSE)
    }
    
    burden_averted_global <- impact_by_region %>%
      group_by(run_id) %>%
      summarise(averted_burden=sum(averted_burden), .groups = 'drop')
    
    burden_averted_global <- quantile(x = burden_averted_global$averted_burden, 
                                      probs = c (0.5, 0.025, 0.975))
    burden_averted_global <- data.table(t(burden_averted_global))
    
    Burden_Averted <- rbind(burden_averted_regional, burden_averted_global)
    
    Burden_Averted <- cbind(data.table(Counts=c(WHOregion, "Global")), Burden_Averted)
    
    Burden_Averted <- Burden_Averted[,c("Counts", "2.5%", "50%", "97.5%")]
    
    return(Burden_Averted)
    
  } # end of function -- aggregate_impact_by_region

# -------------------------------------------------------------------------
  # create function to edit table
  edit_table <- function(avertable_burden){
    
    avertable_burden [,2:4] <- 
      lapply(avertable_burden[,2:4], function(x) comma(x,  format = "d"))
    
    avertable_burden [, burden_averted := 
                        paste(avertable_burden$"50%","(",
                              avertable_burden$"2.5%","-",
                              avertable_burden$"97.5%",")")]
    
    avertable_burden <- avertable_burden[,c(1,5)]
  }# end of function -- edit_table
  
# -----------------------------------------------------------------------------
  # create vaccine avertable burden table
  
  create_avertable_burden_table <- function(Associated_death_averted,
                                            Attributable_death_averted,
                                            Associated_daly_averted,
                                            Attributable_daly_averted,
                                            Associated_death_averted_opt,
                                            Attributable_death_averted_opt,
                                            Associated_daly_averted_opt,
                                            Attributable_daly_averted_opt,
                                            Associated_death_averted_inc,
                                            Attributable_death_averted_inc, 
                                            Associated_daly_averted_inc,
                                            Attributable_daly_averted_inc){
    
    avertable_death <-
      data.table(Associated_death_averted   = edit_table(Associated_death_averted)[,1:2],
                 Attributable_death_averted = edit_table(Attributable_death_averted)[,2])
    
    avertable_daly <-
      data.table(Associated_daly_averted    = edit_table(Associated_daly_averted)[,1:2],
                 Attributable_daly_averted  = edit_table(Attributable_daly_averted)[,2])
    
    avertable_burden_con <- right_join(avertable_death, avertable_daly, 
                                       by = c("Associated_death_averted.Counts"
                                              = "Associated_daly_averted.Counts"))
    
    avertable_burden_con <- 
      avertable_burden_con %>% 
      rename("Counts"                      = "Associated_death_averted.Counts",
             "Associated_death_baseline"   = "Associated_death_averted.burden_averted",
             "Attributable_death_baseline" = "Attributable_death_averted.burden_averted",
             "Associated_daly_baseline"    = "Associated_daly_averted.burden_averted",
             "Attributable_daly_baseline"  = "Attributable_daly_averted.burden_averted")
    
  # -------------------------------------------------------------------------
    avertable_death_opt <-
      data.table(Associated_death_averted_opt   = edit_table(Associated_death_averted_opt)[,1:2],
                 Attributable_death_averted_opt = edit_table(Attributable_death_averted_opt)[,2])
    
    avertable_daly_opt <-
      data.table(Associated_daly_averted_opt    = edit_table(Associated_daly_averted_opt)[,1:2],
                 Attributable_daly_averted_opt  = edit_table(Attributable_daly_averted_opt)[,2])
    
    avertable_burden_opt <- right_join(avertable_death_opt, avertable_daly_opt, 
                                       by=c("Associated_death_averted_opt.Counts"
                                            = "Associated_daly_averted_opt.Counts"))
    avertable_burden_opt <- 
      avertable_burden_opt %>% 
      rename("Counts"                            = "Associated_death_averted_opt.Counts",
             "Associated_death_high-potential"   = "Associated_death_averted_opt.burden_averted",
             "Attributable_death_high-potential" = "Attributable_death_averted_opt.burden_averted",
             "Associated_daly_high-potential"    = "Associated_daly_averted_opt.burden_averted",
             "Attributable_daly_high-potential"  = "Attributable_daly_averted_opt.burden_averted")
    # -------------------------------------------------------------------------
    avertable_death_inc <-
      data.table(Associated_death_averted_inc   = edit_table(Associated_death_averted_inc)[,1:2],
                 Attributable_death_averted_inc = edit_table(Attributable_death_averted_inc)[,2])
    
    avertable_daly_inc <-
      data.table(Associated_daly_averted_inc    = edit_table(Associated_daly_averted_inc)[,1:2],
                 Attributable_daly_averted_inc  = edit_table(Attributable_daly_averted_inc)[,2])
    
    avertable_burden_inc <- right_join(avertable_death_inc, avertable_daly_inc, 
                                       by=c("Associated_death_averted_inc.Counts"
                                            = "Associated_daly_averted_inc.Counts"))
    avertable_burden_inc <- 
      avertable_burden_inc %>% 
      rename("Counts"                         = "Associated_death_averted_inc.Counts",
             "Associated_death_incremental"   = "Associated_death_averted_inc.burden_averted",
             "Attributable_death_incremental" = "Attributable_death_averted_inc.burden_averted",
             "Associated_daly_incremental"    = "Associated_daly_averted_inc.burden_averted",
             "Attributable_daly_incremental"  = "Attributable_daly_averted_inc.burden_averted")
    
    avertable_burden <- left_join(avertable_burden_con, 
                                  avertable_burden_opt,
                                  by = "Counts") %>%
                        left_join(., avertable_burden_inc, by= "Counts")
  
    return(avertable_burden)
    
  }# end of function -- create_avertable_burden_table
  
# ------------------------------------------------------------------------------
# create graph for vaccine impact by WHO region
  
create_burden_averted_by_region_graph  <- function(Attributable_burden_averted,
                                                   Associated_burden_averted,
                                                   ylim_max,
                                                   ylabel,
                                                   title_name){

Associated_burden_averted$Resistance   <- "Associated with resistance"

Attributable_burden_averted$Resistance <- "Attributable to resistance" 
  
burden_averted_by_region <- rbind(Associated_burden_averted, Attributable_burden_averted)

burden_averted_by_region <- burden_averted_by_region %>% rename("lower_value"  = "2.5%",
                                                                "median_value" = "50%",
                                                                "upper_value"  = "97.5%")

burden_averted_by_region  <- burden_averted_by_region %>% filter(Counts != "Global")
                   
ggplot(burden_averted_by_region, 
       aes(x = reorder(Counts, -median_value), 
           y=median_value, fill=Resistance,
           width=ifelse(Resistance == "Associated with resistance", 0.8, 0.6))) +
  geom_bar(stat = "identity", position="dodge") +
  scale_fill_manual(values = c("lightsteelblue3","lightsteelblue4")) +
  labs(x = "WHO region", y = paste(ylabel)) +
  ylim(0, ylim_max) +
  geom_errorbar(aes(ymin=lower_value, ymax=upper_value), width=0.15,
                size=0.5, position=position_dodge(0)) +
  theme_classic() +
  theme(legend.position = c(0.8, 0.9)) +
  ggtitle(title_name) +
  theme(plot.title = element_text(hjust=-0.1, vjust=3, size = 20)) +
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1))

} # end of function -- create_burden_averted_by_region_graph

# ------------------------------------------------------------------------------



# -------------------------------------------------------------------------
# global vaccine avertable deaths and DALYs attributable to and associated with 
# bacterial antimicrobial resistance by infectious syndrome, 2019

aggregate_impact_by_dp <- function(input_data, 
                                   DiseasePresentation,
                                   mode){
  
  if (mode == "incremental") {
    input_data$va_health_burden <- input_data$va_incre
    
  } else if (mode == "high_potential") {
    input_data$va_health_burden <- input_data$va_high
    
  } else {
    input_data$va_health_burden <- input_data$va_base
  }
  
  impact_by_dp <- input_data %>%
      group_by(Disease_presentation, run_id) %>%
      summarise(averted_burden = sum(va_health_burden), .groups = 'drop')
  
  # -------------------------------------------------------------------------
  
  burden_averted_dp <- data.table("50%"=numeric(), "2.5%"=numeric(), "97.5%"=numeric())
  
  for(i in DiseasePresentation){
    dt <- impact_by_dp %>% filter(Disease_presentation == i)
    dt <- quantile(x = dt$averted_burden, probs = c (0.5, 0.025, 0.975))
    dt <- data.table(t(dt))
    burden_averted_dp <- rbindlist (list (burden_averted_dp, dt),
                                   use.names = FALSE)
  }
  
  burden_averted <- cbind(data.table(Counts=DiseasePresentation, burden_averted_dp))
  
  burden_averted <- burden_averted[, c("Counts", "2.5%", "50%", "97.5%")]
  
  return(burden_averted)} # end of function -- aggregate_impact_by_dp

# -------------------------------------------------------------------------



# -------------------------------------------------------------------------
# create graph for vaccine impact by infectious syndrome

create_burden_averted_by_dp_graph <- function(Attributable_burden_averted,
                                              Associated_burden_averted,
                                              ylim_max,
                                              ylabel,
                                              title_name){
  
  Associated_burden_averted$Resistance   <- "Associated with resistance" 
  
  Attributable_burden_averted$Resistance <- "Attributable to resistance" 
  
  
    burden_averted_by_dp <- rbind(Associated_burden_averted, Attributable_burden_averted)
  
    burden_averted_by_dp <-  burden_averted_by_dp %>% rename("lower_value"  = "2.5%",
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
    
    burden_averted_by_dp$Counts <- gsub("Typhoid, paratyphoid, and iNTS", "TF/PF/iNTS",
                                        burden_averted_by_dp$Counts)
    
    burden_averted_by_dp$Counts <- gsub("Gonorrhoea and chlamydia", "GC/CT",
                                        burden_averted_by_dp$Counts)
    
ggplot(burden_averted_by_dp, 
       aes(x = reorder(Counts, -median_value),
           y = median_value, fill = Resistance,
           width=ifelse(Resistance == "Associated with resistance", 0.8, 0.6))) +

  geom_bar(stat = "identity", position="identity") +
    scale_fill_manual(values = c("lightsteelblue3","lightsteelblue4")) +
    labs(x = "Infectious syndrome", y = paste(ylabel)) + 
    ylim(0,ylim_max) +
    geom_errorbar(aes(ymin=lower_value, ymax=upper_value), width=0.15,
                  size=0.5, position=position_dodge(0)) +
    theme_classic() +
    theme(legend.position = c(0.9, 0.9)) +
    ggtitle(title_name) +
    theme(plot.title = element_text(hjust=-0.05, vjust=3, size = 20)) +
    theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1))

  } # end of function -- create_burden_averted_by_dp_graph

# ------------------------------------------------------------------------------
# global vaccine avertable deaths and DALYs attributable to and associated with 
# bacterial antimicrobial resistance by pathogen, 2019

aggregate_impact_by_pathogen <- function(input_data, 
                                         pathogenlist,
                                         mode){
  
  if (mode == "incremental") {
    input_data$va_health_burden <- input_data$va_incre
    
  } else if (mode == "high_potential") {
    input_data$va_health_burden <- input_data$va_high
    
  } else {
    input_data$va_health_burden <- input_data$va_base
  }
  
    impact_by_pathogen <- input_data %>%
      group_by(Pathogen, run_id) %>%
      summarise(averted_burden = sum(va_health_burden), .groups = 'drop')
    
  # -------------------------------------------------------------------------
  burden_averted_pathogen <- data.table("50%"=numeric(), "2.5%"=numeric(), "97.5%"=numeric())
  
  for(i in pathogenlist){
    dt <- impact_by_pathogen %>% filter(Pathogen == i)
    dt <- quantile(x = dt$averted_burden, probs = c (0.5, 0.025, 0.975))
    dt <- data.table(t(dt))
    burden_averted_pathogen <- rbindlist (list (burden_averted_pathogen, dt),
                                          use.names = FALSE)
  }
  
  Burden_Averted <- cbind(data.table(Counts=pathogenlist, burden_averted_pathogen))
  
  Burden_Averted <- Burden_Averted[,c("Counts", "2.5%", "50%", "97.5%")]
  
  return(Burden_Averted)} # end of function -- aggregate_impact_by_pathogen

# -------------------------------------------------------------------------
# create graph for vaccine impact by pathogen

create_burden_averted_by_pathogen_graph <- function(Attributable_burden_averted,
                                                    Associated_burden_averted,
                                                    ylim_max,
                                                    ylabel,
                                                    title_name){
  
  Associated_burden_averted$Resistance    <- "Associated with resistance" 
  
  Attributable_burden_averted$Resistance  <- "Attributable to resistance" 
  
  
  burden_averted_by_pathogen <- rbind(Associated_burden_averted, Attributable_burden_averted)
  
  burden_averted_by_pathogen <-  burden_averted_by_pathogen %>% rename("lower_value"  = "2.5%",
                                                                       "median_value" = "50%",
                                                                       "upper_value"  = "97.5%")
  
  ggplot(burden_averted_by_pathogen, 
       aes(x = reorder(Counts, -median_value), 
           y = median_value, fill = Resistance,
           width = ifelse(Resistance == "Associated with resistance", 0.8, 0.6))) +
    geom_bar(stat = "identity", position = "identity") +
    scale_fill_manual(values = c("lightsteelblue3","lightsteelblue4")) +
    labs(x = "Pathogen", y = paste(ylabel)) + 
    ylim(0,ylim_max) +
    geom_errorbar(aes(ymin = lower_value, ymax = upper_value), width = 0.15,
                  size = 0.5, position = position_dodge(0)) +
    theme_classic() +
    theme(legend.position = c(0.9, 0.9)) +
    ggtitle(title_name) +
    theme(plot.title = element_text(hjust = -0.05, vjust = 3, size = 20)) +
    theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1))
  
} # end of function -- create_burden_averted_by_pathogen_graph
# ------------------------------------------------------------------------------
# vaccine avertable burden of corresponding vaccines

estimate_burden_averted_add <- function(pathogen,
                                        vaccine_type,
                                        data_input,
                                        mode_input){
  
  burden_add <- data.table(data_input) %>%
    filter(Pathogen == pathogen)
  
  burden_add[, Efficacy:= vaccine_type[,"Efficacy"]]
  burden_add[, Coverage:= vaccine_type[,"Coverage"]]
  burden_add[, Duration:= vaccine_type[,"Duration"]]
  burden_add[, DP:= vaccine_type[,"DP"]]
  burden_add[, VC:= vaccine_type[,"VC"]]
  burden_add[, VO:= vaccine_type[,"VO"]]
  burden_add[, MainAnalysis:= vaccine_type[,"MainAnalysis"]]
  burden_add[Pathogen == "Streptococcus pneumoniae" & 
               Efficacy == "0.58" & 
               Disease_presentation == "LRI and thorax infections",
             Efficacy := 0.25]
  burden_add[Pathogen == "Streptococcus pneumoniae" & 
               Efficacy == "0.7" & 
               Disease_presentation == "LRI and thorax infections",
             Efficacy := 0.5]
  
  burden_add$va_base <- 
    estimate_vaccine_impact(data     = burden_add,
                            scenario = "conservative")[, va_health_burden]
  
  burden_add$va_high <- 
    estimate_vaccine_impact(data     = burden_add,
                            scenario = "optimistic")[, va_health_burden]
  
  burden_add$va_incre <- 
    burden_add$va_high - burden_add$va_base
  
  burden_averted_add <- 
    aggregate_impact_by_pathogen(input_data   = burden_add,
                                 mode         = mode_input,
                                 pathogenlist = pathogen)
  
  return(burden_averted_add)} # end of function -- estimate_burden_averted_add

# -------------------------------------------------------------------------
# estimate vaccine impact on the corresponding pathogen

vaccine_imapct_by_vaccine <- function(data_in,
                                      mode_in){

  vaccine_impact_by_vaccine <- bind_rows(list(
  
  estimate_burden_averted_add(pathogen       = "Acinetobacter baumannii",
                              vaccine_type   = vaccine_profile_dt_add[1,],
                              data_input     = data_in,
                              mode_input     = mode_in),
  
  estimate_burden_averted_add(pathogen       = "Acinetobacter baumannii",
                              vaccine_type   = vaccine_profile_dt_add[2,],
                              data_input     = data_in,
                              mode_input     = mode_in),
  
  estimate_burden_averted_add(pathogen       = "Enterococcus faecium",
                              vaccine_type   = vaccine_profile_dt_add[3,],
                              data_input     = data_in,
                              mode_input     = mode_in),
  
  estimate_burden_averted_add(pathogen       = "Escherichia coli",
                              vaccine_type   = vaccine_profile_dt_add[4,],
                              data_input     = data_in,
                              mode_input     = mode_in),
  
  estimate_burden_averted_add(pathogen       = "Escherichia coli",
                              vaccine_type   = vaccine_profile_dt_add[5,],
                              data_input     = data_in,
                              mode_input     = mode_in),
  
  estimate_burden_averted_add(pathogen       = "Escherichia coli",
                              vaccine_type   = vaccine_profile_dt_add[6,],
                              data_input     = data_in,
                              mode_input     = mode_in),
  
  estimate_burden_averted_add(pathogen       = "Group A Streptococcus",
                              vaccine_type   = vaccine_profile_dt_add[7,],
                              data_input     = data_in,
                              mode_input     = mode_in),
  
  estimate_burden_averted_add(pathogen       = "Haemophilus influenzae",
                              vaccine_type   = vaccine_profile_dt_add[8,],
                              data_input     = data_in,
                              mode_input     = mode_in),
  
  estimate_burden_averted_add(pathogen       = "Klebsiella pneumoniae",
                              vaccine_type   = vaccine_profile_dt_add[9,],
                              data_input     = data_in,
                              mode_input     = mode_in),
  
  estimate_burden_averted_add(pathogen       = "Klebsiella pneumoniae",
                              vaccine_type   = vaccine_profile_dt_add[10,],
                              data_input     = data_in,
                              mode_input     = mode_in),
  
  estimate_burden_averted_add(pathogen       = "Mycobacterium tuberculosis",
                              vaccine_type   = vaccine_profile_dt_add[11,],
                              data_input     = data_in,
                              mode_input     = mode_in),
  
  estimate_burden_averted_add(pathogen       = "Mycobacterium tuberculosis",
                              vaccine_type   = vaccine_profile_dt_add[12,],
                              data_input     = data_in,
                              mode_input     = mode_in),
  
  estimate_burden_averted_add(pathogen       = "Neisseria gonorrhoeae",
                              vaccine_type   = vaccine_profile_dt_add[13,],
                              data_input     = data_in,
                              mode_input     = mode_in),
  
  estimate_burden_averted_add(pathogen       = "Non-typhoidal Salmonella",
                              vaccine_type   = vaccine_profile_dt_add[14,],
                              data_input     = data_in,
                              mode_input     = mode_in),
  
  estimate_burden_averted_add(pathogen       = "Pseudomonas aeruginosa",
                              vaccine_type   = vaccine_profile_dt_add[15,],
                              data_input     = data_in,
                              mode_input     = mode_in),
  
  estimate_burden_averted_add(pathogen       = "Salmonella Paratyphi",
                              vaccine_type   = vaccine_profile_dt_add[16,],
                              data_input     = data_in,
                              mode_input     = mode_in),
  
  estimate_burden_averted_add(pathogen       = "Salmonella Typhi",
                              vaccine_type   = vaccine_profile_dt_add[17,],
                              data_input     = data_in,
                              mode_input     = mode_in),
  
  estimate_burden_averted_add(pathogen       = "Shigella spp.",
                              vaccine_type   = vaccine_profile_dt_add[18,],
                              data_input     = data_in,
                              mode_input     = mode_in),
  
  estimate_burden_averted_add(pathogen       = "Staphylococcus aureus",
                              vaccine_type   = vaccine_profile_dt_add[19,],
                              data_input     = data_in,
                              mode_input     = mode_in),
  
  estimate_burden_averted_add(pathogen       = "Streptococcus pneumoniae",
                              vaccine_type   = vaccine_profile_dt_add[20,],
                              data_input     = data_in,
                              mode_input     = mode_in),
  
  estimate_burden_averted_add(pathogen       = "Streptococcus pneumoniae",
                              vaccine_type   = vaccine_profile_dt_add[21,],
                              data_input     = data_in,
                              mode_input     = mode_in)))
  
  vaccine_impact_by_vaccine$Counts <- vaccine_profile_dt_add$Vaccine_pathogen

  return(vaccine_impact_by_vaccine)} # end of function -- vaccine_imapct_by_vaccine

# -------------------------------------------------------------------------
create_burden_averted_by_vp_graph <- function(Attributable_burden_averted,
                                              Associated_burden_averted,
                                              ylim_max,
                                              ylabel,
                                              title_name){
  
  Associated_burden_averted$Resistance    <- "Associated with resistance" 
  
  Attributable_burden_averted$Resistance  <- "Attributable to resistance" 
  
  burden_averted_by_pathogen <- rbind(Associated_burden_averted, Attributable_burden_averted)
  
  burden_averted_by_pathogen <-  burden_averted_by_pathogen %>% rename("lower_value"  = "2.5%",
                                                                       "median_value" = "50%",
                                                                       "upper_value"  = "97.5%")
  
  ggplot(burden_averted_by_pathogen, 
         aes(x = reorder(Counts, -median_value), 
             y = median_value, fill = Resistance,
             width = ifelse(Resistance == "Associated with resistance", 0.8, 0.6))) +
    geom_bar(stat = "identity", position = "identity") +
    scale_fill_manual(values = c("lightsteelblue3","lightsteelblue4")) +
    labs(x = "Vaccine", y = paste(ylabel)) + 
    ylim(0,ylim_max) +
    geom_errorbar(aes(ymin = lower_value, ymax = upper_value), width = 0.15,
                  size = 0.5, position = position_dodge(0)) +
    theme_classic() +
    theme(legend.position = c(0.9, 0.9)) +
    ggtitle(title_name) +
    theme(plot.title = element_text(hjust = -0.05, vjust = 3, size = 20)) +
    theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1))
  
} # end of function -- create_burden_averted_by_vaccine_profile_graph
# ------------------------------------------------------------------------------



# ------------------------------------------------------------------------------
# further analysis -- vaccine avertable burden of existing vaccines

estimate_existing_vaccine_impact <- function(input_data){
  
  va <- input_data %>%
    filter(Pathogen == "Haemophilus influenzae" | 
           Pathogen == "Streptococcus pneumoniae")
  
  va <- data.table(va)
  
  va[, va_age := 0]
  
  # using existing Hib vaccine efficacy
  va[Pathogen  == "Haemophilus influenzae" & Age_group == "PN",
     va_age := burden_psa * (4/48 * hib_vaccine_coverage_2019 * 0.59 +
                             4/48 * hib_vaccine_coverage_2019 * 0.92 + 
                            37/48 * hib_vaccine_coverage_2019 * 0.93) * 0.95]
  
  va[Pathogen  == "Haemophilus influenzae" & Age_group == "1 to 4",
     va_age := burden_psa * 0.93 * hib_vaccine_coverage_2018 * 0.95]
  
  # using existing PCV efficacy
  va[Pathogen  == "Streptococcus pneumoniae" & Age_group == "PN",
     va_age := burden_psa * (4/48 * pcv_vaccine_coverage_2019 * 0.29 + 
                             4/48 * pcv_vaccine_coverage_2019 * 0.58 + 
                            37/48 * pcv_vaccine_coverage_2019 * 0.58)]
  
  va[Pathogen  == "Streptococcus pneumoniae" & Age_group == "PN" &
     Disease_presentation == "LRI and thorax infections" ,
     va_age := burden_psa * (4/48 * pcv_vaccine_coverage_2019 * 0.25 + 
                             4/48 * pcv_vaccine_coverage_2019 * 0.25 + 
                            37/48 * pcv_vaccine_coverage_2019 * 0.25)]
  
  va[Pathogen  == "Streptococcus pneumoniae" & Age_group == "1 to 4",
     va_age := burden_psa * 0.58 * pcv_vaccine_coverage_2018]
  
  va[Pathogen  == "Streptococcus pneumoniae" & Age_group == "1 to 4" & 
     Disease_presentation == "LRI and thorax infections",
     va_age := burden_psa * 0.25 * pcv_vaccine_coverage_2018]
  
  # applying vaccine target infectious syndrome
  vaccine_impact <- va
  
  vaccine_impact[, va_health_burden := 0]
  
  vaccine_impact[DP  == "All" |
                   
                   (DP  == "BSI, CNS infections, Cardiac infections, LRI" & 
                      (Disease_presentation == "BSI" | 
                         Disease_presentation == "CNS infections" | 
                         Disease_presentation == "Cardiac infections" |
                         Disease_presentation == "LRI and thorax infections")),
                 
                 va_health_burden := va_age]
  
  # aggregate impact by run_id
  impact_by_pathogen <- vaccine_impact %>%
    group_by(Pathogen, run_id) %>%
    summarise(averted_burden = sum(va_health_burden), .groups = 'drop')
  
  # aggregate impact by pathogen
  burden_averted_pathogen <- data.table("50%"=numeric(), "2.5%"=numeric(), "97.5%"=numeric())
  
  existingvaccine <- c("Haemophilus influenzae", "Streptococcus pneumoniae")
  
  for(i in existingvaccine){
    dt <- impact_by_pathogen %>% filter(Pathogen == i)
    dt <- quantile(x = dt$averted_burden, probs = c (0.5, 0.025, 0.975))
    dt <- data.table(t(dt))
    burden_averted_pathogen <- rbindlist (list (burden_averted_pathogen, dt),
                                          use.names = FALSE)
  }
  
  Burden_Averted <- cbind(data.table(Counts=existingvaccine, burden_averted_pathogen))
  
  Burden_Averted <- Burden_Averted[,c("Counts", "2.5%", "50%", "97.5%")]
  
  return(Burden_Averted)} # end of function -- estimate_existing_vaccine_impact

# ------------------------------------------------------------------------------
# Appendix -- vaccine avertable burdens by infectious syndrome and pathogen

aggregate_impact_by_dp_pathogen <- function(input_data, 
                                            input_rep,
                                            mode){
  
  input_data$va_health_burden <- input_data$va_base

  impact_by_dp_p <- input_data %>%
      group_by(Disease_presentation, Pathogen, run_id) %>%
      summarise(averted_burden = sum(va_health_burden), .groups = 'drop')

  impact_by_dp_p$number <- rep(input_rep, each=run)
  
  # -------------------------------------------------------------------------
  
  burden_averted_dp    <- data.table("50%"=numeric(), "2.5%"=numeric(), "97.5%"=numeric())
  
  for(i in input_rep){
    dt <- impact_by_dp_p %>% filter(number == i)
    dt <- quantile(x = dt$averted_burden, probs = c (0.5, 0.025, 0.975))
    dt <- data.table(t(dt))
    burden_averted_dp <- rbindlist (list (burden_averted_dp, dt),
                                    use.names = FALSE)
  }
  
  burden_averted <- cbind(data.table(impact_by_dp_p[impact_by_dp_p$run_id == "1",1:2], 
                                     burden_averted_dp))
  
  return(burden_averted)} # end of function -- aggregate_impact_by_dp_pathogen

# ------------------------------------------------------------------------------
# create graph of vaccine impact by infectious syndrome and pathogen

burden_averted_by_dp_pat <- function(data_input, image_file){
  
  data_input <- data_input[, 1:3] 
  
  colnames(data_input) <- 
    c("Disease_presentation", "Pathogen", "median")
  
  data_input <- data_input %>%
    filter(median != 0)
  
  data_input$Disease_presentation <- gsub("BSI", "Bloodstream Infections", data_input$Disease_presentation)
  data_input$Disease_presentation <- gsub("TB", "Tuberculosis", data_input$Disease_presentation)
    
  BSI <- data_input %>%
    filter(Disease_presentation == "Bloodstream Infections") %>%
    arrange(median)
  
  BSI$Pathogen <- gsub("Group A Streptococcus", "Others", BSI$Pathogen)
  
  BSI$Pathogen <- gsub("Non-typhoidal Salmonella", "Others", BSI$Pathogen)
  
  BSI$Pathogen <- gsub("Salmonella Typhi", "Others", BSI$Pathogen)
  
  BSI$Pathogen <- gsub("Enterococcus faecium", "Others", BSI$Pathogen)
  
  BSI$Pathogen <- gsub("Pseudomonas aeruginosa", "Others", BSI$Pathogen)
  
  BSI <- BSI %>% group_by(Disease_presentation, Pathogen) %>%
    summarise(median=sum(median), .groups = 'drop')
  
  LRI_TB <- data_input %>%
    filter(Disease_presentation == "LRI and thorax infections"|
           Disease_presentation == "Tuberculosis")
  
  others <- data.frame("Disease_presentation" = "Others",
                       "Pathogen"             = " ",
                       "median"               = sum(data_input$median) 
                       - sum(BSI$median) - sum(LRI_TB$median))
  
  dp_by_pathogen <- rbind(BSI, LRI_TB, others)
  
  dp_by_pathogen <- dp_by_pathogen %>% arrange(Disease_presentation, median)
  
  png(image_file, 
      width = 10, height = 10, res = 1000, units="in")
  
  PieDonut(dp_by_pathogen, aes(pies = Disease_presentation, donuts = Pathogen, count = median),
           showPieName=FALSE, showRatioThreshold = 0, start = 4.2,
           r0 = 0, r1 = 0.5, r2 = 0.6,
           pieAlpha = 0.7, donutAlpha = 0.7,
           labelpositionThreshold=1,
           titlesize = 4.5)
  
  dev.off()
} # end of function -- burden_averted_by_dp_pat
# ------------------------------------------------------------------------------