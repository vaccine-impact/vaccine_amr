# ------------------------------------------------------------------------------
# population data from UN WPP (World Population Prospects)
#
# 1950 - 2100
# ------------------------------------------------------------------------------
# load libraries
library (dplyr)
library (readxl)
library (reshape2)
library (data.table)

# move to base directory (run code from source directory)
# setwd("~/GitHub/vaccine_amr/code")
source_wd <- getwd ()
setwd ("../")
# ------------------------------------------------------------------------------

# create demography 
create_demography <- function (country_dt,
                               pop_file) {
  # "WPP2019_INT_F03_1_POPULATION_BY_AGE_ANNUAL_BOTH_SEXES.xlsx" data cleaning
  population_file_a <- read_excel("~/GitHub/vaccine_amr/data/WPP2019_INT_F03_1_POPULATION_BY_AGE_ANNUAL_BOTH_SEXES.xlsx", 
                                  col_names = FALSE)
  population_file_b  <- read_excel("~/GitHub/vaccine_amr/data/WPP2019_INT_F03_1_POPULATION_BY_AGE_ANNUAL_BOTH_SEXES.xlsx", 
                                   sheet = "MEDIUM VARIANT", col_names = FALSE)
  population_file <- rbind(population_file_a, population_file_b)
  
  names(population_file) <- lapply(population_file[13, ], as.character)
  population_file <- population_file %>% filter(Type=="Country/Area")
  
  # Removing duplicates (year of 2020)
  population_file <- population_file[!(population_file$Variant == "Medium variant" & population_file$"Reference date (as of 1 July)" == "2020"), ]
  
  # ------------------------------------------------------------------------------
  # "WPP2019_F01_LOCATIONS.XLSX" data cleaning
  WPP_iso3 <- read_excel("~/GitHub/vaccine_amr/data/WPP2019_F01_LOCATIONS.XLSX", 
                         col_names = FALSE)
  names(WPP_iso3) <- lapply(WPP_iso3[13, ], as.character)
  WPP_iso3 <- WPP_iso3[-c(1:13),c(2,5)]
  
  # ------------------------------------------------------------------------------
  # adding ISO3 to WPP2019_POPULATION data frame
  population_iso3 <- left_join(population_file, WPP_iso3, by=c("Region, subregion, country or area *"="Region, subregion, country or area*"))
  
  population_iso3 <- population_iso3 %>% relocate(c("ISO3 Alpha-code"))
  population_iso3 <- population_iso3[,-c(2:8)]
  
  # ------------------------------------------------------------------------------
  # reshaping and rearranging the data frame (201 countries)
  population_iso3 <- melt(as.data.frame(population_iso3), id = c("ISO3 Alpha-code", "Reference date (as of 1 July)"))
  names(population_iso3) <- c("iso3_code","year","age","population")
  population_iso3 <- population_iso3%>%arrange(iso3_code,year,age)
  
  # ------------------------------------------------------------------------------
  # WHO 194 countries only
  WHO194 <- country_dt[,c(1)]
  WHO_poulation <- right_join(population_iso3, WHO194, by=c("iso3_code"="iso3_code"))
  View(WHO_poulation)
  
  fwrite (x = WHO_poulation, 
          file = pop_file)
  
  # go back to source directory
  setwd (source_wd)
  
  return ()
  
} # end of function -- create_demography

# create demography 
create_demography (country_dt = country_dt,
                   pop_file   = "data/population.csv")
