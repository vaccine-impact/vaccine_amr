# ------------------------------------------------------------------------------
# functions.R
#
# functions for analysis to estimate vaccine averted AMR health burden
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# create data table of (194) WHO countries classified by (6) WHO regions 
# and (4) World Bank income levels
create_country_table <- function (regional_classification_file,
                                  country_income_classification_file) {
  
  # open files for WHO country-region and WB country-income classifications
  regional_dt <- fread (file = regional_classification_file)
  income_dt   <- fread (file = country_income_classification_file)
  
  # extract (194) WHO countries
  country_dt <- regional_dt [WHO_region != "Not Classified"]
  
  # extract columns for income group
  income_dt <- income_dt [, c("iso3_code", "World_Bank_income_group")]
  
  # add column for World Bank income level
  country_dt <- merge (x     = country_dt, 
                       y     = income_dt, 
                       by.x  = "iso3_code",
                       by.y  = "iso3_code", 
                       all.x = TRUE,
                       all.y = FALSE)
  
  # return data table of of (194) WHO countries classified by (6) WHO regions 
  # and (4) World Bank income levels
  return (country_dt)
  
} # end of function -- create_country_table
# ------------------------------------------------------------------------------
