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
  
  
} # end of function -- create_country_table
# ------------------------------------------------------------------------------