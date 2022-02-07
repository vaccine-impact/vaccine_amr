# [table in appendix] vaccine avertable health burdens associated with and attributable
# to AMR by WHO region, pathogen, disease presentation, and age group

head(combined_dt)

a <- combined_dt[, c("WHO_region", "Disease_presentation", "Age_group", "Pathogen", 
                     "Efficacy", "Coverage", "Duration", "DiseasePresentation", "Vaccination", "Associated_resistant_mean")]
a <- a %>% rename("burden_psa" = "Associated_resistant_mean")
a <- estimate_vaccine_impact("any", data = a)
a <- a$va_health_burden

b <- combined_dt[, c("WHO_region", "Disease_presentation", "Age_group", "Pathogen", 
                     "Efficacy", "Coverage", "Duration", "DiseasePresentation", "Vaccination", "Associated_resistant_lower")]
b <- b %>% rename("burden_psa" = "Associated_resistant_lower")
b <- estimate_vaccine_impact("any", data = b)
b <- b$va_health_burden

c <- combined_dt[, c("WHO_region", "Disease_presentation", "Age_group", "Pathogen", 
                     "Efficacy", "Coverage", "Duration", "DiseasePresentation", "Vaccination", "Associated_resistant_upper")]
c <- c %>% rename("burden_psa" = "Associated_resistant_upper")
c <- estimate_vaccine_impact("any", data = c)
c <- c$va_health_burden

d <- combined_dt[, c("WHO_region", "Disease_presentation", "Age_group", "Pathogen", 
                     "Efficacy", "Coverage", "Duration", "DiseasePresentation", "Vaccination", "Attributable_resistance_mean")]
d <- d %>% rename("burden_psa" = "Attributable_resistance_mean")
d <- estimate_vaccine_impact("any", data = d)
d <- d$va_health_burden

e <- combined_dt[, c("WHO_region", "Disease_presentation", "Age_group", "Pathogen", 
                     "Efficacy", "Coverage", "Duration", "DiseasePresentation", "Vaccination", "Attributable_resistance_lower")]
e <- e %>% rename("burden_psa" = "Attributable_resistance_lower")
e <- estimate_vaccine_impact("any", data = e)
e <- e$va_health_burden

f <- combined_dt[, c("WHO_region", "Disease_presentation", "Age_group", "Pathogen", 
                     "Efficacy", "Coverage", "Duration", "DiseasePresentation", "Vaccination", "Attributable_resistance_upper")]
f <- f %>% rename("burden_psa" = "Attributable_resistance_upper")
f <- estimate_vaccine_impact("any", data = f)
f <- f$va_health_burden

x <- data.frame("vaccine avertable deaths associated with resistance (mean)" = a,
                "vaccine avertable deaths associated with resistance (lower)" = b,
                "vaccine avertable deaths associated with resistance (upper)" = c,
                "vaccine avertable deaths attributable to resistance (mean)" = d,
                "vaccine avertable deaths attributable to resistance (lower)" = e,
                "vaccine avertable deaths attributable to resistance (upper)" = f)

AMR_burden_data_updated <- cbind(death_burden_dt, x)

AMR_burden_data_updated <- AMR_burden_data_updated[, -c(5:7)]

AMR_burden_data_updated <- AMR_burden_data_updated %>%
  rename("deaths associated with resistance (mean)" =  "Associated_resistant_mean",
         "deaths associated with resistance (lower)"  =  "Associated_resistant_lower",
         "deaths associated with resistance (upper)"  =  "Associated_resistant_upper",
         "deaths attributable to resistance (mean)" =  "Attributable_resistance_mean",
         "deaths attributable to resistance (lower)" =  "Attributable_resistance_lower",
         "deaths attributable to resistance (upper)" =  "Attributable_resistance_upper")

# save as xlsx
fwrite (x    = AMR_burden_data_updated,
        file = "AMR burden data_updated.csv")
