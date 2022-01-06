# ------------------------------------------------------------------------------
# vaccine impact estimates
# ------------------------------------------------------------------------------
# load libraries
library (readr)
library (dplyr)
library (ggplot2)

# clear workspace
rm (list = ls ())
# ------------------------------------------------------------------------------
# IHME data on the burden of drug resistance Tuberculosis
# location: select only countries and territories
# year: 2019
# sex: both
# matrix: number
# multidrug-resistant tuberculosis without extensive drug resistance, extensively drug-resistant Tuberculosis

# import data set
IHME_tuberculosis <- read_csv("~/GitHub/vaccine_amr/data/IHME/IHME-GBD_2019_DATA-c47395f9-1.csv")

IHME_iso3_codes <- read_csv("~/GitHub/vaccine_amr/data/IHME/IHME_iso3_codes.csv")
IHME_iso3_codes <- IHME_iso3_codes[,c("location_id", "ISO3")]

demogrphic <- read_csv("~/GitHub/vaccine_amr/data/population.csv")

# data cleaning
IHME_tuberculosis <- left_join(IHME_tuberculosis, IHME_iso3_codes, by=c("location_id"="location_id"))

IHME_tuberculosis <- IHME_tuberculosis[c(4<IHME_tuberculosis$age_id & IHME_tuberculosis$age_id<22|
                                           IHME_tuberculosis$age_id & IHME_tuberculosis$age_id==28), ] # age group

IHME_tuberculosis <- IHME_tuberculosis[,c("ISO3","measure_name", "age_id", "age_name", "cause_id", "cause_name", "val")]

# ------------------------------------------------------------------------------
# multidrug-resistant tuberculosis without extensive drug resistance + extensively drug-resistant Tuberculosis
resistant_tuberculosis <- IHME_tuberculosis %>%
  group_by(ISO3, measure_name, age_id, age_name) %>%
  summarise(val = sum(val))

resistant_tuberculosis$age_name <- factor(resistant_tuberculosis$age_name, 
                                          levels=c("<1 year","Post Neonatal","1 to 4","5 to 9","10 to 14","15 to 19","20 to 24","25 to 29","30 to 34","35 to 39","40 to 44","45 to 49", 
                                                                 "50 to 54","55 to 59","60 to 64","65 to 69","70 to 74","75 to 79","80 plus"), order=T)
# tuberculosis vaccine number 2
# target population: all ages
tb_vaccine_coverage <- 0.7
tb_vaccine_efficacy <- 0.8
tb_vaccine_duration <- 10
# ------------------------------------------------------------------------------
# vaccine impact on the resistant tuberculosis incidence on 2019
resistant_tuberculosis_incidence <- resistant_tuberculosis[resistant_tuberculosis$measure_name=="Incidence",]

# vaccine averted resistant tuberculosis incidence of all countries
resistant_tuberculosis_incidence <- resistant_tuberculosis_incidence %>% 
                                group_by(age_name) %>%
                                summarise(nonvaccinated_incidence = sum(val)) %>%
                                mutate(vaccinated_incidence = nonvaccinated_incidence * (1 - tb_vaccine_coverage * tb_vaccine_efficacy)) %>%
                                mutate(vaccine_averted_incidence = nonvaccinated_incidence - vaccinated_incidence)

vaccine_averted_incidence_2019 <- sum(resistant_tuberculosis_incidence$vaccine_averted_incidence) 

# graph for the number of resistant tuberculosis incidences by age in 2019.
p1 <- ggplot(resistant_tuberculosis_incidence, aes(x=age_name)) +
  geom_point(aes(y=vaccinated_incidence, group=1, colour="vaccinated")) +
  geom_line(aes(y=vaccinated_incidence, group=1, colour="vaccinated")) +
  geom_point(aes(y=nonvaccinated_incidence, group=1, colour="non-vaccinated")) +
  geom_line(aes(y=nonvaccinated_incidence, group=1, colour="non-vaccinated")) +
  ylab("Number of resistant tuberculosis incidences") + 
  xlab("age") +
  theme_classic() +
  expand_limits(y=0)
p1
# ------------------------------------------------------------------------------
# vaccine impact on the resistant tuberculosis death on 2019
resistant_tuberculosis_death <- resistant_tuberculosis[resistant_tuberculosis$measure_name=="Deaths",]

# graph for the all countries
resistant_tuberculosis_death <- resistant_tuberculosis_death %>% 
  group_by(age_name) %>%
  summarise(nonvaccinated_death = sum(val)) %>%
  mutate(vaccinated_death = nonvaccinated_death * (1 - tb_vaccine_coverage * tb_vaccine_efficacy)) %>%
  mutate(vaccine_averted_death = nonvaccinated_death - vaccinated_death)

vaccine_averted_death_2019 <- sum(resistant_tuberculosis_death$vaccine_averted_death) 

# graph for the number of resistant tuberculosis deaths by age in 2019.
p2 <- ggplot(resistant_tuberculosis_death, aes(x=age_name)) +
  geom_point(aes(y=vaccinated_death, group=1, colour="vaccinated")) +
  geom_line(aes(y=vaccinated_death, group=1, colour="vaccinated")) +
  geom_point(aes(y=nonvaccinated_death, group=1, colour="non-vaccinated")) +
  geom_line(aes(y=nonvaccinated_death, group=1, colour="non-vaccinated")) +
  ylab("Number of resistant tuberculosis deaths") + 
  xlab("age") +
  theme_classic() +
  expand_limits(y=0)
p2
# ------------------------------------------------------------------------------
# vaccine mpact on the resistant tuberculosis DALY on 2019
resistant_tuberculosis_DALY <- resistant_tuberculosis[resistant_tuberculosis$measure_name=="DALYs (Disability-Adjusted Life Years)",]

# graph for the all countries
resistant_tuberculosis_DALY <- resistant_tuberculosis_DALY %>% 
  group_by(age_name) %>%
  summarise(nonvaccinated_DALY = sum(val)) %>%
  mutate(vaccinated_DALY = nonvaccinated_DALY * (1 - tb_vaccine_coverage * tb_vaccine_efficacy)) %>%
  mutate(vaccine_averted_DALY = nonvaccinated_DALY - vaccinated_DALY)

vaccine_averted_DALY_2019 <- sum(resistant_tuberculosis_DALY$vaccine_averted_DALY) 

# graph for the number of resistant tuberculosis DALY by age in 2019.
p3 <- ggplot(resistant_tuberculosis_DALY, aes(x=age_name)) +
  geom_point(aes(y=vaccinated_DALY, group=1, colour="vaccinated")) +
  geom_line(aes(y=vaccinated_DALY, group=1, colour="vaccinated")) +
  geom_point(aes(y=nonvaccinated_DALY, group=1, colour="non-vaccinated")) +
  geom_line(aes(y=nonvaccinated_DALY, group=1, colour="non-vaccinated")) +
  ylab("Number of resistant tuberculosis DALY") + 
  xlab("age") +
  theme_classic() +
  expand_limits(y=0)
p3