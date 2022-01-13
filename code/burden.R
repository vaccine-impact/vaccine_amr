# ------------------------------------------------------------------------------
# vaccine impact estimates
# ------------------------------------------------------------------------------
# load libraries
library (readr)
library (dplyr)
library (ggplot2)
library (reshape2)
library (tidyverse)

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
# target population: all ages
tb_vaccine_coverage <- 0.7
tb_vaccine_efficacy <- 0.8
tb_vaccine_duration <- 10

# data frame for the AMR tuberculosis health burden 
resistant_tb <- IHME_tuberculosis %>%
  group_by(ISO3, measure_name, age_id, age_name) %>%
  summarise(val = sum(val)) %>%   # multidrug-resistant tuberculosis + extensively drug-resistant Tuberculosis
  mutate(vaccinated_val = val * (1 - tb_vaccine_coverage * tb_vaccine_efficacy)) %>%
  mutate(vaccine_averted_val = val - vaccinated_val)


# Chage age to 1 to 99
resistant_tb$age_name <- factor(resistant_tb$age_name, 
                                          levels=c("<1 year","1 to 4","5 to 9","10 to 14","15 to 19","20 to 24","25 to 29","30 to 34","35 to 39","40 to 44","45 to 49", 
                                                                 "50 to 54","55 to 59","60 to 64","65 to 69","70 to 74","75 to 79","80 plus"), order=T)

resistant_tb <- resistant_tb %>% arrange(ISO3, measure_name, age_name) 

resistant_tb <- mutate(resistant_tb, frequency = if(age_name == "<1 year"){frequency = 1}
                 else if (age_name == "1 to 4"){frequency = 4}
                 else if (age_name == "80 plus"){frequency = 20}
                 else {frequency = 5})

resistant_tb <- resistant_tb[rep(row.names(resistant_tb), resistant_tb$frequency), 1:8]

resistant_tb$age <- rep(0:99, 612)

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# example graph: vaccine averted resistant burden - Sierra Leone
resistant_tb_SLE <- resistant_tb %>%
  filter(ISO3 == "SLE" & measure_name == "Incidence")

# vaccine averted resistant burden
ggplot(resistant_tb_SLE, aes(x=age)) +
  geom_line(aes(y=vaccinated_val, group=1, colour="vaccinated")) +
  geom_line(aes(y=val, group=1, colour="non-vaccinated")) +
  ylab("Number of vaccine averted incidence per 100,000") + 
  xlab("age") +
  theme_classic() +
  expand_limits(y=0)
# ------------------------------------------------------------------------------
# saving cvs file: vaccine averted resistant burden
resistant_tb_burden <- resistant_tb %>% group_by(measure_name, ISO3) %>%
  summarise(vaccine_averted_val = sum(vaccine_averted_val))

resistant_tb_incidence <- resistant_tb_burden[resistant_tb_burden$measure_name=="Incidence",]
resistant_tb_deaths <- resistant_tb_burden[resistant_tb_burden$measure_name=="Deaths",]
resistant_tb_DALYs <- resistant_tb_burden[resistant_tb_burden$measure_name=="DALYs (Disability-Adjusted Life Years)",]

fwrite (x    = resistant_tb_incidence, 
        file = "resistant_tb_incidence.csv")

resistant_tb_incidence$vaccine_averted_val

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

# Create map of incidence averted (AMR tuberculosis)

# load libraries
library (ggplot2)
library (sf)
library (rnaturalearth)
library (rnaturalearthdata)
library (rgeos)
library (data.table)
library (ggforce)
library (ggpubr)
library (stats)

# ------------------------------------------------------------------------------

create_map <- function (burden_file, 
                        map_file) {
  
  # read results for dalys averted
  table <- fread (file = burden_file)
  
  # map tutorial
  world <- ne_countries (scale       = "medium", 
                         returnclass = "sf")
  setDT (world)
  setkey (world, iso_a3)
  
  # combine tables to add geometry
  dt <- merge (x    = table, 
               y    = world, 
               by.x = "ISO3", 
               by.y = "iso_a3", 
               all  = F ) # why it doesn't work when all=T or all.x=T
  
  # generate map of dalys averted per 1000 vaccinated individuals
  plot <- ggplot (data = dt) +
    geom_sf (aes (fill = vaccine_averted_val, geometry = geometry)) + 
    scale_fill_viridis_c (option = "viridis", direction = -1, na.value = "grey90") +
    # scale_fill_viridis_c (option = "plasma", direction = -1, na.value = "grey90") +
    ggtitle ("DALYS averted per 1000 fully vaccinated individuals") + 
    theme (legend.title = element_blank()) + 
    theme (axis.text.x = element_blank(), axis.ticks = element_blank()) + 
    theme (axis.text.y = element_blank(), axis.ticks = element_blank()) + 
    theme (panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    theme (plot.title = element_text(size = 10))
  
  # save plot
  ggsave (filename = map_file, 
          width = 6, 
          height = 3, 
          dpi = 600)
  
  return ()
  
} # end of function -- create_map
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# create map of daly results
create_map (burden_file = "resistant_tb_incidence.csv", 
            map_file  = "IncidenceAvertedPer100,000.png")
# ------------------------------------------------------------------------------