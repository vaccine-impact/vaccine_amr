# ------------------------------------------------------------------------------
# vaccine impact estimates using IHME AMR burden (Deaths)
# ------------------------------------------------------------------------------
# load libraries
library (readr)
library (readxl)
library (dplyr)
library (ggplot2)
library (reshape2)
library (tidyverse)
library (data.table)

# clear workspace
rm (list = ls ())

# set working directory
setwd("~/GitHub/vaccine_amr/data")
# ------------------------------------------------------------------------------
## IHME data on the burden ##
# WHO region: Africa, Americas, Eastern Mediterranean, Europe, South-East Asia, 
#             Western Pacific, Unclassified 
# year: ?
# sex: both
# ------------------------------------------------------------------------------
# Import IME AMR burden

IHME_AMR_burden <- read_excel("IHME AMR burden.xlsx", 
                              col_names = FALSE)
# Data cleaning
names(IHME_AMR_burden) <- c("WHO_region", "Disease_presentation","Age_group","Pathogen", 
                            "Associated_suscetible_mean","Associated_suscetible_lower",
                            "Associated_suscetible_upper","Associated_resistant_mean", 
                            "Associated_resistant_lower","Associated_resistant_upper",
                            "Attributable_resistance_mean","Attributable_resistance_lower", 
                            "Attributable_resistance_upper")


IHME_AMR_burden <- IHME_AMR_burden[-c(1:2),]

IHME_AMR_burden[,5:13] <- lapply(IHME_AMR_burden[,5:13], as.numeric)

IHME_AMR_burden$Age_group <- factor(IHME_AMR_burden$Age_group, 
                                levels=c(as.vector(unlist(IHME_AMR_burden[1:23,"Age_group"]))), order=T)

levels(IHME_AMR_burden$Age_group) <- c("EN", "LN", "PN",  "1 to 4", "5 to 9","10 to 14",
                                       "15 to 19", "20 to 24"  ,"25 to 29",  "30 to 34",
                                       "35 to 39", "40 to 44", "45 to 49",  "50 to 54" ,      
                                       "55 to 59", "60 to 64", "65 to 69", "70 to 74",
                                       "75 to 79", "80 to 84", "85 to 89", "90 to 94", "95 plus")      
# ------------------------------------------------------------------------------
# Import WHO vaccine profile
Vaccine_profile <- read_excel("Vaccine profile for IHME burden.xlsx", 
                              sheet = "Vaccine profile assumptions")

Vaccine_profile <- Vaccine_profile[,c(3, 4, 5, 6, 9, 10,12)]

Vaccine_profile <- Vaccine_profile %>%
  rename("Efficacy" = "Efficacy (%)",
         "Coverage" = "Coverage in target group",
         "Duration" = "Duration of protection (year)", 
         "DiseasePresentation" = "Disease presentation")

Vaccine_profile$Efficacy <- Vaccine_profile$Efficacy * 1/100
Vaccine_profile$Coverage <- Vaccine_profile$Coverage * 1/100

Vaccine_profile <- Vaccine_profile %>% filter(Selection == "Yes")

fwrite (x    = Vaccine_profile, 
        file = "Vaccine_profile.csv")
# ------------------------------------------------------------------------------
# Key element vectors

reference <- IHME_AMR_burden %>% count(Pathogen, Disease_presentation)
reference

pathogenlist <- IHME_AMR_burden %>% count(Pathogen)
pathogenlist <- as.vector(unlist(pathogenlist[,"Pathogen"]))

burden <- IHME_AMR_burden %>% count(Disease_presentation)
burden <- as.vector(unlist(burden[,"Disease_presentation"]))

WHO_region <- IHME_AMR_burden %>% count(WHO_region)
WHO_region <- as.vector(unlist(WHO_region[,"WHO_region"]))

burden_type <- names(IHME_AMR_burden)[5:13]
# ------------------------------------------------------------------------------

## the death trend across all age groups

GlobalTrend <- function(pathogen){
  x <- IHME_AMR_burden[IHME_AMR_burden$Pathogen == pathogen,]
  
  x <- left_join(x, Vaccine_profile, by=c("Pathogen" = "Pathogen"))
  
  x <- x %>% 
    group_by(Pathogen, Efficacy, Coverage, Age_group) %>%
    summarise(sum_Attributable_resistance_mean=sum(Attributable_resistance_mean), .groups = 'drop') 
  
  x <- x %>% 
    mutate(v_Attributable_resistance_mean = sum_Attributable_resistance_mean * Efficacy * Coverage,
           va_Attributable_resistance_mean = v_Attributable_resistance_mean-sum_Attributable_resistance_mean)
  
  ggplot(x, aes(x=Age_group)) +
    geom_line(aes(y=sum_Attributable_resistance_mean, group=1, colour="non-vaccinated")) +
    ylab("Number of Vaccine Attributable Deaths") + 
    xlab("Age") +
    ggtitle(paste(pathogen,"(global)")) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5)) +
    expand_limits(y=0)
  
  # save plot
  ggsave (filename = paste(pathogen,"_global.png"),
          width = 15, 
          height = 6, 
          dpi = 600
          )
}

lapply(pathogenlist, GlobalTrend)

# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------

burden_vaccine_profile <- left_join(IHME_AMR_burden, Vaccine_profile, by=c("Pathogen" = "Pathogen"))
  
burden_vaccine_profile <- burden_vaccine_profile %>%
    mutate(v_Attributable_resistance_mean 
           = ifelse(Vaccination=="4 week (effective at 6 week)" & Duration=="5" & Age_group == "Post Neonatal", Attributable_resistance_mean * Efficacy * Coverage * 47/48,
                    ifelse(Vaccination=="4 week (effective at 6 week)" & Duration=="5" & Age_group == "1 to 4", Attributable_resistance_mean * Efficacy * Coverage,
                    ifelse(Vaccination=="4 week (effective at 6 week)" & Duration=="5" & Age_group == "5 to 9", Attributable_resistance_mean * Efficacy * Coverage * 1/5*5/52,
                           
             ifelse(Vaccination=="4 week (effective at 6 week)" & Duration=="2" & Age_group == "Post Neonatal", Attributable_resistance_mean * Efficacy * Coverage * 47/48,
                    ifelse(Vaccination=="4 week (effective at 6 week)" & Duration=="2" & Age_group == "1 to 4", Attributable_resistance_mean * Efficacy * Coverage * 1/4*5/52,
                           
             ifelse(Vaccination=="4 week (effective at 6 week), 60 or above" & Duration=="2" & Age_group == "Post Neonatal", Attributable_resistance_mean * Efficacy * Coverage * 47/48,
                    ifelse(Vaccination=="4 week (effective at 6 week), 60 or above" & Duration=="2" & Age_group == "1 to 4", Attributable_resistance_mean * Efficacy * Coverage * 1/4*5/52,
                    ifelse(Vaccination=="4 week (effective at 6 week), 60 or above" & Duration=="2" & Age_group == "60 to 64", Attributable_resistance_mean * Efficacy * Coverage * (52+5)/(52*5),

              ifelse(Vaccination=="12 weeks" & Duration=="20" & Age_group == "Post Neonatal", Attributable_resistance_mean * Efficacy * Coverage * 39/48,
                                  ifelse(Vaccination=="12 weeks" & Duration=="20" & Age_group == "1 to 4", Attributable_resistance_mean * Efficacy * Coverage,
                                         ifelse(Vaccination=="12 weeks" & Duration=="20" & Age_group == "5 to 9", Attributable_resistance_mean * Efficacy * Coverage,
                                                ifelse(Vaccination=="12 weeks" & Duration=="20" & Age_group == "10 to 14", Attributable_resistance_mean * Efficacy * Coverage,
                                                       ifelse(Vaccination=="12 weeks" & Duration=="20" & Age_group == "15 to 19", Attributable_resistance_mean * Efficacy * Coverage,
                                                              ifelse(Vaccination=="12 weeks" & Duration=="20" & Age_group == "20 to 24", Attributable_resistance_mean * Efficacy * Coverage * 1/5 * 13/52,
                                                                     
             ifelse(Vaccination=="10 or above" & Duration=="10" & Age_group == "10 to 14", Attributable_resistance_mean * Efficacy * Coverage,
                    ifelse(Vaccination=="10 or above" & Duration=="10" & Age_group == "15 to 19", Attributable_resistance_mean * Efficacy * Coverage,
                                         
             ifelse(Vaccination=="60 or above" & Duration=="10" & Age_group == "60 to 64", Attributable_resistance_mean * Efficacy * Coverage,
                    ifelse(Vaccination=="60 or above" & Duration=="10" & Age_group == "65 to 69", Attributable_resistance_mean * Efficacy * Coverage,
                           
             ifelse(Vaccination=="60 or above" & Duration=="5" & Age_group == "60 to 64", Attributable_resistance_mean * Efficacy * Coverage,
                    
             ifelse(Vaccination=="60 or above" & Duration=="2" & Age_group == "60 to 64", Attributable_resistance_mean * Efficacy * Coverage * 2/5,
                    
             ifelse(Vaccination=="60 or above" & Duration=="1" & Age_group == "60 to 64", Attributable_resistance_mean * Efficacy * Coverage * 1/5,
 
                        Attributable_resistance_mean))))))))))))))))))))),
             va_Attributable_resistance_mean = Attributable_resistance_mean - v_Attributable_resistance_mean)
  
BurdenAverted <- function(pathogen){
  x <- burden_vaccine_profile[burden_vaccine_profile$Pathogen == pathogen,]
  
  x <- x %>%
    group_by(WHO_region) %>%
    summarise(averted_burden=sum(va_Attributable_resistance_mean), .groups = 'drop')

  x <- x[, 2]

  colnames(x) <- pathogen

  return(x)}

BurdenAverted_df <- lapply(pathogenlist, BurdenAverted)

BurdenAverted_df <- do.call(cbind.data.frame, BurdenAverted_df)

BurdenAverted_df$WHO_region <- WHO_region
# ------------------------------------------------------------------------------
# Table 2
Death_Averted <- BurdenAverted_df %>% 
  mutate(sum = rowSums(.[1:14]))

Death_Averted <- Death_Averted[-6,c("WHO_region", "sum")]

fwrite (x    = Death_Averted,
        file = "Death.Averted.csv")

# ------------------------------------------------------------------------------
# Figure 1
ggplot(Death_Averted) +
  aes(x = WHO_region, weight = sum) +
  geom_bar(fill = "#8DA5B6") +
  labs(x = "WHO region", y = "Vaccine Averted Deaths") +
  theme_classic()

# ------------------------------------------------------------------------------
# Figure 2

# ------------------------------------------------------------------------------
# Figure 3
