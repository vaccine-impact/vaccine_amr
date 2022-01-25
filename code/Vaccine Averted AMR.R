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
library (ggplot2)

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

## Death trend across all age groups: the outputs were used to decide the age of vaccination

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
# vaccine avertable attributable resistance mean
burden_vaccine_profile <- left_join(IHME_AMR_burden, Vaccine_profile, by=c("Pathogen" = "Pathogen"))
  
burden_vaccine_profile <- burden_vaccine_profile %>%
    mutate(v_Attributable_resistance_mean 
           = ifelse(Vaccination=="4 week (effective at 6 week)" & Duration=="5" & Age_group == "Post Neonatal", Attributable_resistance_mean * Efficacy * Coverage * 47/48,
                    ifelse(Vaccination=="4 week (effective at 6 week)" & Duration=="5" & Age_group == "1 to 4", Attributable_resistance_mean * Efficacy * Coverage,
                    ifelse(Vaccination=="4 week (effective at 6 week)" & Duration=="5" & Age_group == "5 to 9", Attributable_resistance_mean * Efficacy * Coverage * 1/5*5/52,
                           
             ifelse(Vaccination=="4 week (effective at 6 week)" & Duration=="2" & Age_group == "Post Neonatal" & Pathogen != "Streptococcus pneumoniae", Attributable_resistance_mean * Efficacy * Coverage * 47/48,
                    ifelse(Vaccination=="4 week (effective at 6 week)" & Duration=="2" & Age_group == "1 to 4" & Pathogen != "Streptococcus pneumoniae", Attributable_resistance_mean * Efficacy * Coverage * 1/4*5/52,
 
             ifelse(Vaccination=="4 week (effective at 6 week)" & Duration=="2" & Age_group == "Post Neonatal" & Pathogen == "Streptococcus pneumoniae" & Disease_presentation == "LRI and thorax infections",  Attributable_resistance_mean * 0.5 * Coverage * 47/48,
                    ifelse(Vaccination=="4 week (effective at 6 week)" & Duration=="2" & Age_group == "1 to 4" & Pathogen == "Streptococcus pneumoniae" & Disease_presentation == "LRI and thorax infections", Attributable_resistance_mean * 0.5 * Coverage * 1/4*5/52,

             ifelse(Vaccination=="4 week (effective at 6 week)" & Duration=="2" & Age_group == "Post Neonatal" & Pathogen == "Streptococcus pneumoniae" & Disease_presentation != "LRI and thorax infections",  Attributable_resistance_mean * Efficacy * Coverage * 47/48,
                    ifelse(Vaccination=="4 week (effective at 6 week)" & Duration=="2" & Age_group == "1 to 4" & Pathogen == "Streptococcus pneumoniae" & Disease_presentation != "LRI and thorax infections", Attributable_resistance_mean * Efficacy * Coverage * 1/4*5/52,
                                         
                                                                    
             ifelse(Vaccination=="4 week (effective at 6 week), 60 or above" & Duration=="2" & Age_group == "Post Neonatal", Attributable_resistance_mean * Efficacy * Coverage * 47/48,
                    ifelse(Vaccination=="4 week (effective at 6 week), 60 or above" & Duration=="2" & Age_group == "1 to 4", Attributable_resistance_mean * Efficacy * Coverage * 1/4*5/52,
                    ifelse(Vaccination=="4 week (effective at 6 week), 60 or above" & Duration=="2" & Age_group == "60 to 64", Attributable_resistance_mean * Efficacy * Coverage * 2/5,

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
 
                        Attributable_resistance_mean))))))))))))))))))))))))),
             va_Attributable_resistance_mean = Attributable_resistance_mean - v_Attributable_resistance_mean)

Disease_presentation <- function(pathogen){
  x <- burden_vaccine_profile[burden_vaccine_profile$Pathogen == pathogen, ]
  
  x <- if(x$DiseasePresentation[1]=="BSI") {
    x[x$Disease_presentation == "BSI",]
  } else if(x$DiseasePresentation[1]=="Diarrhoea") {
    x[x$Disease_presentation == "Diarrhoea",]
  } else if(x$DiseasePresentation[1]=="UTI") {
    x[x$Disease_presentation == "UTI",]
  } else if(x$DiseasePresentation[1]=="BSI, LRI and thorax infections") {
    x[x$Disease_presentation == "BSI"|x$Disease_presentation == "LRI and thorax infections",]
  } else {
    x
  }
  
  return(x)
}
  
BurdenAverted_df <- lapply(pathogenlist, Disease_presentation)

BurdenAverted_df <- do.call(rbind.data.frame, BurdenAverted_df)

BurdenAverted_df <- BurdenAverted_df %>%
  filter(BurdenAverted_df$WHO_region != "unclassified")
# ------------------------------------------------------------------------------
# vaccine avertable associated resistance mean
burden_vaccine_profile_a <- left_join(IHME_AMR_burden, Vaccine_profile, by=c("Pathogen" = "Pathogen"))

burden_vaccine_profile_a <- burden_vaccine_profile_a %>%
  mutate(v_Associated_resistant_mean 
         = ifelse(Vaccination=="4 week (effective at 6 week)" & Duration=="5" & Age_group == "Post Neonatal", Associated_resistant_mean * Efficacy * Coverage * 47/48,
                  ifelse(Vaccination=="4 week (effective at 6 week)" & Duration=="5" & Age_group == "1 to 4", Associated_resistant_mean * Efficacy * Coverage,
                  ifelse(Vaccination=="4 week (effective at 6 week)" & Duration=="5" & Age_group == "5 to 9", Associated_resistant_mean * Efficacy * Coverage * 1/5*5/52,
                                
           ifelse(Vaccination=="4 week (effective at 6 week)" & Duration=="2" & Age_group == "Post Neonatal" & Pathogen != "Streptococcus pneumoniae", Associated_resistant_mean * Efficacy * Coverage * 47/48,
                  ifelse(Vaccination=="4 week (effective at 6 week)" & Duration=="2" & Age_group == "1 to 4" & Pathogen != "Streptococcus pneumoniae", Associated_resistant_mean * Efficacy * Coverage * 1/4*5/52,
                                              
           ifelse(Vaccination=="4 week (effective at 6 week)" & Duration=="2" & Age_group == "Post Neonatal" & Pathogen == "Streptococcus pneumoniae" & Disease_presentation == "LRI and thorax infections",  Associated_resistant_mean * 0.5 * Coverage * 47/48,
                  ifelse(Vaccination=="4 week (effective at 6 week)" & Duration=="2" & Age_group == "1 to 4" & Pathogen == "Streptococcus pneumoniae" & Disease_presentation == "LRI and thorax infections", Associated_resistant_mean * 0.5 * Coverage * 1/4*5/52,
                                                            
           ifelse(Vaccination=="4 week (effective at 6 week)" & Duration=="2" & Age_group == "Post Neonatal" & Pathogen == "Streptococcus pneumoniae" & Disease_presentation != "LRI and thorax infections",  Associated_resistant_mean * Efficacy * Coverage * 47/48,
                  ifelse(Vaccination=="4 week (effective at 6 week)" & Duration=="2" & Age_group == "1 to 4" & Pathogen == "Streptococcus pneumoniae" & Disease_presentation != "LRI and thorax infections", Associated_resistant_mean * Efficacy * Coverage * 1/4*5/52,
                                                                          
                                                                          
           ifelse(Vaccination=="4 week (effective at 6 week), 60 or above" & Duration=="2" & Age_group == "Post Neonatal", Associated_resistant_mean * Efficacy * Coverage * 47/48,
                  ifelse(Vaccination=="4 week (effective at 6 week), 60 or above" & Duration=="2" & Age_group == "1 to 4", Associated_resistant_mean * Efficacy * Coverage * 1/4*5/52,
                  ifelse(Vaccination=="4 week (effective at 6 week), 60 or above" & Duration=="2" & Age_group == "60 to 64", Associated_resistant_mean * Efficacy * Coverage * 2/5,
                                                                                               
           ifelse(Vaccination=="12 weeks" & Duration=="20" & Age_group == "Post Neonatal", Associated_resistant_mean * Efficacy * Coverage * 39/48,
                  ifelse(Vaccination=="12 weeks" & Duration=="20" & Age_group == "1 to 4", Associated_resistant_mean * Efficacy * Coverage,
                  ifelse(Vaccination=="12 weeks" & Duration=="20" & Age_group == "5 to 9", Associated_resistant_mean * Efficacy * Coverage,
                  ifelse(Vaccination=="12 weeks" & Duration=="20" & Age_group == "10 to 14", Associated_resistant_mean * Efficacy * Coverage,
                  ifelse(Vaccination=="12 weeks" & Duration=="20" & Age_group == "15 to 19", Associated_resistant_mean * Efficacy * Coverage,
                  ifelse(Vaccination=="12 weeks" & Duration=="20" & Age_group == "20 to 24", Associated_resistant_mean * Efficacy * Coverage * 1/5 * 13/52,
                                                                                                                                         
           ifelse(Vaccination=="10 or above" & Duration=="10" & Age_group == "10 to 14", Associated_resistant_mean * Efficacy * Coverage,
                  ifelse(Vaccination=="10 or above" & Duration=="10" & Age_group == "15 to 19", Associated_resistant_mean * Efficacy * Coverage,
                                                                                                                                                       
           ifelse(Vaccination=="60 or above" & Duration=="10" & Age_group == "60 to 64", Associated_resistant_mean * Efficacy * Coverage,
                  ifelse(Vaccination=="60 or above" & Duration=="10" & Age_group == "65 to 69", Associated_resistant_mean * Efficacy * Coverage,
                                                                                                                                                                     
           ifelse(Vaccination=="60 or above" & Duration=="5" & Age_group == "60 to 64", Associated_resistant_mean * Efficacy * Coverage,
                                                                                                                                                                            
           ifelse(Vaccination=="60 or above" & Duration=="2" & Age_group == "60 to 64", Associated_resistant_mean * Efficacy * Coverage * 2/5,
                                                                                                                                                                                   
           ifelse(Vaccination=="60 or above" & Duration=="1" & Age_group == "60 to 64", Associated_resistant_mean * Efficacy * Coverage * 1/5,
                                                                                                                                                                                          
                      Associated_resistant_mean))))))))))))))))))))))))),
         va_Associated_resistant_mean = Associated_resistant_mean - v_Associated_resistant_mean)

Disease_presentation_a <- function(pathogen){
  x <- burden_vaccine_profile_a[burden_vaccine_profile_a$Pathogen == pathogen, ]
  
  x <- if(x$DiseasePresentation[1]=="BSI") {
    x[x$Disease_presentation == "BSI",]
  } else if(x$DiseasePresentation[1]=="Diarrhoea") {
    x[x$Disease_presentation == "Diarrhoea",]
  } else if(x$DiseasePresentation[1]=="UTI") {
    x[x$Disease_presentation == "UTI",]
  } else if(x$DiseasePresentation[1]=="BSI, LRI and thorax infections") {
    x[x$Disease_presentation == "BSI"|x$Disease_presentation == "LRI and thorax infections",]
  } else {
    x
  }
  
  return(x)
}

BurdenAverted_associated_df <- lapply(pathogenlist, Disease_presentation_a)

BurdenAverted_associated_df <- do.call(rbind.data.frame, BurdenAverted_associated_df)

BurdenAverted_associated_df <- BurdenAverted_associated_df %>%
  filter(BurdenAverted_associated_df$WHO_region != "unclassified")
# ------------------------------------------------------------------------------
# Uncertainty analysis: using log normal

?qnorm

# ------------------------------------------------------------------------------
# Table 2

AttributableResistance <-  BurdenAverted_df%>%
  group_by(WHO_region) %>%
  summarise(averted_burden=sum(va_Attributable_resistance_mean), .groups = 'drop')

AssociatedResistant <- BurdenAverted_associated_df %>%
  group_by(WHO_region) %>%
  summarise(averted_burden=sum(va_Associated_resistant_mean), .groups = 'drop')

Averted_Death <- left_join(AttributableResistance, AssociatedResistant, by=c("WHO_region" = "WHO_region"))

Averted_Death <- Averted_Death %>% 
  rename("Attributable to Resistance" = "averted_burden.x",
         "Associted to Resistance" = "averted_burden.y")

fwrite (x    = Averted_Death,
        file = "Averted_Death.csv")

# ------------------------------------------------------------------------------
# Figure 1
DeathAverted_Region <- rbind(AttributableResistance %>% mutate(Resistance = "Atributable to resistance"), 
      AssociatedResistant %>% mutate(Resistance = "Associated with resistance"))

ggplot(DeathAverted_Region, aes(x = WHO_region, y=averted_burden, fill=Resistance)) +
 geom_bar(stat = "identity", position="dodge") +
 scale_fill_manual(values = c("#D4E3FF","#054C70")) +
 labs(x = "WHO region", y = "Vaccine Avertable Deaths") + 
 theme_bw()

# ------------------------------------------------------------------------------
# Figure 2
DeathAverted_Syndrome <- rbind(
 BurdenAverted_df %>%
  mutate(Resistance = "Atributable to resistance") %>%
  group_by(Disease_presentation, Resistance) %>%
  summarise(averted_burden=sum(va_Attributable_resistance_mean),.groups = 'drop'), 
 BurdenAverted_associated_df %>%
  mutate(Resistance = "Associated with resistance") %>%
  group_by(Disease_presentation, Resistance) %>%
  summarise(averted_burden=sum(va_Associated_resistant_mean),.groups = 'drop'))

ggplot(DeathAverted_Syndrome, aes(x = reorder(Disease_presentation,-averted_burden), y=averted_burden, fill=Resistance)) +
  geom_bar(stat = "identity", position="dodge") +
  scale_fill_manual(values = c("#D4E3FF","#054C70")) +
  labs(x = "Infectious syndrome", y = "Vaccine Avertable Deaths") + 
  theme_bw(base_size=8.15)

# ------------------------------------------------------------------------------
# Figure 3
DeathAverted_Pathogen <- rbind(
  BurdenAverted_df %>%
    mutate(Resistance = "Atributable to resistance") %>%
    group_by(Pathogen, Resistance) %>%
    summarise(averted_burden=sum(va_Attributable_resistance_mean),.groups = 'drop'),
  BurdenAverted_associated_df %>%
    mutate(Resistance = "Associated with resistance") %>%
    group_by(Pathogen, Resistance) %>%
    summarise(averted_burden=sum(va_Associated_resistant_mean),.groups = 'drop'))

ggplot(DeathAverted_Pathogen, aes(x = reorder(Pathogen,-averted_burden), y=averted_burden, fill=Resistance)) +
  geom_bar(stat = "identity", position="dodge") +
  scale_fill_manual(values = c("#D4E3FF","#054C70")) +
  labs(x = "Pathogen", y = "Vaccine Avertable Deaths") +
  theme_bw(base_size=8.13)
