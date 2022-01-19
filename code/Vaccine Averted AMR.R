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

names(IHME_AMR_burden) <- c("WHO_region", "Disease_presentation","Age_group","Pathogen", 
                            "Associated_suscetible_mean","Associated_suscetible_lower",
                            "Associated_suscetible_upper","Associated_resistant_mean", 
                            "Associated_resistant_lower","Associated_resistant_upper",
                            "Attributable_resistance_mean","Attributable_resistance_lower", 
                            "Attributable_resistance_upper")

# Data cleaning
IHME_AMR_burden <- IHME_AMR_burden[-c(1:2),]

IHME_AMR_burden[,5:13] <- lapply(IHME_AMR_burden[,5:13], as.numeric)

IHME_AMR_burden$Age_group <- factor(IHME_AMR_burden$Age_group, 
                                levels=c(as.vector(unlist(IHME_AMR_burden[1:23,"Age_group"]))), order=T)
# ------------------------------------------------------------------------------
# Import WHO vaccine profile
Vaccine_profile <- read_excel("Vaccine profile for IHME burden.xlsx", 
                              sheet = "Vaccine profile assumptions")

Vaccine_profile <- Vaccine_profile[,c(3, 4, 5, 6, 9, 10)]

Vaccine_profile <- Vaccine_profile %>%
  rename("Efficacy" = "Efficacy (%)",
         "Coverage" = "Coverage in target group",
         "Duration" = "Duration of protection (year)", 
         "DiseasePresentation" = "Disease presentation")

Vaccine_profile$Efficacy <- Vaccine_profile$Efficacy * 1/100
Vaccine_profile$Coverage <- Vaccine_profile$Coverage * 1/100

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
##"Group A Streptococcus"

Group_A_Streptococcus <- IHME_AMR_burden[IHME_AMR_burden$Pathogen == "Group A Streptococcus",]

Group_A_Streptococcus <- left_join(Group_A_Streptococcus, Vaccine_profile, by=c("Pathogen" = "Pathogen"))

Group_A_Streptococcus <- Group_A_Streptococcus %>% 
  group_by(WHO_region, Pathogen, Efficacy, Coverage, Age_group) %>%
  summarise(sum_Attributable_resistance_mean=sum(Attributable_resistance_mean), .groups = 'drop') 

Group_A_Streptococcus <- Group_A_Streptococcus %>% 
  mutate(v_Attributable_resistance_mean = ifelse(Age_group >= "Post Neonatal" & Age_group < "5 to 9",
                                                 sum_Attributable_resistance_mean * Efficacy * Coverage,
                                                 sum_Attributable_resistance_mean),
        va_Attributable_resistance_mean = v_Attributable_resistance_mean-sum_Attributable_resistance_mean)

burdenplot <- function(region) {
  
y  <- Group_A_Streptococcus[Group_A_Streptococcus$WHO_region == region,] 

ggplot(y, aes(x=Age_group)) +
  geom_line(aes(y=sum_Attributable_resistance_mean, group=1, colour="non-vaccinated")) +
  geom_line(aes(y=v_Attributable_resistance_mean, group=1, colour="vaccinated")) +
  ylab("Number of Vaccine Attributable Deaths") + 
  xlab("Age") +
  ggtitle(paste("Group A Streptococcus","(", region,")")) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  expand_limits(y=0)
}

lapply(WHO_region, burdenplot)

Vaccine_averted_burden <- Group_A_Streptococcus %>%
  group_by(WHO_region) %>%
  summarise(averted_burden=sum(sum_Attributable_resistance_mean))

# ------------------------------------------------------------------------------
