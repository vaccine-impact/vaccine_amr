# ------------------------------------------------------------------------------
# vaccine impact estimates using IHME AMR burden (Deaths)
# ------------------------------------------------------------------------------
# load libraries
library (readr)
library(readxl)
library (dplyr)
library (ggplot2)
library (reshape2)
library (tidyverse)

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



# ------------------------------------------------------------------------------
reference <- IHME_AMR_burden %>% count(Pathogen, Disease_presentation)
reference

pathogen <- IHME_AMR_burden %>% count(Pathogen)
pathogen

plot <- function(df, pathogen){
  
  y <- df[df$Pathogen == pathogen,]
  
  y <- y %>% 
    group_by(Pathogen, Age_group) %>%
    summarise(sum=sum(Attributable_resistance_mean))
    
  ggplot(y, aes(x=Age_group)) +
    geom_line(aes(y=sum, group=1, colour="vaccinated")) +
    ylab("Number of Vaccine Averted Deaths") + 
    xlab("Age") +
    ggtitle(pathogen) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5)) +
    expand_limits(y=0)
  
  }

plot(df = IHME_AMR_burden, pathogen = "Group A Streptococcus")






