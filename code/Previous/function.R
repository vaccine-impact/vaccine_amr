# ------------------------------------------------------------------------------
# function.R
#
# functions for analysis to estimate vaccine averted AMR health burden
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# create and clean up the data frame of AMR burden (deaths)

create_death_burden_table <- function(AMR_death_burden){
  names(AMR_death_burden) <- c("WHO_region", "Disease_presentation","Age_group","Pathogen", 
                            "Associated_suscetible_mean","Associated_suscetible_lower",
                            "Associated_suscetible_upper","Associated_resistant_mean", 
                            "Associated_resistant_lower","Associated_resistant_upper",
                            "Attributable_resistance_mean","Attributable_resistance_lower", 
                            "Attributable_resistance_upper")


  AMR_death_burden <- AMR_death_burden[-c(1:2),]

  AMR_death_burden[,5:13] <- lapply(AMR_death_burden[,5:13], as.numeric)
  
  AMR_death_burden$Age_group <- factor(AMR_death_burden$Age_group, 
                                levels=c(as.vector(unlist(AMR_death_burden[1:23,"Age_group"]))), order=T)

  levels(AMR_death_burden$Age_group) <- c("EN", "LN", "PN",  "1 to 4", "5 to 9","10 to 14",
                                       "15 to 19", "20 to 24"  ,"25 to 29",  "30 to 34",
                                       "35 to 39", "40 to 44", "45 to 49",  "50 to 54" ,      
                                       "55 to 59", "60 to 64", "65 to 69", "70 to 74",
                                       "75 to 79", "80 to 84", "85 to 89", "90 to 94", "95 plus")
  
  AMR_death_burden <- AMR_death_burden %>%
    filter(AMR_death_burden$WHO_region != "unclassified")

  fwrite (x    = AMR_death_burden, 
          file = "AMR_death_burden.csv")
  
  return(AMR_death_burden)
  }

# ------------------------------------------------------------------------------
# create and clean up the data frame of WHO vaccine profile
  
create_vaccine_profile_table <- function(vaccine_profile){
    
  vaccine_profile <- vaccine_profile[,c(3, 4, 5, 6, 9, 10,12)]

  vaccine_profile <- vaccine_profile %>%
  rename("Efficacy" = "Efficacy (%)",
         "Coverage" = "Coverage in target group",
         "Duration" = "Duration of protection (year)", 
         "DiseasePresentation" = "Disease presentation")

  vaccine_profile$Efficacy <- vaccine_profile$Efficacy * 1/100
  vaccine_profile$Coverage <- vaccine_profile$Coverage * 1/100

  vaccine_profile <- vaccine_profile %>% filter(Selection == "Yes")

  fwrite (x    = vaccine_profile, 
        file = "Vaccine_profile.csv")

  return(vaccine_profile)
  }

# ------------------------------------------------------------------------------
# create graph of death trend by pathogen across all age groups

create_death_by_pathogen_graph <- function(pathogen){
  dt <- death_burden_dt[death_burden_dt$Pathogen == pathogen,]
  
  dt <- left_join(dt, vaccine_profile_dt, by=c("Pathogen" = "Pathogen"))
  
  dt <- dt %>% 
    group_by(Pathogen, Efficacy, Coverage, Age_group) %>%
    summarise(sum_Attributable_resistance_mean=sum(Attributable_resistance_mean), .groups = 'drop') 
  
  dt <- dt %>% 
    mutate(v_Attributable_resistance_mean = sum_Attributable_resistance_mean * Efficacy * Coverage,
           va_Attributable_resistance_mean = v_Attributable_resistance_mean-sum_Attributable_resistance_mean)
  
  ggplot(dt, aes(x=Age_group)) +
    geom_line(aes(y=sum_Attributable_resistance_mean, group=1, colour="non-vaccinated")) +
    ylab("Number of Vaccine Attributable Deaths") + 
    xlab("Age") +
    ggtitle(paste(pathogen,"(global)")) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5)) +
    expand_limits(y=0)
  
  # save plot
  ggsave (filename = paste(pathogen,"_global.png"),
          width = 15, 
          height = 6, 
          dpi = 600
          )
}

# ------------------------------------------------------------------------------



# ------------------------------------------------------------------------------
# vaccine avertable deaths attributable to resistance

# take into account vaccine target age group
estimate_vaccine_impact_a <- function(){
  
  vaccine_avertable_deaths_a <- left_join(death_burden_dt, vaccine_profile_dt, by=c("Pathogen" = "Pathogen"))
  
  vaccine_avertable_deaths_a <- vaccine_avertable_deaths_a %>%
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

# take into account vaccine target pathogen
  vaccine_avertable_deaths_b <- vaccine_avertable_deaths_a [0, ]

  for(i in pathogenlist){
    dt <- vaccine_avertable_deaths_a[vaccine_avertable_deaths_a$Pathogen == i, ]
  
    dt <- if(dt$DiseasePresentation[1]=="BSI") {
      dt[dt$Disease_presentation == "BSI",]
    } else if(dt$DiseasePresentation[1]=="Diarrhoea") {
      dt[dt$Disease_presentation == "Diarrhoea",]
    } else if(dt$DiseasePresentation[1]=="UTI") {
      dt[dt$Disease_presentation == "UTI",]
    } else if(dt$DiseasePresentation[1]=="BSI, LRI and thorax infections") {
      dt[dt$Disease_presentation == "BSI"|dt$Disease_presentation == "LRI and thorax infections",]
    } else {
      dt
    }
  vaccine_avertable_deaths_b <- rbindlist (list (vaccine_avertable_deaths_b, dt),
                                   use.names = TRUE)
  }

  vaccine_avertable_deaths_b <- vaccine_avertable_deaths_b %>%
    filter(vaccine_avertable_deaths_b$WHO_region != "unclassified")
  
  return(vaccine_avertable_deaths_b)
  }

# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# vaccine avertable associated resistance mean

estimate_vaccine_impact_b <- function() {
  
  vaccine_avertable_deaths_b <- left_join(death_burden_dt, vaccine_profile_dt, by=c("Pathogen" = "Pathogen"))

vaccine_avertable_deaths_b <- vaccine_avertable_deaths_b %>%
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
  x <- vaccine_avertable_deaths_b[vaccine_avertable_deaths_b$Pathogen == pathogen, ]
  
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

vaccine_avertable_deaths_b <- vaccine_avertable_deaths_b %>%
  filter(vaccine_avertable_deaths_b$WHO_region != "unclassified")

return(vaccine_avertable_deaths_b)}

# ------------------------------------------------------------------------------




# ------------------------------------------------------------------------------
# Table 2

AttributableResistance <-  BurdenAverted_attributable_dt%>%
  group_by(WHO_region) %>%
  summarise(averted_burden=sum(va_Attributable_resistance_mean), .groups = 'drop')

AssociatedResistant <- BurdenAverted_associated_dt %>%
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
 theme_classic()

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
  theme_classic(base_size=8.15)

# ------------------------------------------------------------------------------
# Figure 3
DeathAverted_Pathogen <- rbind(
  BurdenAverted_attributable_dt %>%
    mutate(Resistance = "Atributable to resistance") %>%
    group_by(Pathogen, Resistance) %>%
    summarise(averted_burden=sum(va_Attributable_resistance_mean),.groups = 'drop'),
  BurdenAverted_associated_dt %>%
    mutate(Resistance = "Associated with resistance") %>%
    group_by(Pathogen, Resistance) %>%
    summarise(averted_burden=sum(va_Associated_resistant_mean),.groups = 'drop'))

ggplot(DeathAverted_Pathogen, aes(x = reorder(Pathogen,-averted_burden), y=averted_burden, fill=Resistance)) +
  geom_bar(stat = "identity", position="dodge") +
  scale_fill_manual(values = c("#D4E3FF","#054C70")) +
  labs(x = "Pathogen", y = "Vaccine Avertable Deaths") +
  theme_classic(base_size=8.13)
