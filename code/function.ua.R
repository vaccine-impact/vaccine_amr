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
# log-normal distribution
# error in row 505

vaccine_avertable_deaths_a <- left_join(death_burden_dt, vaccine_profile_dt, by=c("Pathogen" = "Pathogen"))
vaccine_avertable_deaths_a <- data.table(vaccine_avertable_deaths_a)

vaccine_avertable_deaths_a$Attributable_resistance_lower <- vaccine_avertable_deaths_a$Attributable_resistance_lower + 0.5
vaccine_avertable_deaths_a$Attributable_resistance_mean  <- vaccine_avertable_deaths_a$Attributable_resistance_mean + 0.5
vaccine_avertable_deaths_a$Attributable_resistance_upper <- vaccine_avertable_deaths_a$Attributable_resistance_upper + 0.5

vaccine_avertable_deaths_a [, attributable_sd_log := suppressMessages (suppressWarnings (get.lnorm.par (p = c (0.025, 0.5, 0.975),
                                                                                                        q = c (Attributable_resistance_lower, Attributable_resistance_mean, Attributable_resistance_upper),
                                                                                                        show.output = FALSE,
                                                                                                        plot        = FALSE) ["sdlog"]) ),
                            by = .(WHO_region,Disease_presentation,Age_group,Pathogen)]

vaccine_avertable_deaths_a [, attributable_mean_log := suppressMessages (suppressWarnings (get.lnorm.par (p = c (0.025, 0.5, 0.975),
                                                                                                          q = c (Attributable_resistance_lower, Attributable_resistance_mean, Attributable_resistance_upper),
                                                                                                          show.output = FALSE,
                                                                                                          plot        = FALSE) ["meanlog"]) ),
                            by = .(WHO_region,Disease_presentation,Age_group,Pathogen)]


# ------------------------------------------------------------------------------
# initialise psa data table

vaccine_avertable_deaths_a [, run_id := 0]

vaccine_impact_psa <- data.table(vaccine_avertable_deaths_a [0, ])

psa <- 100

for (i in 1:psa) {
  # copy data table
  dt <- copy (vaccine_avertable_deaths_a)
  # set run id
  dt [, run_id := i]
  # add combined vaccine impact estimates to benefit risk table
  vaccine_impact_psa <- rbindlist (list (vaccine_impact_psa, dt),
                                   use.names = TRUE)
}
# update estimates for deaths averted by vaccination
# note: mean_log is same value for a given ISO_code and Vaccine (similary for sd_log)
# which is why mean (mean_log) is done to get one value
# rlnorm will generate 'psa' number of values
# ----------------------------------------------------------------------------
vaccine_impact_psa [, random_attributable :=
                      rlnorm ( n      = psa,
                              meanlog = mean (attributable_mean_log),
                              sdlog   = mean (attributable_sd_log)),
                    by = .(WHO_region,Disease_presentation,Age_group,Pathogen) ]

# take into account vaccine target age group
estimate_vaccine_imapct_a <- function(i){

  vaccine_impact <- vaccine_impact_psa %>%
    filter(run_id == i)
  
  vaccine_impact <- vaccine_impact %>%
    mutate(v_health_burden 
           = ifelse(Vaccination=="4 week (effective at 6 week)" & Duration=="5" & Age_group == "Post Neonatal", random_attributable * Efficacy * Coverage * 47/48,
                    ifelse(Vaccination=="4 week (effective at 6 week)" & Duration=="5" & Age_group == "1 to 4", random_attributable * Efficacy * Coverage,
                    ifelse(Vaccination=="4 week (effective at 6 week)" & Duration=="5" & Age_group == "5 to 9", random_attributable * Efficacy * Coverage * 1/5*5/52,
                           
             ifelse(Vaccination=="4 week (effective at 6 week)" & Duration=="2" & Age_group == "Post Neonatal" & Pathogen != "Streptococcus pneumoniae", random_attributable * Efficacy * Coverage * 47/48,
                    ifelse(Vaccination=="4 week (effective at 6 week)" & Duration=="2" & Age_group == "1 to 4" & Pathogen != "Streptococcus pneumoniae", random_attributable * Efficacy * Coverage * 1/4*5/52,
 
             ifelse(Vaccination=="4 week (effective at 6 week)" & Duration=="2" & Age_group == "Post Neonatal" & Pathogen == "Streptococcus pneumoniae" & Disease_presentation == "LRI and thorax infections",  random_attributable * 0.5 * Coverage * 47/48,
                    ifelse(Vaccination=="4 week (effective at 6 week)" & Duration=="2" & Age_group == "1 to 4" & Pathogen == "Streptococcus pneumoniae" & Disease_presentation == "LRI and thorax infections", random_attributable * 0.5 * Coverage * 1/4*5/52,

             ifelse(Vaccination=="4 week (effective at 6 week)" & Duration=="2" & Age_group == "Post Neonatal" & Pathogen == "Streptococcus pneumoniae" & Disease_presentation != "LRI and thorax infections",  random_attributable * Efficacy * Coverage * 47/48,
                    ifelse(Vaccination=="4 week (effective at 6 week)" & Duration=="2" & Age_group == "1 to 4" & Pathogen == "Streptococcus pneumoniae" & Disease_presentation != "LRI and thorax infections", random_attributable * Efficacy * Coverage * 1/4*5/52,
                                         
                                                                    
             ifelse(Vaccination=="4 week (effective at 6 week), 60 or above" & Duration=="2" & Age_group == "Post Neonatal", random_attributable * Efficacy * Coverage * 47/48,
                    ifelse(Vaccination=="4 week (effective at 6 week), 60 or above" & Duration=="2" & Age_group == "1 to 4", random_attributable * Efficacy * Coverage * 1/4*5/52,
                    ifelse(Vaccination=="4 week (effective at 6 week), 60 or above" & Duration=="2" & Age_group == "60 to 64", random_attributable * Efficacy * Coverage * 2/5,

             ifelse(Vaccination=="12 weeks" & Duration=="20" & Age_group == "Post Neonatal", random_attributable * Efficacy * Coverage * 39/48,
                    ifelse(Vaccination=="12 weeks" & Duration=="20" & Age_group == "1 to 4", random_attributable * Efficacy * Coverage,
                    ifelse(Vaccination=="12 weeks" & Duration=="20" & Age_group == "5 to 9", random_attributable * Efficacy * Coverage,
                    ifelse(Vaccination=="12 weeks" & Duration=="20" & Age_group == "10 to 14", random_attributable * Efficacy * Coverage,
                    ifelse(Vaccination=="12 weeks" & Duration=="20" & Age_group == "15 to 19", random_attributable * Efficacy * Coverage,
                    ifelse(Vaccination=="12 weeks" & Duration=="20" & Age_group == "20 to 24", random_attributable * Efficacy * Coverage * 1/5 * 13/52,
                                                                     
             ifelse(Vaccination=="10 or above" & Duration=="10" & Age_group == "10 to 14", random_attributable * Efficacy * Coverage,
                    ifelse(Vaccination=="10 or above" & Duration=="10" & Age_group == "15 to 19", random_attributable * Efficacy * Coverage,
                                         
             ifelse(Vaccination=="60 or above" & Duration=="10" & Age_group == "60 to 64", random_attributable * Efficacy * Coverage,
                    ifelse(Vaccination=="60 or above" & Duration=="10" & Age_group == "65 to 69", random_attributable * Efficacy * Coverage,
                           
             ifelse(Vaccination=="60 or above" & Duration=="5" & Age_group == "60 to 64", random_attributable * Efficacy * Coverage,
                    
             ifelse(Vaccination=="60 or above" & Duration=="2" & Age_group == "60 to 64", random_attributable * Efficacy * Coverage * 2/5,
                    
             ifelse(Vaccination=="60 or above" & Duration=="1" & Age_group == "60 to 64", random_attributable * Efficacy * Coverage * 1/5,
 
                        random_attributable))))))))))))))))))))))))),
             va_health_burden = random_attributable - v_health_burden)

# take into account vaccine target pathogen
  vaccine_impact_b <- vaccine_impact [0, ]

  for(p in pathogenlist){
    dt <- vaccine_impact[vaccine_impact$Pathogen == p, ]
  
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
    
  vaccine_impact_b <- rbindlist (list (vaccine_impact_b, dt),
                                   use.names = TRUE)
  } 
  return(vaccine_impact_b)
  }

vaccine_impact_by_region <- data.table(WHO_region=NA, averted_burden=NA, run_id=NA)

for(i in 1:psa){
# copy data table
dt <- estimate_vaccine_imapct_a(i)
# set run id
dt <-  dt %>%
  group_by(WHO_region, run_id) %>%
  summarise(averted_burden=sum(va_health_burden), .groups = 'drop')
# add combined vaccine impact estimates to benefit risk table
vaccine_impact_by_region <- rbindlist (list (vaccine_impact_by_region  , dt),
                                       use.names = TRUE) 
}

# ------------------------------------------------------------------------------

a <- function(region){
  # copy data table
  dt <- vaccine_impact_by_region %>% filter(WHO_region == region)
  # add combined vaccine impact estimates to benefit risk table
  dt <- quantile (x = dt$averted_burden, probs = c (0.5, 0.025, 0.975))
  print(dt)
  }

a("Africa")

impact_by_region <- data.table(WHO_region=WHO_region)

b <- sapply(WHO_region, FUN=a)
b <- t(b) 
b <- data.table(b)
b <- b[, c("2.5%", "50%", "97.5%")]

Death_Averted <- cbind(impact_by_region,b)

fwrite (x    = Death_Averted,
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
  theme_classic(base_size=8.13)
