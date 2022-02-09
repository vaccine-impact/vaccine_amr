# ------------------------------------------------------------------------------
# functions.R
#
# functions for analysis to estimate vaccine averted AMR health burden
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# create and clean up the data frame of AMR burden (deaths)

create_death_burden_table <- function(AMR_death_burden,
                                      death_burden_file){
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
          file = death_burden_file)
  
  return(AMR_death_burden)
} # end of function -- create_death_burden_table

# ------------------------------------------------------------------------------
# create and clean up the data frame of WHO vaccine profile

create_vaccine_profile_table <- function(vaccine_profile,
                                         vaccine_profile_file){
  
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
          file = vaccine_profile_file)
  
  return(vaccine_profile)
} # end of function -- create_vaccine_profile_table

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
          path = "figures",
          width = 15, 
          height = 6, 
          dpi = 600)
} # end of function -- create_death_by_pathogen_graph

# ------------------------------------------------------------------------------



# ------------------------------------------------------------------------------
# create combined table: disease burden + vaccine profile

create_combined_table <- function(death_burden_dt, 
                                  vaccine_profile_dt,
                                  attributable_burden_file,
                                  associated_burden_file){
  
# create combined table  
  combined_table <- left_join(death_burden_dt, vaccine_profile_dt, by=c("Pathogen" = "Pathogen"))
  combined_table <- data.table(combined_table)
  
# separate burden attributable to AMR and associated with AMR  
  
  attributable_burden <- combined_table[, -c(5:10)]
  
  attributable_burden <- attributable_burden %>% 
    rename("burden_lower_value" = "Attributable_resistance_lower",
           "burden_mean_value"  = "Attributable_resistance_mean",
           "burden_upper_value" = "Attributable_resistance_upper")
  
   fwrite (x    = attributable_burden,
           file = attributable_burden_file)
  
  associated_burden <- combined_table[, -c(5:7, 11:13)] 
  associated_burden <- associated_burden %>%
    rename("burden_lower_value" = "Associated_resistant_lower",
           "burden_mean_value"  = "Associated_resistant_mean",
           "burden_upper_value" = "Associated_resistant_upper")
  
  fwrite (x    = associated_burden,
          file = associated_burden_file)
  
  return(combined_table)
  
} # end of function -- create_combined_table

# ------------------------------------------------------------------------------
# log-normal distribution
# error in row 505
uncertainty_analysis_baseline <- function(psa, 
                                          data){ 
  
  burden_dt <- data.table(data)
  
  
  # + 0.5 for the value which is 0 (attributable to resistance)
  
  burden_dt$burden_lower_value <- burden_dt$burden_lower_value + 0.5
  
  burden_dt$burden_mean_value  <- burden_dt$burden_mean_value + 0.5
  
  burden_dt$burden_upper_value <- burden_dt$burden_upper_value + 0.5
  
  # ----------------------------------------------------------------------------
  # estimate mean log & sd log (attributable to resistance)
  burden_dt [, burden_sd_log := suppressMessages (suppressWarnings (get.lnorm.par (p           = c (0.025, 0.5, 0.975),
                                                                                   q           = c (burden_lower_value, burden_mean_value, burden_upper_value),
                                                                                   show.output = FALSE,
                                                                                   plot        = FALSE,
                                                                                   tol         = 0.0016) ["sdlog"]) ),
             by = .(WHO_region,Disease_presentation,Age_group,Pathogen)]
  
  burden_dt [, burden_mean_log := suppressMessages (suppressWarnings (get.lnorm.par (p           = c (0.025, 0.5, 0.975),
                                                                                     q           = c (burden_lower_value, burden_mean_value, burden_upper_value),
                                                                                     show.output = FALSE,
                                                                                     plot        = FALSE,
                                                                                     tol         = 0.0016) ["meanlog"]) ),
             by = .(WHO_region,Disease_presentation,Age_group,Pathogen)]
  
  
  # ----------------------------------------------------------------------------
  # create the table for psa
  
  burden_dt [, run_id := 0]
  
  vaccine_impact_psa <- data.table(burden_dt [0, ])
  
  for (i in 1:psa) {
    # copy data table
    dt <- copy (burden_dt)
    # set run id
    dt [, run_id := i]
    # add combined different scenarios
    vaccine_impact_psa <- rbindlist (list (vaccine_impact_psa, dt),
                                     use.names = TRUE)
  }
# ------------------------------------------------------------------------------
  # rlnorm will generate 'psa' number of values (attributable to resistance)
  
  vaccine_impact_psa [, burden_psa :=
                        rlnorm ( n       = psa,
                                 meanlog = mean (burden_mean_log),
                                 sdlog   = mean (burden_sd_log)),
                      by = .(WHO_region,Disease_presentation,Age_group,Pathogen) ]
  
  
  return(vaccine_impact_psa)
  
} # end of function -- uncertainty_analysis_baseline

# ------------------------------------------------------------------------------
# applying vaccine impact on vaccine target population

# applying vaccine target age group
estimate_vaccine_impact <- function(i, data){
  
  if(is.numeric(i) == TRUE) {
    vaccine_target_age <- data %>%
      filter(run_id == i)
  } else {
    vaccine_target_age <- data
  }
  
  vaccine_target_age <- vaccine_target_age %>%
    mutate(va_health_burden_age =
            ifelse(Vaccination=="4 week (effective at 6 week)" & Duration=="5" & Age_group == "PN", burden_psa * Efficacy * Coverage * 47/48,
            ifelse(Vaccination=="4 week (effective at 6 week)" & Duration=="5" & Age_group == "1 to 4", burden_psa * Efficacy * Coverage,
            ifelse(Vaccination=="4 week (effective at 6 week)" & Duration=="5" & Age_group == "5 to 9", burden_psa * Efficacy * Coverage * 1/5*5/52,
                                        
            ifelse(Vaccination=="4 week (effective at 6 week)" & Duration=="2" & Age_group == "PN" & Pathogen != "Streptococcus pneumoniae", burden_psa * Efficacy * Coverage * 47/48,
            ifelse(Vaccination=="4 week (effective at 6 week)" & Duration=="2" & Age_group == "1 to 4" & Pathogen != "Streptococcus pneumoniae", burden_psa * Efficacy * Coverage * 1/4*5/52,
                                                      
            ifelse(Vaccination=="4 week (effective at 6 week)" & Duration=="2" & Age_group == "PN" & Pathogen == "Streptococcus pneumoniae" & Disease_presentation == "LRI and thorax infections", burden_psa * 0.5 * Coverage * 47/48,
            ifelse(Vaccination=="4 week (effective at 6 week)" & Duration=="2" & Age_group == "1 to 4" & Pathogen == "Streptococcus pneumoniae" & Disease_presentation == "LRI and thorax infections", burden_psa * 0.5 * Coverage * 1/4*5/52,
                                                                    
            ifelse(Vaccination=="4 week (effective at 6 week)" & Duration=="2" & Age_group == "PN" & Pathogen == "Streptococcus pneumoniae" & Disease_presentation != "LRI and thorax infections", burden_psa * Efficacy * Coverage * 47/48,
            ifelse(Vaccination=="4 week (effective at 6 week)" & Duration=="2" & Age_group == "1 to 4" & Pathogen == "Streptococcus pneumoniae" & Disease_presentation != "LRI and thorax infections", burden_psa * Efficacy * Coverage * 1/4*5/52,
                                                                                  
                                                                                  
            ifelse(Vaccination=="4 week (effective at 6 week), 60 or above" & Duration=="2" & Age_group == "PN", burden_psa * Efficacy * Coverage * 47/48,
            ifelse(Vaccination=="4 week (effective at 6 week), 60 or above" & Duration=="2" & Age_group == "1 to 4", burden_psa * Efficacy * Coverage * 1/4*5/52,
            ifelse(Vaccination=="4 week (effective at 6 week), 60 or above" & Duration=="2" & Age_group == "60 to 64", burden_psa * Efficacy * Coverage * 2/5,
                                                                                                       
            ifelse(Vaccination=="12 weeks" & Duration=="20" & Age_group == "PN", burden_psa * Efficacy * Coverage * 39/48,
            ifelse(Vaccination=="12 weeks" & Duration=="20" & Age_group == "1 to 4", burden_psa * Efficacy * Coverage,
            ifelse(Vaccination=="12 weeks" & Duration=="20" & Age_group == "5 to 9", burden_psa * Efficacy * Coverage,
            ifelse(Vaccination=="12 weeks" & Duration=="20" & Age_group == "10 to 14", burden_psa * Efficacy * Coverage,
            ifelse(Vaccination=="12 weeks" & Duration=="20" & Age_group == "15 to 19", burden_psa * Efficacy * Coverage,
            ifelse(Vaccination=="12 weeks" & Duration=="20" & Age_group == "20 to 24", burden_psa * Efficacy * Coverage * 1/5 * 13/52,
                                                                                                                                                 
            ifelse(Vaccination=="10 or above" & Duration=="10" & Age_group == "10 to 14", burden_psa * Efficacy * Coverage,
            ifelse(Vaccination=="10 or above" & Duration=="10" & Age_group == "15 to 19", burden_psa * Efficacy * Coverage,
                                                                                                                                                               
            ifelse(Vaccination=="60 or above" & Duration=="10" & Age_group == "60 to 64", burden_psa * Efficacy * Coverage,
            ifelse(Vaccination=="60 or above" & Duration=="10" & Age_group == "65 to 69", burden_psa * Efficacy * Coverage,
                                                                                                                                                                             
            ifelse(Vaccination=="60 or above" & Duration=="5" & Age_group == "60 to 64", burden_psa * Efficacy * Coverage,
                                                                                                                                                                                    
            ifelse(Vaccination=="60 or above" & Duration=="2" & Age_group == "60 to 64", burden_psa * Efficacy * Coverage * 2/5,
                                                                                                                                                                                           
            ifelse(Vaccination=="60 or above" & Duration=="1" & Age_group == "60 to 64", burden_psa * Efficacy * Coverage * 1/5,
                   
                   0))))))))))))))))))))))))))

# applying vaccine target pathogen
  vaccine_impact <- vaccine_target_age %>%
    mutate(va_health_burden =
            ifelse(DiseasePresentation  == "All" | DiseasePresentation  == "All (50% for LRI and throx infections)", va_health_burden_age,
            ifelse(DiseasePresentation  == "BSI" & Disease_presentation == "BSI", va_health_burden_age, 
            ifelse(DiseasePresentation  == "Diarrhoea" & Disease_presentation == "Diarrhoea", va_health_burden_age,
            ifelse(DiseasePresentation  == "UTI" & Disease_presentation == "UTI", va_health_burden_age,
            ifelse(DiseasePresentation  == "BSI, LRI and thorax infections" & 
                  (Disease_presentation == "BSI" | Disease_presentation =="LRI and thorax infections"), va_health_burden_age,
                    0))))))
  
return(vaccine_impact)

 } # end of function -- estimate_vaccine_impact

# ------------------------------------------------------------------------------
# [table 2] Deaths and DALYs associated with and attributable to bacterial antimicrobial resistance
# globally and by WHO_region, 2019

aggregate_impact_by_region <- function(impact_by_region){

WHOregion <- unique(death_burden_dt$WHO_region)

death_averted_regional <- data.table("50%"=numeric(), "2.5%"=numeric(), "97.5%"=numeric())
 
impact_by_region_dt <- impact_by_region

for(i in WHOregion){
  dt <- impact_by_region_dt %>% filter(WHO_region == i)
  dt <- quantile(x = dt$averted_burden, probs = c (0.5, 0.025, 0.975))
  dt <- data.table(t(dt))
  death_averted_regional <- rbindlist (list (death_averted_regional, dt),
                                  use.names = FALSE)
}

death_averted_global <- impact_by_region_dt %>%
  group_by(run_id) %>%
  summarise(averted_burden=sum(averted_burden), .groups = 'drop')

death_averted_global <- quantile(x = death_averted_global$averted_burden, probs = c (0.5, 0.025, 0.975))
death_averted_global <- data.table(t(death_averted_global))

Death_Averted <- rbind(death_averted_global, death_averted_regional)

Death_Averted <- cbind(data.table(Counts=c("Global", WHOregion)), Death_Averted)

Death_Averted <- Death_Averted[,c("Counts", "2.5%", "50%", "97.5%")]

return(Death_Averted)} # end of function -- aggregate_impact_by_region

# ------------------------------------------------------------------------------
# Figure 1 -- Create vaccine averted burden by region

create_burden_averted_by_region_graph <- function(Attributable_burden_averted,
                                                  Associated_burden_averted){

Associated_burden_averted$Resistance <- "Associated with resistance" 

Attributable_burden_averted$Resistance <- "Attributable to resistance" 
  
  
burden_averted_by_region <- rbind(Associated_burden_averted, Attributable_burden_averted)

burden_averted_by_region <- burden_averted_by_region %>% rename("lower_value" = "2.5%",
                                                                "median_value"  = "50%",
                                                                "upper_value" = "97.5%")

burden_averted_by_region  <- burden_averted_by_region %>% filter(Counts != "Global")
                   
ggplot(burden_averted_by_region, aes(x = reorder(Counts, -median_value), y=median_value, fill=Resistance)) +
  geom_bar(stat = "identity", position="dodge") +
  scale_fill_manual(values = c("#D4E3FF","#054C70")) +
  labs(x = "WHO region", y = "Vaccine Avertable Deaths") +
  ylim(0,80000) +
  geom_errorbar(aes(ymin=lower_value, ymax=upper_value), width=0.25,
                size=0.5, position=position_dodge(0.9)) +
  theme_classic() +
  theme(legend.position = c(0.9, 0.9))

ggsave (filename = "Figure 1.png",
        path = "figures",
        width = 15, 
        height = 6, 
        dpi = 600)

} # end of function -- create_burden_averted_by_region_graph

# ------------------------------------------------------------------------------
# Figure 2 Create vaccine averted burden by disease presentation

aggregate_impact_by_dp <- function(impact_by_dp){
  
  DiseasePresentation <- unique(death_burden_dt$Disease_presentation)
    
  death_averted_dp <- data.table("50%"=numeric(), "2.5%"=numeric(), "97.5%"=numeric())
  
  impact_by_dp_dt <- impact_by_dp
  
  for(i in DiseasePresentation){
    dt <- impact_by_dp_dt %>% filter(Disease_presentation == i)
    dt <- quantile(x = dt$averted_burden, probs = c (0.5, 0.025, 0.975))
    dt <- data.table(t(dt))
    death_averted_dp <- rbindlist (list (death_averted_dp, dt),
                                     use.names = FALSE)
  }
  
  Death_Aveted <- cbind(data.table(Counts=DiseasePresentation, death_averted_dp))
  
  Death_Aveted <- Death_Aveted[,c("Counts", "2.5%", "50%", "97.5%")]
  
  return(Death_Aveted)} # end of function -- aggregate_impact_by_dp


# create Figure 2

create_burden_averted_by_dp_graph <- function(Attributable_burden_averted,
                                              Associated_burden_averted){
  
  Associated_burden_averted$Resistance <- "Associated with resistance" 
  
  Attributable_burden_averted$Resistance <- "Attributable to resistance" 
  
  
    burden_averted_by_dp <- rbind(Associated_burden_averted, Attributable_burden_averted)
  
    burden_averted_by_dp<-  burden_averted_by_dp %>% rename("lower_value" = "2.5%",
                                                            "median_value"  = "50%",
                                                            "upper_value" = "97.5%")
  
  ggplot(  burden_averted_by_dp, aes(x = reorder(Counts, -median_value), y=median_value, fill=Resistance)) +
    geom_bar(stat = "identity", position="dodge") +
    scale_fill_manual(values = c("#D4E3FF","#054C70")) +
    labs(x = "Infectious syndrome", y = "Vaccine Avertable Deaths") + 
    ylim(0,80000) +
    geom_errorbar(aes(ymin=lower_value, ymax=upper_value), width=0.25,
                  size=0.5, position=position_dodge(0.9)) +
    theme_classic(base_size=9) +
    theme(legend.position = c(0.9, 0.9))
  
  ggsave (filename = "Figure 2.png",
          path = "figures",
          width = 15, 
          height = 6, 
          dpi = 600)
  } # end of function -- create_burden_averted_by_dp_graph

# ------------------------------------------------------------------------------
# Figure 3 Create vaccine averted burden by pathogen

aggregate_impact_by_pathogen <- function(impact_by_pathogen){
 
  death_averted_pathogen <- data.table("50%"=numeric(), "2.5%"=numeric(), "97.5%"=numeric())
  
  impact_by_pathogen_dt <- impact_by_pathogen
  
  for(i in pathogenlist){
    dt <- impact_by_pathogen_dt %>% filter(Pathogen == i)
    dt <- quantile(x = dt$averted_burden, probs = c (0.5, 0.025, 0.975))
    dt <- data.table(t(dt))
    death_averted_pathogen <- rbindlist (list (death_averted_pathogen, dt),
                                         use.names = FALSE)
  }
  
  Death_Aveted <- cbind(data.table(Counts=pathogenlist, death_averted_pathogen))
  
  Death_Aveted <- Death_Aveted[,c("Counts", "2.5%", "50%", "97.5%")]
  
  return(Death_Aveted)} # end of function -- aggregate_impact_by_pathogen


# create Figure 3

create_burden_averted_by_pathogen_graph <- function(Attributable_burden_averted,
                                                    Associated_burden_averted){
  
  Associated_burden_averted$Resistance <- "Associated with resistance" 
  
  Attributable_burden_averted$Resistance <- "Attributable to resistance" 
  
  
  burden_averted_by_pathogen <- rbind(Associated_burden_averted, Attributable_burden_averted)
  
  burden_averted_by_pathogen<-  burden_averted_by_pathogen %>% rename("lower_value" = "2.5%",
                                                                      "median_value" = "50%",
                                                                      "upper_value" = "97.5%")
  
  ggplot(burden_averted_by_pathogen, aes(x = reorder(Counts, -median_value), y=median_value, fill=Resistance)) +
    geom_bar(stat = "identity", position="dodge") +
    scale_fill_manual(values = c("#D4E3FF","#054C70")) +
    labs(x = "Pathogen", y = "Vaccine Avertable Deaths") + 
    ylim(0,60000) +
    geom_errorbar(aes(ymin = lower_value, ymax = upper_value), width=0.25,
                  size=0.5, position=position_dodge(0.9)) +
    theme_classic(base_size=7.5) +
    theme(legend.position = c(0.9, 0.9))
  
  
  ggsave (filename = "Figure 3.png",
          path = "figures",
          width = 15, 
          height = 6, 
          dpi = 600)
  
} # end of function -- create_burden_averted_by_pathogen_graph
# ------------------------------------------------------------------------------

