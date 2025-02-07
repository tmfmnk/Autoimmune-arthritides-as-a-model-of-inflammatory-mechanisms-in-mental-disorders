#Libraries

library(tidyverse)
library(survival)
library(broom)
library(patchwork)

#Import data

load(file = "path/AZV_rheuma_risk_of_psych_hosp/Data/data_main_analysis.RData")

#Absolute risk

imap_dfr(data_main_analysis,
         function(data_rheuma_dg, rheuma_cohort_names) {
           
           imap_dfr(data_rheuma_dg,
                    function(data_rheuma_and_psychiatric_dg, outcome_names) {
                      
                      data_rheuma_and_psychiatric_dg %>%
                        mutate(age_group = case_when(VEK <= 40 ~ "40 or less",
                                                     VEK >= 41 & VEK <= 59 ~ "41 to 59",
                                                     VEK >= 60 ~ "60 or more")) %>%
                        group_by(age_group, exposure) %>%
                        summarise(abs_risk = paste(formatC(sum(outcome == 1), big.mark = ","),
                                                   paste0("(",
                                                          formatC(round(mean(outcome) * 100, 2), format = "f", digits = 2),
                                                          ")"))) %>%
                        ungroup() %>%
                        mutate(cohort = rheuma_cohort_names,
                               outcome = outcome_names) 
                      
                    })
         }) %>%
  pivot_wider(names_from = c("cohort", "exposure", "age_group"),
              values_from = "abs_risk") %>%
  mutate(order = case_when(outcome == "Organic disorders" ~ 1,
                           outcome == "Alzheimer’s disease" ~ 2,
                           outcome == "Substance use disorders" ~ 3,
                           outcome == "Alcohol use disorders" ~ 4,
                           outcome == "Drug use disorders" ~ 5,
                           outcome == "Psychotic disorders" ~ 6,
                           outcome == "Schizophrenia" ~ 7,
                           outcome == "Mood disorders" ~ 8,
                           outcome == "Bipolar disorder" ~ 9,
                           outcome == "Depression" ~ 10,
                           outcome == "Moderate or severe depression" ~ 11,
                           outcome == "Anxiety disorders" ~ 12,
                           outcome == "Other anxiety disorders" ~ 13,
                           outcome == "Reaction to severe stress" ~ 14,
                           outcome == "Somatoform disorders" ~ 15,
                           outcome == "Other neurotic disorders" ~ 16,
                           outcome == "Behavioural syndromes" ~ 17)) %>%
  arrange(order) %>%
  select(-order) %>%
  write.csv(file = "path/AZV_rheuma_risk_of_psych_hosp/Results/Abs_risk_age_stratified.csv",
            row.names = FALSE)

#Fitting models
#Up to 40 years

mfull_up_to_40y <- imap_dfr(data_main_analysis,
                  function(data_rheuma_dg, rheuma_cohort_names) {
                    
                    imap_dfr(data_rheuma_dg,
                             function(data_rheuma_and_psychiatric_dg, outcome_names) {
                               
                               tidy(coxph(Surv(years_until_outcome_or_censoring, outcome) ~ exposure + VEK + POHL + month_discharge + year_discharge + strata(RODCIS2_exposed), 
                                          data = data_rheuma_and_psychiatric_dg,
                                          subset = VEK <= 40),
                                    conf.int = TRUE,
                                    exponentiate = TRUE) %>%
                                 mutate(cohort = rheuma_cohort_names,
                                        outcome = outcome_names) %>%
                                 filter(term == "exposureexposed")
                             })
                  })

#41 to 59 years

mfull_41_59y <- imap_dfr(data_main_analysis,
                            function(data_rheuma_dg, rheuma_cohort_names) {
                              
                              imap_dfr(data_rheuma_dg,
                                       function(data_rheuma_and_psychiatric_dg, outcome_names) {
                                         
                                         tidy(coxph(Surv(years_until_outcome_or_censoring, outcome) ~ exposure + VEK + POHL + month_discharge + year_discharge + strata(RODCIS2_exposed), 
                                                    data = data_rheuma_and_psychiatric_dg,
                                                    subset = VEK >= 41 & VEK <= 59),
                                              conf.int = TRUE,
                                              exponentiate = TRUE) %>%
                                           mutate(cohort = rheuma_cohort_names,
                                                  outcome = outcome_names) %>%
                                           filter(term == "exposureexposed")
                                       })
                            })

#60 years or more 

mfull_60y_or_more <- imap_dfr(data_main_analysis,
                         function(data_rheuma_dg, rheuma_cohort_names) {
                           
                           imap_dfr(data_rheuma_dg,
                                    function(data_rheuma_and_psychiatric_dg, outcome_names) {
                                      
                                      tidy(coxph(Surv(years_until_outcome_or_censoring, outcome) ~ exposure + VEK + POHL + month_discharge + year_discharge + strata(RODCIS2_exposed), 
                                                 data = data_rheuma_and_psychiatric_dg,
                                                 subset = VEK >= 60),
                                           conf.int = TRUE,
                                           exponentiate = TRUE) %>%
                                        mutate(cohort = rheuma_cohort_names,
                                               outcome = outcome_names) %>%
                                        filter(term == "exposureexposed")
                                    })
                         })
#Table

mfull_up_to_40y %>%
  mutate(age_group = "40 years or less") %>%
  bind_rows(mfull_41_59y %>%
              mutate(age_group = "41 to 59 years")) %>%
  bind_rows(mfull_60y_or_more %>%
              mutate(age_group = "60 years or more")) %>%
  mutate(age_group = factor(age_group,
                            levels = c("40 years or less",
                                       "41 to 59 years",
                                       "60 years or more"))) %>%
  mutate(cohort = factor(cohort,
                         levels = c("Rheumatoid arthritis or ankylosing spondylitis",
                                    "Rheumatoid arthritis",
                                    "Ankylosing spondylitis"))) %>%
  mutate(order = case_when(outcome == "Organic disorders" ~ 1,
                           outcome == "Alzheimer’s disease" ~ 2,
                           outcome == "Substance use disorders" ~ 3,
                           outcome == "Alcohol use disorders" ~ 4,
                           outcome == "Drug use disorders" ~ 5,
                           outcome == "Psychotic disorders" ~ 6,
                           outcome == "Schizophrenia" ~ 7,
                           outcome == "Mood disorders" ~ 8,
                           outcome == "Bipolar disorder" ~ 9,
                           outcome == "Depression" ~ 10,
                           outcome == "Moderate or severe depression" ~ 11,
                           outcome == "Anxiety disorders" ~ 12,
                           outcome == "Other anxiety disorders" ~ 13,
                           outcome == "Reaction to severe stress" ~ 14,
                           outcome == "Somatoform disorders" ~ 15,
                           outcome == "Other neurotic disorders" ~ 16,
                           outcome == "Behavioural syndromes" ~ 17)) %>%
  arrange(cohort, age_group, order) %>%
  transmute(cohort,
            age_group,
            outcome,
            order,
            estimate = paste(formatC(round(estimate, 2), format = "f", digits = 2),
                             paste0("(", 
                                    formatC(round(conf.low, 2), format = "f", digits = 2),
                                    "; ",
                                    formatC(round(conf.high, 2), format = "f", digits = 2),
                                    ")"))) %>%
  pivot_longer(-c(cohort, age_group, outcome, order)) %>%
  pivot_wider(names_from = c("cohort"),
              values_from = "value") %>%
  select(-name,
         -order) %>%
  write.csv(file = "path/AZV_rheuma_risk_of_psych_hosp/Results/HR_age_stratified_long.csv",
            row.names = FALSE)
