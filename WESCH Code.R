##Packages
library(readxl)
library(tidyverse)
library(janitor)
library(brms)

##Data and seed
wesch_data <- read_xlsx("WESCH Data.xlsx") %>% 
  clean_names() %>% 
  mutate(condition_string = as.factor(condition_string),
         age_group_string = as.factor(age_group_string),
         location_string = as.factor(location_string),
         gender_string = as.factor(gender_string))

wesch_data$age_group_string <- relevel(wesch_data$age_group_string, ref = "Younger")
wesch_data$condition_string <- relevel(wesch_data$condition_string, ref = "You")

set.seed(31)

#Priors
prior_tendencies <- c(
  set_prior("normal(0,1.25)", class = "Intercept"), #median = exp(0) = 1 baseline odds; -2sd = exp(-2.5) = 0.08 baseline odds; +2sd = 12.2 baseline odds
  set_prior("normal(0,1)", class = "b")) #median = exp(0) = 1 odds ratio; -2sd = exp(-2) = 0.13 odds ratio; +2sd = 7.39 odds ratio

qqplot_commitment <- wesch_data %>% 
  filter(commitment_latency < 61)

qqplot_helping <- wesch_data %>% 
  filter(helping_latency < 98)

ggplot(qqplot_commitment, aes(sample = commitment_latency)) +
  stat_qq(distribution = stats::qlnorm) +
  stat_qq_line() +
  ggtitle("Commitment latency, normal quantiles") #looks approximately lognormal

ggplot(qqplot_helping, aes(sample = helping_latency)) +
  stat_qq(distribution = stats::qnorm) +
  stat_qq_line() +
  ggtitle("Helping latency, normal quantiles") #looks approximately lognormal

prior_latencies <- c(
  set_prior("normal(3,0.5)", class = "Intercept"), #median = exp(3) = 20 seconds; -2sd = exp(3)*exp(1) = 55 seconds; +2sd = exp(3)*exp(-1) = 7 seconds
  set_prior("normal(0,0.5)", class = "b"))

prior_survival <- c(
  set_prior("normal(0,1)", class = "b"))
exp(0) # = median hazard ratio of non-intercept predictors; a hazard ratio of 1 = no change between levels of predictor
exp(1) # = 1SD of hazard ratio of non-intercept predictors; coefficient estimates must be exponentiated to get the hazard ratio

##Commitment
#Data
commit_tendency <- wesch_data %>% 
  filter(approach_n_y < 9)

commit_tendency_we <- wesch_data %>% 
  filter(condition_string != "You")

commit_tendency_you <- wesch_data %>% 
  filter(condition_string != "We")

commit_tendency_3yo <- wesch_data %>% 
  filter(age_years < 3.5)

commit_tendency_4yo <- wesch_data %>% 
  filter(age_years > 3.5)

commit_latency <- wesch_data %>% 
  mutate(censored = ifelse(commitment_latency > 60, "right", "none"),
         commitment_latency = ifelse(censored == "right", 60, commitment_latency)) %>% 
  select(censored, commitment_latency, age_stndrd, age_group_string, condition_string, location_string, first, gender_string) %>% 
  relocate(censored, .after = commitment_latency) %>% 
  filter(censored == "none")

commit_latency_surv <- wesch_data %>% 
  mutate(censored = ifelse(commitment_latency > 60, "right", "none"),
         commitment_latency = ifelse(censored == "right", 60, commitment_latency)) %>% 
  select(censored, commitment_latency, age_stndrd, age_group_string, condition_string, location_string, first, gender_string) %>% 
  relocate(censored, .after = commitment_latency)

leave_taking <- wesch_data %>%
  filter(leave_taking_type < 4)

commit_leavetake_we <- leave_taking %>% 
  filter(condition_string != "you")

commit_leavetake_you <- leave_taking %>% 
  filter(condition_string != "we")

commit_leavetake_3yo <- leave_taking %>% 
  filter(age_years < 3.5)

commit_leavetake_4yo <- leave_taking %>% 
  filter(age_years > 3.5)

#Tendency to abandon
commit_glmm1 <- brm(approach_n_y ~ age_group_string + location_string + first + gender_string +
                      (1|location_string),
                    data = commit_tendency,
                    family = bernoulli,
                    prior = prior_tendencies,
                    iter = 10000,
                    warmup = 2000,
                    control = list(adapt_delta = 0.9999, max_treedepth = 50),
                    save_pars = save_pars(all = TRUE),
                    cores = 4)

commit_glmm2 <- brm(approach_n_y ~ condition_string + age_group_string + location_string + first + gender_string+
                      (1|location_string),
                    data = commit_tendency,
                    family = bernoulli,
                    prior = prior_tendencies,
                    iter = 10000,
                    warmup = 2000,
                    control = list(adapt_delta = 0.9999, max_treedepth = 50),
                    save_pars = save_pars(all = TRUE),
                    cores = 4)

commit_glmm3 <- brm(approach_n_y ~ condition_string * age_group_string + location_string + first + gender_string +
                      (1|location_string),
                    data = commit_tendency,
                    family = bernoulli,
                    prior = prior_tendencies,
                    iter = 10000,
                    warmup = 2000,
                    control = list(adapt_delta = 0.9999, max_treedepth = 50),
                    save_pars = save_pars(all = TRUE),
                    cores = 4)

plot(commit_glmm1)
plot(commit_glmm2)
plot(commit_glmm3)

pp_check(commit_glmm3, nsamples = 100)

post_prob_commit123 <- post_prob(commit_glmm1, commit_glmm2, commit_glmm3)

BF_commit_abandon <- 0.5349811/0.2489203

#Abandoning tendency: Looking at age within we-framing
commit_glmm4 <- brm(approach_n_y ~ age_group_string +
                      (1|location_string),
                    data = commit_tendency_we,
                    family = bernoulli,
                    prior = prior_tendencies,
                    iter = 10000,
                    warmup = 2000,
                    control = list(adapt_delta = 0.9999, max_treedepth = 50),
                    save_pars = save_pars(all = TRUE),
                    cores = 4)

plot(commit_glmm4)
pp_check(commit_glmm4, nsamples = 100)

#Abandoning tendency: Looking at age within you-framing
commit_glmm5 <- brm(approach_n_y ~ age_group_string +
                      (1|location_string),
                    data = commit_tendency_you,
                    family = bernoulli,
                    prior = prior_tendencies,
                    iter = 10000,
                    warmup = 2000,
                    control = list(adapt_delta = 0.9999, max_treedepth = 50),
                    save_pars = save_pars(all = TRUE),
                    cores = 4)

plot(commit_glmm5)
pp_check(commit_glmm5, nsamples = 100)

#Abandoning tendency: Looking at condition within age [2.5,3.5]
commit_glmm6 <- brm(approach_n_y ~ condition_string +
                      (1|location_string),
                    data = commit_tendency_3yo,
                    family = bernoulli,
                    prior = prior_tendencies,
                    iter = 10000,
                    warmup = 2000,
                    control = list(adapt_delta = 0.9999, max_treedepth = 50),
                    save_pars = save_pars(all = TRUE),
                    cores = 4)

plot(commit_glmm6)
pp_check(commit_glmm6, nsamples = 100) #giving some strange posterior predictions

#Abandoning tendency: Looking at condition within age [3.5,4.5]
commit_glmm7 <- brm(approach_n_y ~ condition_string +
                      (1|location_string),
                    data = commit_tendency_4yo,
                    family = bernoulli,
                    prior = prior_tendencies,
                    iter = 10000,
                    warmup = 2000,
                    control = list(adapt_delta = 0.9999, max_treedepth = 50),
                    save_pars = save_pars(all = TRUE),
                    cores = 4)

plot(commit_glmm7)
pp_check(commit_glmm7, nsamples = 100)

#Latency to abandon
commit_glmm8 <- brm(commitment_latency ~ age_group_string + location_string + first + gender_string +
                      (1|location_string),
                    data = commit_latency,
                    family = lognormal,
                    prior = prior_latencies,
                    iter = 10000,
                    warmup = 2000,
                    control = list(adapt_delta = 0.9999, max_treedepth = 60),
                    save_pars = save_pars(all = TRUE),
                    cores = 4)

commit_glmm9 <- brm(commitment_latency ~ age_group_string + condition_string + location_string + first + gender_string +
                      (1|location_string),
                    data = commit_latency,
                    family = lognormal,
                    prior = prior_latencies,
                    iter = 10000,
                    warmup = 2000,
                    control = list(adapt_delta = 0.9999, max_treedepth = 60),
                    save_pars = save_pars(all = TRUE),
                    cores = 4)

commit_glmm10 <- brm(commitment_latency ~ age_group_string * condition_string + location_string + first + gender_string +
                       (1|location_string),
                     data = commit_latency,
                     family = lognormal,
                     prior = prior_latencies,
                     iter = 10000,
                     warmup = 2000,
                     control = list(adapt_delta = 0.9999, max_treedepth = 60),
                     save_pars = save_pars(all = TRUE),
                     cores = 4)

plot(commit_glmm8)
plot(commit_glmm9)
plot(commit_glmm10)

loo(commit_glmm8)
loo(commit_glmm9)
loo(commit_glmm10)

pp_check(commit_glmm10, nsamples = 100)

post_prob_commit8910 <- post_prob(commit_glmm8, commit_glmm9, commit_glmm10)

BF_commit_latency <- 0.4762492/0.2789535

commit_glmm8_surv <- brm(commitment_latency ~ age_group_string + location_string + first + gender_string,
                         data = commit_latency_surv,
                         family = cox,
                         prior = prior_survival,
                         iter = 10000,
                         warmup = 2000,
                         control = list(adapt_delta = 0.9999, max_treedepth = 60),
                         save_pars = save_pars(all = TRUE),
                         cores = 4)

commit_glmm9_surv <- brm(commitment_latency ~ age_group_string + condition_string + location_string + first + gender_string,   
                         data = commit_latency_surv,
                         family = cox,
                         prior = prior_survival,
                         iter = 10000,
                         warmup = 2000,
                         control = list(adapt_delta = 0.9999, max_treedepth = 60),
                         save_pars = save_pars(all = TRUE),
                         cores = 4)

commit_glmm10_surv <- brm(commitment_latency ~ age_group_string * condition_string + location_string + first + gender_string,
                          data = commit_latency_surv,
                          family = cox,
                          prior = prior_survival,
                          iter = 10000,
                          warmup = 2000,
                          control = list(adapt_delta = 0.9999, max_treedepth = 60),
                          save_pars = save_pars(all = TRUE),
                          cores = 4)

plot(commit_glmm8_surv)
plot(commit_glmm9_surv)
plot(commit_glmm10_surv)

post_prob_commit8910_surv <- post_prob(commit_glmm8_surv, commit_glmm9_surv, commit_glmm10_surv)

BF_commit_surv <- 0.76339170/0.15687151

#Leave-taking
commit_glmm11 <- brm(leave_taking_type ~ age_group_string * condition_string + location_string + first + gender_string +
                       (1|location_string),
                     data = leave_taking,
                     family = categorical,
                     iter = 10000,
                     warmup = 2000,
                     control = list(adapt_delta = 0.9999, max_treedepth = 50),
                     cores = 4)

plot(commit_glmm11) ##trace plots look weird, likely too little data coming from likelihood term

commit_glmm12 <- brm(leave_taking_type_recoded ~ age_group_string + location_string + first + gender_string +
                       (1|location_string),
                     data = leave_taking,
                     family = bernoulli,
                     prior = prior_tendencies,
                     iter = 10000,
                     warmup = 2000,
                     control = list(adapt_delta = 0.9999, max_treedepth = 50),
                     save_pars = save_pars(all = TRUE),
                     cores = 4)

commit_glmm13 <- brm(leave_taking_type_recoded ~ age_group_string + condition_string + location_string + first + gender_string +
                       (1|location_string),
                     data = leave_taking,
                     family = bernoulli,
                     prior = prior_tendencies,
                     iter = 10000,
                     warmup = 2000,
                     control = list(adapt_delta = 0.9999, max_treedepth = 50),
                     save_pars = save_pars(all = TRUE),
                     cores = 4)

commit_glmm14 <- brm(leave_taking_type_recoded ~ age_group_string * condition_string + location_string + first + gender_string +
                       (1|location_string),
                     data = leave_taking,
                     family = bernoulli,
                     prior = prior_tendencies,
                     iter = 10000,
                     warmup = 2000,
                     control = list(adapt_delta = 0.9999, max_treedepth = 50),
                     save_pars = save_pars(all = TRUE),
                     cores = 4)

plot(commit_glmm12)
plot(commit_glmm13)
plot(commit_glmm14)

pp_check(commit_glmm12, nsamples = 500)

post_prob_commit121314 <- post_prob(commit_glmm12, commit_glmm13, commit_glmm14)

BF_commit_takeleave <- 0.4305622/0.4180324

hypothesis(commit_glmm14, "condition_stringWe > 0")

posterior_samples_takeleave <- posterior_samples(commit_glmm14)
HPD_takeleave_diff <- mean(posterior_samples_takeleave$b_condition_stringWe < 0) - 
  mean(posterior_samples_takeleave$b_condition_stringWe < -0.08)
HPD_takeleave_ratio <- (HPD_takeleave_diff/.95)

table(leave_taking$first, leave_taking$leave_taking_type_recoded)

leave_taking_noorder3 <- leave_taking %>% 
  filter(first != 5)

commit_glmm14_noorder3 <- brm(leave_taking_type_recoded ~ age_group_string * condition_string + location_string + first + gender_string +
                       (1|location_string),
                       data = leave_taking_noorder3,
                       family = bernoulli,
                       prior = prior_tendencies,
                       iter = 10000,
                       warmup = 2000,
                       control = list(adapt_delta = 0.9999, max_treedepth = 50),
                       save_pars = save_pars(all = TRUE),
                       cores = 4)

#Leave-taking tendency: Looking at age within we-framing
commit_glmm15 <- brm(leave_taking_type_recoded ~ age_group_string +
                       (1|location_string),
                    data = commit_leavetake_we,
                    family = bernoulli,
                    prior = prior_tendencies,
                    iter = 10000,
                    warmup = 2000,
                    control = list(adapt_delta = 0.9999, max_treedepth = 50),
                    save_pars = save_pars(all = TRUE),
                    cores = 4)

plot(commit_glmm15)
pp_check(commit_glmm15, nsamples = 100)

#Leave-taking tendency: Looking at age within you-framing
commit_glmm16 <- brm(leave_taking_type_recoded ~ age_group_string +
                       (1|location_string),
                    data = commit_leavetake_you,
                    family = bernoulli,
                    prior = prior_tendencies,
                    iter = 10000,
                    warmup = 2000,
                    control = list(adapt_delta = 0.9999, max_treedepth = 50),
                    save_pars = save_pars(all = TRUE),
                    cores = 4)

plot(commit_glmm16)
pp_check(commit_glmm16, nsamples = 100)

#Leave-taking tendency: Looking at condition within age [2.5,3.5]
commit_glmm17 <- brm(leave_taking_type_recoded ~ condition_string +
                       (1|location_string),
                    data = commit_leavetake_3yo,
                    family = bernoulli,
                    prior = prior_tendencies,
                    iter = 10000,
                    warmup = 2000,
                    control = list(adapt_delta = 0.9999, max_treedepth = 50),
                    save_pars = save_pars(all = TRUE),
                    cores = 4)

plot(commit_glmm17)
pp_check(commit_glmm17, nsamples = 100)

hypothesis(commit_glmm17, "condition_stringWe > 0")

#Leave-taking tendency: Looking at condition within age [3.5,4.5]
commit_glmm18 <- brm(leave_taking_type_recoded ~ condition_string +
                       (1|location_string),
                    data = commit_leavetake_4yo,
                    family = bernoulli,
                    prior = prior_tendencies,
                    iter = 10000,
                    warmup = 2000,
                    control = list(adapt_delta = 0.9999, max_treedepth = 50),
                    save_pars = save_pars(all = TRUE),
                    cores = 4)

plot(commit_glmm18)
pp_check(commit_glmm18, nsamples = 100)

##Sharing
#Data
sharing <- wesch_data %>% 
  filter(number_of_erasers_shared < 7)

#Resource distribution
sharing_glmm1 <- brm(number_of_erasers_shared ~ age_group_string * condition_string + location_string + first + gender_string,
                     data = sharing,
                     family = categorical,
                     iter = 10000,
                     warmup = 2000,
                     control = list(adapt_delta = 0.9999, max_treedepth = 50),
                     cores = 4)

plot(sharing_glmm1) ##trace plots look weird, likely too little data coming from likelihood term


sharing_glmm2 <- brm(number_of_erasers_shared_recoded ~ age_group_string + location_string + first + gender_string +
                       (1|location_string),
                     data = sharing,
                     family = bernoulli,
                     prior = prior_tendencies,
                     iter = 10000,
                     warmup = 2000,
                     control = list(adapt_delta = 0.9999, max_treedepth = 50),
                     save_pars = save_pars(all = TRUE),
                     cores = 4)

sharing_glmm3 <- brm(number_of_erasers_shared_recoded ~ age_group_string + condition_string + location_string + first + gender_string +
                       (1|location_string),
                     data = sharing,
                     family = bernoulli,
                     prior = prior_tendencies,
                     iter = 10000,
                     warmup = 2000,
                     control = list(adapt_delta = 0.9999, max_treedepth = 50),
                     save_pars = save_pars(all = TRUE),
                     cores = 4)

sharing_glmm4 <- brm(number_of_erasers_shared_recoded ~ age_group_string * condition_string + location_string + first + gender_string +
                       (1|location_string),
                     data = sharing,
                     family = bernoulli,
                     prior = prior_tendencies,
                     iter = 10000,
                     warmup = 2000,
                     control = list(adapt_delta = 0.9999, max_treedepth = 50),
                     save_pars = save_pars(all = TRUE),
                     cores = 4)

plot(sharing_glmm2)
plot(sharing_glmm3)
plot(sharing_glmm4)

pp_check(sharing_glmm2, nsamples = 100)

post_prob_sharing234 <- post_prob(sharing_glmm2, sharing_glmm3, sharing_glmm4)

BF_share <- 0.4007365/0.3440198

#Looking at proportion of equal sharers
sharing_equal <- sharing %>% 
  mutate(equal_share = ifelse(number_of_erasers_shared == 3, 1, 0))

table(sharing_equal$equal_share, sharing_equal$condition_string)

sharing_generous <- sharing %>% 
  mutate(generous_share = ifelse(number_of_erasers_shared > 3, 1, 0))

table(sharing_generous$generous_share, sharing_generous$condition_string)

##Helping
#Data
helping <- wesch_data %>% 
  filter(help_n_y < 5) %>% 
  mutate(zero_inflated_help = ifelse(help_n_y == 0, 1, 0)) %>% 
  relocate(zero_inflated_help, .after = help_n_y)

helping_tendency_we <- wesch_data %>% 
  filter(condition_string != "You",
         help_n_y < 6)

helping_tendency_you <- wesch_data %>% 
  filter(condition_string != "We",
         help_n_y < 6)

helping_tendency_3yo <- wesch_data %>% 
  filter(age_years < 3.5,
         help_n_y < 6)

helping_tendency_4yo <- wesch_data %>% 
  filter(age_years > 3.5,
         help_n_y < 6)

help_latency <- helping %>% 
  mutate(censored = ifelse(helping_latency > 60, "right", "none"),
         helping_latency = ifelse(censored == "right", 60, helping_latency)) %>% 
  select(censored, age_group_string, helping_latency, age_stndrd, condition_string, location_string, first, gender_string) %>% 
  relocate(censored, .after = helping_latency) %>% 
  filter(helping_latency != 60)
46/47
help_latency_surv <- helping %>% 
  mutate(censored = ifelse(helping_latency > 60, "right", "none"),
         helping_latency = ifelse(censored == "right", 60, helping_latency)) %>% 
  select(censored, age_group_string, helping_latency, age_stndrd, condition_string, location_string, first, gender_string) %>% 
  relocate(censored, .after = helping_latency)

#Tendency to help
helping_glmm1 <- brm(help_n_y ~ age_group_string + location_string + first + gender_string +
                       (1|location_string),
                     data = helping,
                     family = bernoulli,
                     prior = prior_tendencies,
                     iter = 10000,
                     warmup = 2000,
                     control = list(adapt_delta = 0.9999, max_treedepth = 50),
                     save_pars = save_pars(all = TRUE),
                     cores = 4)

helping_glmm2 <- brm(help_n_y ~ condition_string + age_group_string + location_string + first + gender_string +
                       (1|location_string),
                     data = helping,
                     family = bernoulli,
                     prior = prior_tendencies,
                     iter = 10000,
                     warmup = 2000,
                     control = list(adapt_delta = 0.9999, max_treedepth = 50),
                     save_pars = save_pars(all = TRUE),
                     cores = 4)

helping_glmm3 <- brm(help_n_y ~ condition_string * age_group_string + location_string + first + gender_string +
                       (1|location_string),
                     data = helping,
                     family = bernoulli,
                     prior = prior_tendencies,
                     iter = 10000,
                     warmup = 2000,
                     control = list(adapt_delta = 0.9999, max_treedepth = 50),
                     save_pars = save_pars(all = TRUE),
                     cores = 4)

plot(helping_glmm1)
plot(helping_glmm2)
plot(helping_glmm3)

pp_check(helping_glmm1, nsamples = 100)

post_prob_help123 <- post_prob(helping_glmm1, helping_glmm2, helping_glmm3)

BF_help_tendency <- 0.4159858/0.3345935
  
#Looking at age within we-framing
helping_glmm4 <- brm(help_n_y ~ age_group_string +
                       (1|location_string),
                     data = helping_tendency_we,
                     family = bernoulli,
                     prior = prior_tendencies,
                     iter = 10000,
                     warmup = 2000,
                     control = list(adapt_delta = 0.9999, max_treedepth = 50),
                     save_pars = save_pars(all = TRUE),
                     cores = 4)

plot(helping_glmm4)
pp_check(helping_glmm4, nsamples = 100)

#Looking at age within you-framing
helping_glmm5 <- brm(help_n_y ~ age_group_string +
                       (1|location_string),
                     data = helping_tendency_you,
                     family = bernoulli,
                     prior = prior_tendencies,
                     iter = 10000,
                     warmup = 2000,
                     control = list(adapt_delta = 0.9999, max_treedepth = 50),
                     save_pars = save_pars(all = TRUE),
                     cores = 4)

plot(helping_glmm5)
pp_check(helping_glmm5, nsamples = 100)

#Looking at condition within age [2.5,3.5]
helping_glmm6 <- brm(help_n_y ~ condition_string +
                       (1|location_string),
                     data = helping_tendency_3yo,
                     family = bernoulli,
                     prior = prior_tendencies,
                     iter = 10000,
                     warmup = 2000,
                     control = list(adapt_delta = 0.9999, max_treedepth = 50),
                     save_pars = save_pars(all = TRUE),
                     cores = 4)

plot(helping_glmm6)
pp_check(helping_glmm6, nsamples = 100)

#Looking at condition within age [3.5,4.5]
helping_glmm7 <- brm(help_n_y ~ condition_string +
                       (1|location_string),
                     data = helping_tendency_4yo,
                     family = bernoulli,
                     prior = prior_tendencies,
                     iter = 10000,
                     warmup = 2000,
                     control = list(adapt_delta = 0.9999, max_treedepth = 50),
                     save_pars = save_pars(all = TRUE),
                     cores = 4)

plot(helping_glmm7)
pp_check(helping_glmm7, nsamples = 100)

#Latency to help
helping_glmm8 <- brm(helping_latency ~ age_group_string + location_string + first + gender_string +
                       (1|location_string),
                     data = help_latency,
                     family = lognormal,
                     prior = prior_latencies,
                     iter = 10000,
                     warmup = 2000,
                     control = list(adapt_delta = 0.99999, max_treedepth = 60),
                     save_pars = save_pars(all = TRUE),
                     cores = 4)

helping_glmm9 <- brm(helping_latency ~ age_group_string + condition_string + location_string + first + gender_string +
                       (1|location_string),
                     data = help_latency,
                     family = lognormal,
                     prior = prior_latencies,
                     iter = 10000,
                     warmup = 2000,
                     control = list(adapt_delta = 0.99999, max_treedepth = 60),
                     save_pars = save_pars(all = TRUE),
                     cores = 4)

helping_glmm10 <- brm(helping_latency ~ age_group_string * condition_string + location_string + first + gender_string +
                        (1|location_string),
                      data = help_latency,
                      family = lognormal,
                      prior = prior_latencies,
                      iter = 10000,
                      warmup = 2000,
                      control = list(adapt_delta = 0.99999, max_treedepth = 60),
                      save_pars = save_pars(all = TRUE),
                      cores = 4)

plot(helping_glmm8)
plot(helping_glmm9)
plot(helping_glmm10)

pp_check(helping_glmm10, nsamples = 100)

post_prob_help8910 <- post_prob(helping_glmm8, helping_glmm9, helping_glmm10)

BF_help_latency <- 0.6692714/0.1179463

helping_glmm8_surv <- brm(helping_latency ~ age_group_string + location_string + first + gender_string,
                          data = help_latency_surv,
                          family = cox,
                          prior = prior_survival,
                          iter = 10000,
                          warmup = 2000,
                          control = list(adapt_delta = 0.9999, max_treedepth = 60),
                          save_pars = save_pars(all = TRUE),
                          cores = 4)

helping_glmm9_surv <- brm(helping_latency ~ age_group_string + condition_string + location_string + first + gender_string,   
                          data = help_latency_surv,
                          family = cox,
                          prior = prior_survival,
                          iter = 10000,
                          warmup = 2000,
                          control = list(adapt_delta = 0.9999, max_treedepth = 60),
                          save_pars = save_pars(all = TRUE),
                          cores = 4)

helping_glmm10_surv <- brm(helping_latency ~ age_group_string * condition_string + location_string + first + gender_string,
                           data = help_latency_surv,
                           family = cox,
                           prior = prior_survival,
                           iter = 10000,
                           warmup = 2000,
                           control = list(adapt_delta = 0.9999, max_treedepth = 60),
                           save_pars = save_pars(all = TRUE),
                           cores = 4)

plot(helping_glmm8_surv)
plot(helping_glmm9_surv)
plot(helping_glmm10_surv)

post_prob_help8910_surv <- post_prob(helping_glmm8_surv, helping_glmm9_surv, helping_glmm10_surv)

BF_help_surv <- 0.6008700/0.2573649

##Table 4
table_commit_glmm1 <- fixed(commit_glmm1)
table_commit_glmm8 <- fixef(commit_glmm8)
table_commit_glmm14 <- fixef(commit_glmm14)
table_sharing_glmm2 <- fixef(sharing_glmm2)
table_help_glmm1 <- fixef(helping_glmm1)
table_help_glmm8 <- fixef(helping_glmm8)

table_fixef <- rbind(table_commit_glmm1,
                     table_commit_glmm8,
                     table_commit_glmm14,
                     table_sharing_glmm2,
                     table_help_glmm1,
                     table_help_glmm8)

write.csv(table_fixef, "table_fixef.csv")

#Figure 2
ggplot(commit_tendency, aes(x = age_group_string, fill = factor(approach_n_y, levels = c(1,0)))) +
  geom_bar(position = "fill") +
  facet_grid(cols = vars(condition_string)) +
  labs(x = "Age group", y = "Proportion of participants", fill = "Participant decision") + 
  scale_fill_manual(values = c("lightblue", "yellowgreen"), labels = c("Abandon E1","Remain with E1"), limits = c("1", "0")) +
  geom_text(stat = "count", aes(y = ..count.., label = ..count..), position = "fill", vjust = 2, size = 6) +
  theme(text = element_text(size = 24),
        axis.title.y = element_text(vjust = 2))
ggsave("Figure 2.png", height = 10, width = 10)

#Figure 3
ggplot(leave_taking, aes(x = condition_string, fill = factor(leave_taking_type_recoded, levels = c(0,1)))) +
  geom_bar(position = "fill") +
  facet_grid(cols = vars(age_group_string)) +
  labs(x = "Age group", y = "Proportion of participants", fill = "Participant decision") + 
  scale_fill_manual(values = c("navajowhite2","aquamarine3"), labels = c("Spontaneously leave", "Take leave"), limits = c("0", "1")) +  geom_text(stat = "count", aes(label = ..count..), position = "fill", vjust = 2,  size = 6) +
  theme(text = element_text(size = 24),
        axis.title.y = element_text(vjust = 2))
ggsave("Figure 3.png", height = 10, width = 10)

#Supplementary Table 1
commit_glmm8_surv_df <- fixef(commit_glmm8_surv) %>% 
  as.data.frame %>% 
  mutate("exp(Estimate)" = exp(Estimate),
         "exp(Error)" = exp(Est.Error),
         "exp(Q2.5)" = exp(Q2.5),
         "exp(Q97.5)" = exp(Q97.5))

helping_glmm8_surv_df <- fixef(helping_glmm8_surv) %>% 
  as.data.frame %>% 
  mutate("exp(Estimate)" = exp(Estimate),
         "exp(Error)" = exp(Est.Error),
         "exp(Q2.5)" = exp(Q2.5),
         "exp(Q97.5)" = exp(Q97.5))
