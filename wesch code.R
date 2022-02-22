##Packages and seed
library(readxl)
library(tidyverse)
library(janitor)
library(brms)

set.seed(31)

##Data
wesch_data <- read_csv("wesch data.csv") %>% 
  clean_names() %>%
  mutate(condition_string = as.factor(condition_string),
         age_group_string = as.factor(age_group_string),
         location_string = as.factor(location_string),
         gender_string = as.factor(gender_string),
         first = as.factor(first)) %>% 
  rename(order = first) %>% 
  select(-c(x23, "x25":"x28")) %>% 
  na.omit()

wesch_data$age_group_string <- relevel(wesch_data$age_group_string, ref = "Younger")
wesch_data$condition_string <- relevel(wesch_data$condition_string, ref = "You")

##Priors
prior_tendencies <- c(
  set_prior("normal(0,1.25)", class = "Intercept"), #median = exp(0) = 1 baseline odds; -2sd = exp(-2.5) = 0.08 baseline odds; +2sd = 12.2 baseline odds
  set_prior("normal(0,1)", class = "b")) #median = exp(0) = 1 odds ratio; -2sd = exp(-2) = 0.13 odds ratio; +2sd = 7.39 odds ratio

qqplot_commitment <- wesch_data %>% 
  filter(commitment_latency < 61)

qqplot_helping <- wesch_data %>% 
  filter(helping_latency < 98)

ggplot(qqplot_commitment, aes(sample = commitment_latency)) +
  stat_qq(distribution = qnorm) +
  stat_qq_line() +
  ggtitle("Commitment latency, normal quantiles")

ggplot(qqplot_helping, aes(sample = helping_latency)) +
  stat_qq(distribution = qlnorm) +
  stat_qq_line() +
  ggtitle("Helping latency, normal quantiles")

prior_latencies <- c(
  set_prior("normal(3.25,0.25)", class = "Intercept"), #median = exp(3.25) = 25.8s; -2sd = exp(2.75) = 15.6s; +2sd = exp(3.75) = 42.5s
  set_prior("normal(0,0.4)", class = "b")) #median = exp(3.25)*exp(0) = 25.8s; -2sd = exp(3.25)*exp(-0.8) = 11.6s; +2sd = exp(3.25)*exp(-0.8) = 57.4s

prior_survival <- c(
  set_prior("normal(0,1)", class = "b")) #exp(0) = median hazard ratio of non-intercept predictors; hazard ratio of 1 = no change between levels of predictor;
                                         #exp(1) # = 1SD of hazard ratio of non-intercept predictors; coefficient estimates must be exponentiated to get the hazard ratio

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
  filter(approach_n_y != 0)

leave_taking <- wesch_data %>%
  filter(leave_taking_type < 4) %>% 
  mutate(leave_taking_verbal_and_nonverbal = ifelse(leave_taking_type %in% (1:2), 1, 0))

commit_leavetake_we <- leave_taking %>% 
  filter(condition_string != "you")

commit_leavetake_you <- leave_taking %>% 
  filter(condition_string != "we")

commit_leavetake_3yo <- leave_taking %>% 
  filter(age_years < 3.5)

commit_leavetake_4yo <- leave_taking %>% 
  filter(age_years > 3.5)

#Tendency to abandon
commit_glmm1 <- brm(approach_n_y ~ age_group_string + gender_string + (1|location_string) + (1|order),
                    data = commit_tendency,
                    family = bernoulli,
                    prior = prior_tendencies,
                    iter = 10000,
                    warmup = 2000,
                    control = list(adapt_delta = 0.9999, max_treedepth = 50),
                    save_pars = save_pars(all = TRUE),
                    cores = 4)

commit_glmm2 <- brm(approach_n_y ~ condition_string + age_group_string + gender_string + (1|location_string) + (1|order),
                    data = commit_tendency,
                    family = bernoulli,
                    prior = prior_tendencies,
                    iter = 10000,
                    warmup = 2000,
                    control = list(adapt_delta = 0.9999, max_treedepth = 50),
                    save_pars = save_pars(all = TRUE),
                    cores = 4)

commit_glmm3 <- brm(approach_n_y ~ condition_string * age_group_string + gender_string + (1|location_string) + (1|order),
                    data = commit_tendency,
                    family = bernoulli,
                    prior = prior_tendencies,
                    iter = 10000,
                    warmup = 2000,
                    control = list(adapt_delta = 0.9999, max_treedepth = 50),
                    save_pars = save_pars(all = TRUE),
                    cores = 4)

conditional_effects(commit_glmm1)
conditional_effects(commit_glmm2)
conditional_effects(commit_glmm3)

plot(commit_glmm1)
plot(commit_glmm2)
plot(commit_glmm3)

pp_check(commit_glmm1, ndraws = 100)
pp_check(commit_glmm2, ndraws = 100)
pp_check(commit_glmm3, ndraws = 100)

postprob_commit_tendency <- tibble("Posterior probability" = post_prob(commit_glmm1, commit_glmm2, commit_glmm3)) %>% 
  mutate("Model" = c("Null", "Reduced", "Full")) %>% 
  relocate("Model",.before = "Posterior probability")

BF_commit_tendency <- tibble("Bayes factor" = sort(postprob_commit_tendency$"Posterior probability",T)[1]/sort(postprob_commit_tendency$"Posterior probability",T)[2])

greaterthan0_committendency <- hypothesis(commit_glmm1,
                                          c("Intercept > 0",
                                            "age_group_stringOlder > 0",
                                            "gender_stringMale > 0"))

#Robustness analysis
postprob_commit_tendency0 <- tibble("Posterior probability" = post_prob(commit_glmm1, commit_glmm2, commit_glmm3, prior_prob = c(0.5, 0.25, 0.25))) %>% 
  mutate("Model" = c("Null", "Reduced", "Full")) %>% 
  relocate("Model",.before = "Posterior probability")

postprob_commit_tendency1 <- tibble("Posterior probability" = post_prob(commit_glmm1, commit_glmm2, commit_glmm3, prior_prob = c(0.25, 0.5, 0.25))) %>% 
  mutate("Model" = c("Null", "Reduced", "Full")) %>% 
  relocate("Model",.before = "Posterior probability")

postprob_commit_tendency2 <- tibble("Posterior probability" = post_prob(commit_glmm1, commit_glmm2, commit_glmm3, prior_prob = c(0.25, 0.25, 0.5))) %>% 
  mutate("Model" = c("Null", "Reduced", "Full")) %>% 
  relocate("Model",.before = "Posterior probability")

#Abandoning tendency: Looking at age within we-framing
commit_glmm4 <- brm(approach_n_y ~ age_group_string + (1|location_string) + (1|order),
                    data = commit_tendency_we,
                    family = bernoulli,
                    prior = prior_tendencies,
                    iter = 10000,
                    warmup = 2000,
                    control = list(adapt_delta = 0.9999, max_treedepth = 50),
                    save_pars = save_pars(all = TRUE),
                    cores = 4)

plot(commit_glmm4)
pp_check(commit_glmm4, ndraws = 100)

hypothesis(commit_glmm4, 
           "age_group_stringOlder < 0")

#Abandoning tendency: Looking at age within you-framing
commit_glmm5 <- brm(approach_n_y ~ age_group_string + (1|location_string) + (1|order),
                    data = commit_tendency_you,
                    family = bernoulli,
                    prior = prior_tendencies,
                    iter = 10000,
                    warmup = 2000,
                    control = list(adapt_delta = 0.9999, max_treedepth = 50),
                    save_pars = save_pars(all = TRUE),
                    cores = 4)

plot(commit_glmm5)
pp_check(commit_glmm5, ndraws = 100)

#Abandoning tendency: Looking at condition within age [2.5,3.5]
commit_glmm6 <- brm(approach_n_y ~ condition_string + (1|location_string) + (1|order),
                    data = commit_tendency_3yo,
                    family = bernoulli,
                    prior = prior_tendencies,
                    iter = 10000,
                    warmup = 2000,
                    control = list(adapt_delta = 0.9999, max_treedepth = 50),
                    save_pars = save_pars(all = TRUE),
                    cores = 4)

plot(commit_glmm6)
pp_check(commit_glmm6, ndraws = 100)

#Abandoning tendency: Looking at condition within age [3.5,4.5]
commit_glmm7 <- brm(approach_n_y ~ condition_string + (1|location_string) + (1|order),
                    data = commit_tendency_4yo,
                    family = bernoulli,
                    prior = prior_tendencies,
                    iter = 10000,
                    warmup = 2000,
                    control = list(adapt_delta = 0.9999, max_treedepth = 50),
                    save_pars = save_pars(all = TRUE),
                    cores = 4)

plot(commit_glmm7)
pp_check(commit_glmm7, ndraws = 100)

#Latency to abandon (Removed location and order random intercepts, produced divergences; included order fixed effect)
commit_glmm8 <- brm(commitment_latency|trunc(ub = 60) ~ age_group_string + gender_string + location_string,
                    data = commit_latency,
                    family = lognormal,
                    prior = prior_latencies,
                    iter = 10000,
                    warmup = 2000,
                    control = list(adapt_delta = 0.9999, max_treedepth = 50),
                    save_pars = save_pars(all = TRUE),
                    cores = 4)

commit_glmm9 <- brm(commitment_latency|trunc(ub = 60) ~ condition_string + age_group_string + gender_string + location_string,
                    data = commit_latency,
                    family = lognormal,
                    prior = prior_latencies,
                    iter = 10000,
                    warmup = 2000,
                    control = list(adapt_delta = 0.9999, max_treedepth = 50),
                    save_pars = save_pars(all = TRUE),
                    cores = 4)

commit_glmm10 <- brm(commitment_latency|trunc(ub = 60) ~ condition_string * age_group_string + gender_string + location_string,
                     data = commit_latency,
                     family = lognormal,
                     prior = prior_latencies,
                     iter = 10000,
                     warmup = 2000,
                     control = list(adapt_delta = 0.9999, max_treedepth = 50),
                     save_pars = save_pars(all = TRUE),
                     cores = 4)

conditional_effects(commit_glmm8)
conditional_effects(commit_glmm9)
conditional_effects(commit_glmm10)

plot(commit_glmm8)
plot(commit_glmm9)
plot(commit_glmm10)

pp_check(commit_glmm8, ndraws = 100)
pp_check(commit_glmm9, ndraws = 100)
pp_check(commit_glmm10, ndraws = 100)

postprob_commit_latency <- tibble("Posterior probability" = post_prob(commit_glmm8, commit_glmm9, commit_glmm10)) %>% 
  mutate("Model" = c("Null", "Reduced", "Full")) %>% 
  relocate("Model",.before = "Posterior probability")
BF_commit_latency <- tibble("Bayes factor" = sort(postprob_commit_latency$"Posterior probability",T)[1]/sort(postprob_commit_latency$"Posterior probability",T)[2])

greaterthan0_commitlatency <- hypothesis(commit_glmm8,
                                         c("Intercept > 0",
                                           "age_group_stringOlder > 0",
                                           "gender_stringMale > 0",
                                           "location_stringMarbles > 0",
                                           "location_stringMOLS > 0"))

#Robustness analysis
postprob_commit_latency0 <- tibble("Posterior probability" = post_prob(commit_glmm8, commit_glmm9, commit_glmm10, prior_prob = c(0.5, 0.25, 0.25))) %>% 
  mutate("Model" = c("Null", "Reduced", "Full")) %>% 
  relocate("Model",.before = "Posterior probability")

postprob_commit_latency1 <- tibble("Posterior probability" = post_prob(commit_glmm8, commit_glmm9, commit_glmm10, prior_prob = c(0.25, 0.5, 0.25))) %>% 
  mutate("Model" = c("Null", "Reduced", "Full")) %>% 
  relocate("Model",.before = "Posterior probability")

postprob_commit_latency2 <- tibble("Posterior probability" = post_prob(commit_glmm8, commit_glmm9, commit_glmm10, prior_prob = c(0.25, 0.25, 0.5))) %>% 
  mutate("Model" = c("Null", "Reduced", "Full")) %>% 
  relocate("Model",.before = "Posterior probability")

#Leave-taking
commit_glmm11 <- brm(leave_taking_verbal_and_nonverbal ~ age_group_string + gender_string + (1|location_string) + (1|order),
                     data = leave_taking,
                     family = bernoulli,
                     prior = prior_tendencies,
                     iter = 10000,
                     warmup = 2000,
                     control = list(adapt_delta = 0.9999, max_treedepth = 50),
                     save_pars = save_pars(all = TRUE),
                     cores = 4)

commit_glmm12 <- brm(leave_taking_verbal_and_nonverbal ~ condition_string + age_group_string + gender_string + (1|location_string) + (1|order),
                     data = leave_taking,
                     family = bernoulli,
                     prior = prior_tendencies,
                     iter = 10000,
                     warmup = 2000,
                     control = list(adapt_delta = 0.9999, max_treedepth = 50),
                     save_pars = save_pars(all = TRUE),
                     cores = 4)

commit_glmm13 <- brm(leave_taking_verbal_and_nonverbal ~ condition_string * age_group_string + gender_string + (1|location_string) + (1|order),
                     data = leave_taking,
                     family = bernoulli,
                     prior = prior_tendencies,
                     iter = 10000,
                     warmup = 2000,
                     control = list(adapt_delta = 0.9999, max_treedepth = 50),
                     save_pars = save_pars(all = TRUE),
                     cores = 4)

conditional_effects(commit_glmm11)
conditional_effects(commit_glmm12)
conditional_effects(commit_glmm13)

plot(commit_glmm11)
plot(commit_glmm12)
plot(commit_glmm13)

pp_check(commit_glmm11, ndraws = 100)
pp_check(commit_glmm12, ndraws = 100)
pp_check(commit_glmm13, ndraws = 100)

postprob_commit_takeleave <- tibble("Posterior probability" = post_prob(commit_glmm11, commit_glmm12, commit_glmm13)) %>% 
  mutate("Model" = c("Null", "Reduced", "Full")) %>% 
  relocate("Model", .before = "Posterior probability")
BF_commit_takeleave <- tibble("Bayes factor" = sort(postprob_commit_takeleave$"Posterior probability",T)[1]/sort(postprob_commit_takeleave$"Posterior probability",T)[2])

greaterthan0_takeleave <- hypothesis(commit_glmm13,
                                     c("Intercept > 0",
                                       "condition_stringWe > 0",
                                       "age_group_stringOlder > 0",
                                       "gender_stringMale > 0",
                                       "condition_stringWe:age_group_stringOlder > 0"))

#Robustness analysis
postprob_commit_takeleave0 <- tibble("Posterior probability" = post_prob(commit_glmm11, commit_glmm12, commit_glmm13, prior_prob = c(0.5, 0.25, 0.25))) %>% 
  mutate("Model" = c("Null", "Reduced", "Full")) %>% 
  relocate("Model",.before = "Posterior probability")

postprob_commit_takeleave1 <- tibble("Posterior probability" = post_prob(commit_glmm11, commit_glmm12, commit_glmm13, prior_prob = c(0.25, 0.5, 0.25))) %>% 
  mutate("Model" = c("Null", "Reduced", "Full")) %>% 
  relocate("Model",.before = "Posterior probability")

postprob_commit_takeleave2 <- tibble("Posterior probability" = post_prob(commit_glmm11, commit_glmm12, commit_glmm13, prior_prob = c(0.25, 0.25, 0.5))) %>% 
  mutate("Model" = c("Null", "Reduced", "Full")) %>% 
  relocate("Model",.before = "Posterior probability")

#Leave-taking tendency: Looking at condition within age [2.5,3.5]
commit_glmm14 <- brm(leave_taking_verbal_and_nonverbal ~ condition_string + (1|location_string) + (1|order),
                     data = commit_leavetake_3yo,
                     family = bernoulli,
                     prior = prior_tendencies,
                     iter = 10000,
                     warmup = 2000,
                     control = list(adapt_delta = 0.9999, max_treedepth = 50),
                     save_pars = save_pars(all = TRUE),
                     cores = 4)

plot(commit_glmm14)
pp_check(commit_glmm14, ndraws = 100)

hypothesis(commit_glmm14, "condition_stringWe > 0")

#Leave-taking tendency: Looking at condition within age [3.5,4.5]
commit_glmm15 <- brm(leave_taking_verbal_and_nonverbal ~ condition_string + (1|location_string) + (1|order),
                     data = commit_leavetake_4yo,
                     family = bernoulli,
                     prior = prior_tendencies,
                     iter = 10000,
                     warmup = 2000,
                     control = list(adapt_delta = 0.9999, max_treedepth = 50),
                     save_pars = save_pars(all = TRUE),
                     cores = 4)

plot(commit_glmm15)
pp_check(commit_glmm15, ndraws = 100)

hypothesis(commit_glmm15, "condition_stringWe > 0")

#Leave-taking tendency: Looking at age within we-framing
commit_glmm16 <- brm(leave_taking_verbal_and_nonverbal ~ age_group_string + (1|location_string) + (1|order),
                    data = commit_leavetake_we,
                    family = bernoulli,
                    prior = prior_tendencies,
                    iter = 10000,
                    warmup = 2000,
                    control = list(adapt_delta = 0.9999, max_treedepth = 50),
                    save_pars = save_pars(all = TRUE),
                    cores = 4)

plot(commit_glmm16)
pp_check(commit_glmm16, ndraws = 100)

hypothesis(commit_glmm16, "age_group_stringOlder > 0")

#Leave-taking tendency: Looking at age within you-framing
commit_glmm17 <- brm(leave_taking_verbal_and_nonverbal ~ age_group_string + (1|location_string) + (1|order),
                    data = commit_leavetake_you,
                    family = bernoulli,
                    prior = prior_tendencies,
                    iter = 10000,
                    warmup = 2000,
                    control = list(adapt_delta = 0.9999, max_treedepth = 50),
                    save_pars = save_pars(all = TRUE),
                    cores = 4)

plot(commit_glmm17)
pp_check(commit_glmm17, ndraws = 100)

hypothesis(commit_glmm17, "age_group_stringOlder > 0")

##Sharing
#Data
sharing <- wesch_data %>% 
  filter(number_of_erasers_shared < 7) %>% 
  mutate(equal_share = ifelse(number_of_erasers_shared == 3, 1, 0),
         generous_share = ifelse(number_of_erasers_shared > 3, 1, 0),
         equal_or_generous_share = ifelse(number_of_erasers_shared %in% (3:6), 1, 0))

share_equal_condition <- table(sharing$equal_share, sharing$condition_string)
share_generous_condition <- table(sharing$generous_share, sharing$condition_string)
share_equalorgenerous_condition <- table(sharing$equal_or_generous_share, sharing$condition_string)

#Resource distribution
sharing_glmm1 <- brm(equal_or_generous_share ~ age_group_string + gender_string + (1|location_string) + (1|order),
                     data = sharing,
                     family = bernoulli,
                     prior = prior_tendencies,
                     iter = 10000,
                     warmup = 2000,
                     control = list(adapt_delta = 0.9999, max_treedepth = 50),
                     save_pars = save_pars(all = TRUE),
                     cores = 4)

sharing_glmm2 <- brm(equal_or_generous_share ~ condition_string + age_group_string + gender_string + (1|location_string) + (1|order),
                     data = sharing,
                     family = bernoulli,
                     prior = prior_tendencies,
                     iter = 10000,
                     warmup = 2000,
                     control = list(adapt_delta = 0.9999, max_treedepth = 50),
                     save_pars = save_pars(all = TRUE),
                     cores = 4)

sharing_glmm3 <- brm(equal_or_generous_share ~ condition_string * age_group_string + gender_string + (1|location_string) + (1|order),
                     data = sharing,
                     family = bernoulli,
                     prior = prior_tendencies,
                     iter = 10000,
                     warmup = 2000,
                     control = list(adapt_delta = 0.9999, max_treedepth = 50),
                     save_pars = save_pars(all = TRUE),
                     cores = 4)

conditional_effects(sharing_glmm1)
conditional_effects(sharing_glmm2)
conditional_effects(sharing_glmm3)

plot(sharing_glmm1)
plot(sharing_glmm2)
plot(sharing_glmm3)

pp_check(sharing_glmm1, ndraws = 100)
pp_check(sharing_glmm2, ndraws = 100)
pp_check(sharing_glmm3, ndraws = 100)

postprob_share <- tibble("Posterior probability" = post_prob(sharing_glmm1, sharing_glmm2, sharing_glmm3)) %>% 
  mutate("Model" = c("Null", "Reduced", "Full")) %>% 
  relocate("Model",.before = "Posterior probability")
BF_share <- tibble("Bayes factor" = sort(postprob_share$"Posterior probability",T)[1]/sort(postprob_share$"Posterior probability",T)[2])

greaterthan0_sharing <- hypothesis(sharing_glmm1,
                                   c("Intercept > 0",
                                     "age_group_stringOlder > 0",
                                     "gender_stringMale > 0"))

#Robustness analysis
postprob_share0 <- tibble("Posterior probability" = post_prob(sharing_glmm1, sharing_glmm2, sharing_glmm3, prior_prob = c(0.5, 0.25, 0.25))) %>% 
  mutate("Model" = c("Null", "Reduced", "Full")) %>% 
  relocate("Model",.before = "Posterior probability")

postprob_share1 <- tibble("Posterior probability" = post_prob(sharing_glmm1, sharing_glmm2, sharing_glmm3, prior_prob = c(0.25, 0.5, 0.25))) %>% 
  mutate("Model" = c("Null", "Reduced", "Full")) %>% 
  relocate("Model",.before = "Posterior probability")

postprob_share2 <- tibble("Posterior probability" = post_prob(sharing_glmm1, sharing_glmm2, sharing_glmm3, prior_prob = c(0.25, 0.25, 0.5))) %>% 
  mutate("Model" = c("Null", "Reduced", "Full")) %>% 
  relocate("Model",.before = "Posterior probability")

#Sharing: Descriptive statistics
number_erasers_shared_you_3 <- sharing %>% 
  filter(condition_string != "We",
         age_group_string != "Older")
stats_number_erasers_shared_you_3 <- tibble(
   mean(number_erasers_shared_you_3$number_of_erasers_shared),
   sd(number_erasers_shared_you_3$number_of_erasers_shared),
   length(number_erasers_shared_you_3$number_of_erasers_shared)) %>% 
  rename("Mean" = "mean(number_erasers_shared_you_3$number_of_erasers_shared)",
         "SD" = "sd(number_erasers_shared_you_3$number_of_erasers_shared)",
         "N" = "length(number_erasers_shared_you_3$number_of_erasers_shared)")

number_erasers_shared_you_4 <- sharing %>% 
  filter(condition_string != "We",
         age_group_string != "Younger")
stats_number_erasers_shared_you_4 <- tibble(
  mean(number_erasers_shared_you_4$number_of_erasers_shared),
  sd(number_erasers_shared_you_4$number_of_erasers_shared),
  length(number_erasers_shared_you_4$number_of_erasers_shared)) %>% 
  rename("Mean" = "mean(number_erasers_shared_you_4$number_of_erasers_shared)",
         "SD" = "sd(number_erasers_shared_you_4$number_of_erasers_shared)",
         "N" = "length(number_erasers_shared_you_4$number_of_erasers_shared)")

number_erasers_shared_we_3 <- sharing %>% 
  filter(condition_string != "You",
         age_group_string != "Older")
stats_number_erasers_shared_we_3 <- tibble(
  mean(number_erasers_shared_we_3$number_of_erasers_shared),
  sd(number_erasers_shared_we_3$number_of_erasers_shared),
  length(number_erasers_shared_we_3$number_of_erasers_shared)) %>% 
  rename("Mean" = "mean(number_erasers_shared_we_3$number_of_erasers_shared)",
         "SD" = "sd(number_erasers_shared_we_3$number_of_erasers_shared)",
         "N" = "length(number_erasers_shared_we_3$number_of_erasers_shared)")

number_erasers_shared_we_4 <- sharing %>% 
  filter(condition_string != "You",
         age_group_string != "Younger")
stats_number_erasers_shared_we_4 <- tibble(
  mean(number_erasers_shared_we_4$number_of_erasers_shared),
  sd(number_erasers_shared_we_4$number_of_erasers_shared),
  length(number_erasers_shared_we_4$number_of_erasers_shared)) %>% 
  rename("Mean" = "mean(number_erasers_shared_we_4$number_of_erasers_shared)",
         "SD" = "sd(number_erasers_shared_we_4$number_of_erasers_shared)",
         "N" = "length(number_erasers_shared_we_4$number_of_erasers_shared)")

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
  select(censored, age_group_string, helping_latency, age_stndrd, condition_string, location_string, order, gender_string) %>% 
  relocate(censored, .after = helping_latency) %>% 
  filter(helping_latency != 60)

#Tendency to help
helping_glmm1 <- brm(help_n_y ~ age_group_string + gender_string + (1|location_string) + (1|order),
                     data = helping,
                     family = bernoulli,
                     prior = prior_tendencies,
                     iter = 10000,
                     warmup = 2000,
                     control = list(adapt_delta = 0.9999, max_treedepth = 50),
                     save_pars = save_pars(all = TRUE),
                     cores = 4)

helping_glmm2 <- brm(help_n_y ~ condition_string + age_group_string + gender_string + (1|location_string) + (1|order),
                     data = helping,
                     family = bernoulli,
                     prior = prior_tendencies,
                     iter = 10000,
                     warmup = 2000,
                     control = list(adapt_delta = 0.9999, max_treedepth = 50),
                     save_pars = save_pars(all = TRUE),
                     cores = 4)

helping_glmm3 <- brm(help_n_y ~ condition_string * age_group_string + gender_string + (1|location_string) + (1|order),
                     data = helping,
                     family = bernoulli,
                     prior = prior_tendencies,
                     iter = 10000,
                     warmup = 2000,
                     control = list(adapt_delta = 0.9999, max_treedepth = 50),
                     save_pars = save_pars(all = TRUE),
                     cores = 4)

conditional_effects(helping_glmm1)
conditional_effects(helping_glmm2)
conditional_effects(helping_glmm3)

plot(helping_glmm1)
plot(helping_glmm2)
plot(helping_glmm3)

pp_check(helping_glmm1, ndraws = 100)
pp_check(helping_glmm2, ndraws = 100)
pp_check(helping_glmm3, ndraws = 100)

postprob_help_tendency <- tibble("Posterior probability" = post_prob(helping_glmm1, helping_glmm2, helping_glmm3)) %>% 
  mutate("Model" = c("Null", "Reduced", "Full")) %>% 
  relocate("Model",.before = "Posterior probability")
BF_help_tendency <- tibble("Bayes factor" = sort(postprob_help_tendency$"Posterior probability",T)[1]/sort(postprob_help_tendency$"Posterior probability",T)[2])

greaterthan0_helptendency <- hypothesis(helping_glmm1,
                                        c("Intercept > 0",
                                          "age_group_stringOlder > 0",
                                          "gender_stringMale > 0"))

#Robustness analysis
postprob_help_tendency0 <- tibble("Posterior probability" = post_prob(helping_glmm1, helping_glmm2, helping_glmm3, prior_prob = c(0.5, 0.25, 0.25))) %>% 
  mutate("Model" = c("Null", "Reduced", "Full")) %>% 
  relocate("Model",.before = "Posterior probability")

postprob_help_tendency1 <- tibble("Posterior probability" = post_prob(helping_glmm1, helping_glmm2, helping_glmm3, prior_prob = c(0.25, 0.5, 0.25))) %>% 
  mutate("Model" = c("Null", "Reduced", "Full")) %>% 
  relocate("Model",.before = "Posterior probability")

postprob_help_tendency2 <- tibble("Posterior probability" = post_prob(helping_glmm1, helping_glmm2, helping_glmm3, prior_prob = c(0.25, 0.25, 0.5))) %>% 
  mutate("Model" = c("Null", "Reduced", "Full")) %>% 
  relocate("Model",.before = "Posterior probability")

help_tendency_age <- table(helping$age_group_string, helping$help_n_y)
help_tendency_condition <- table(helping$condition_string, helping$help_n_y)
help_tendency_gender <- table(helping$gender_string, helping$help_n_y)

#Latency to help (Removed location and order random intercepts, produced divergences; included order fixed effect)
helping_glmm4 <- brm(helping_latency|trunc(ub = 60) ~ age_group_string + gender_string + location_string,
                     data = help_latency,
                     family = lognormal,
                     prior = prior_latencies,
                     iter = 10000,
                     warmup = 2000,
                     control = list(adapt_delta = 0.9999, max_treedepth = 50),
                     save_pars = save_pars(all = TRUE),
                     cores = 4)

helping_glmm5 <- brm(helping_latency|trunc(ub = 60) ~ condition_string + age_group_string + gender_string + location_string,
                     data = help_latency,
                     family = lognormal,
                     prior = prior_latencies,
                     iter = 10000,
                     warmup = 2000,
                     control = list(adapt_delta = 0.9999, max_treedepth = 50),
                     save_pars = save_pars(all = TRUE),
                     cores = 4)

helping_glmm6 <- brm(helping_latency|trunc(ub = 60) ~ condition_string * age_group_string + gender_string + location_string,
                      data = help_latency,
                      family = lognormal,
                      prior = prior_latencies,
                      iter = 10000,
                      warmup = 2000,
                      control = list(adapt_delta = 0.9999, max_treedepth = 50),
                      save_pars = save_pars(all = TRUE),
                      cores = 4)

conditional_effects(helping_glmm4)
conditional_effects(helping_glmm5)
conditional_effects(helping_glmm6)

plot(helping_glmm4)
plot(helping_glmm5)
plot(helping_glmm6)

pp_check(helping_glmm4, ndraws = 100)
pp_check(helping_glmm5, ndraws = 100)
pp_check(helping_glmm6, ndraws = 100)

postprob_help_latency <- tibble("Posterior probability" = post_prob(helping_glmm4, helping_glmm5, helping_glmm6)) %>% 
  mutate("Model" = c("Null", "Reduced", "Full")) %>% 
  relocate("Model",.before = "Posterior probability")
BF_help_latency <- tibble("Bayes factor" = sort(postprob_help_latency$"Posterior probability",T)[1]/sort(postprob_help_latency$"Posterior probability",T)[2])

greaterthan0_helplatency <- hypothesis(helping_glmm4,
                                       c("Intercept > 0",
                                         "age_group_stringOlder > 0",
                                         "gender_stringMale > 0",
                                         "location_stringMarbles > 0",
                                         "location_stringMOLS > 0"))

#Robustness analysis
postprob_help_latency0 <- tibble("Posterior probability" = post_prob(helping_glmm4, helping_glmm5, helping_glmm6, prior_prob = c(0.5, 0.25, 0.25))) %>% 
  mutate("Model" = c("Null", "Reduced", "Full")) %>% 
  relocate("Model",.before = "Posterior probability")

postprob_help_latency1 <- tibble("Posterior probability" = post_prob(helping_glmm4, helping_glmm5, helping_glmm6, prior_prob = c(0.25, 0.5, 0.25))) %>% 
  mutate("Model" = c("Null", "Reduced", "Full")) %>% 
  relocate("Model",.before = "Posterior probability")

postprob_help_latency2 <- tibble("Posterior probability" = post_prob(helping_glmm4, helping_glmm5, helping_glmm6, prior_prob = c(0.25, 0.25, 0.5))) %>% 
  mutate("Model" = c("Null", "Reduced", "Full")) %>% 
  relocate("Model",.before = "Posterior probability")

##Main Text: Table 4
table_fixef <- rbind(fixef(commit_glmm1),
                     fixef(commit_glmm8),
                     fixef(commit_glmm13),
                     fixef(sharing_glmm1),
                     fixef(helping_glmm1),
                     fixef(helping_glmm4)) %>% 
  round(digits = 2)

write.csv(table_fixef, "Table 4.csv")

##Main Text: Figure 2
fig2_commit_tendency <- commit_tendency %>% 
  mutate("asterisk" = ifelse(condition_string == "We", "*", NA),
         condition_string = fct_recode(condition_string,
                                       "We-framing" = "We",
                                       "You-framing" = "You"))

ggplot(fig2_commit_tendency, aes(x = age_group_string, fill = factor(approach_n_y, levels = c(1,0)))) +
  geom_bar(position = "fill") +
  facet_grid(cols = vars(condition_string)) +
  labs(x = "Age group", y = "Proportion of participants", fill = "Participant decision") + 
  scale_fill_manual(values = c("lightblue", "limegreen"), labels = c("Abandon E1","Remain with E1"), limits = c("1", "0")) +
  geom_text(stat = "count", aes(y = ..count.., label = ..count..), position = "fill", vjust = 1.5, size = 8) +
  geom_text(aes(x = 1.5,  y = 1.03, label = asterisk), size = 8) + 
  theme(text = element_text(size = 30),
        axis.title.y = element_text(vjust = 2),
        axis.title.x = element_text(vjust = 0.25))

ggsave("figures/figure_2.png", height = 14, width = 14)

##Main Text: Figure 3
fig3_leave_taking <- leave_taking %>% 
  mutate("dagger" = ifelse(age_group_string == "Younger", "*", NA),
         condition_string = fct_recode(condition_string,
                                       "We-framing" = "We",
                                       "You-framing" = "You"))

ggplot(fig3_leave_taking, aes(x = condition_string, fill = factor(leave_taking_verbal_and_nonverbal, levels = c(0,1)))) +
  geom_bar(position = "fill") +
  facet_grid(cols = vars(age_group_string)) +
  labs(x = "Condition", y = "Proportion of participants", fill = "Participant decision") + 
  scale_fill_manual(values = c("navajowhite2","palegreen4"), labels = c("Spontaneously leave", "Take leave"), limits = c("0", "1")) +
  geom_text(stat = "count", aes(label = ..count..), position = "fill", vjust = 1.5,  size = 8) +
  geom_text(aes(x = 1.5,  y = 1.03, label = dagger), size = 8) + 
  theme(text = element_text(size = 30),
        axis.title.y = element_text(vjust = 2),
        axis.title.x = element_text(vjust = 0.25))

ggsave("figures/figure_3.png", height = 14, width = 14)

 ##Supplementary Material
#Data
commit_latency_surv <- wesch_data %>% 
  mutate(censored = ifelse(commitment_latency > 60, "right", "none"),
         commitment_latency = ifelse(censored == "right", 60, commitment_latency)) %>% 
  select(censored, commitment_latency, age_stndrd, age_group_string, condition_string, location_string, order, gender_string) %>% 
  relocate(censored, .after = commitment_latency)

help_latency_surv <- helping %>% 
  mutate(censored = ifelse(helping_latency > 60, "right", "none"),
         helping_latency = ifelse(censored == "right", 60, helping_latency)) %>% 
  select(censored, age_group_string, helping_latency, age_stndrd, condition_string, location_string, order, gender_string) %>% 
  relocate(censored, .after = helping_latency)

#Commitment: Cox PH models (Removed order and location random intercepts; included location fixed effect but not order fixed effect because too many levels of order)
commit_glmm1_surv <- brm(commitment_latency | cens(censored) ~ age_group_string + gender_string + location_string,
                         data = commit_latency_surv,
                         family = cox,
                         prior = prior_survival,
                         iter = 10000,
                         warmup = 2000,
                         control = list(adapt_delta = 0.9999, max_treedepth = 50),
                         save_pars = save_pars(all = TRUE),
                         cores = 4)

commit_glmm2_surv <- brm(commitment_latency | cens(censored) ~ condition_string + age_group_string + gender_string + location_string,
                         data = commit_latency_surv,
                         family = cox,
                         prior = prior_survival,
                         iter = 10000,
                         warmup = 2000,
                         control = list(adapt_delta = 0.9999, max_treedepth = 50),
                         save_pars = save_pars(all = TRUE),
                         cores = 4)

commit_glmm3_surv <- brm(commitment_latency | cens(censored) ~ condition_string * age_group_string + gender_string + location_string,
                          data = commit_latency_surv,
                          family = cox,
                          prior = prior_survival,
                          iter = 10000,
                          warmup = 2000,
                          control = list(adapt_delta = 0.9999, max_treedepth = 50),
                          save_pars = save_pars(all = TRUE),
                          cores = 4)

plot(commit_glmm1_surv)
plot(commit_glmm2_surv)
plot(commit_glmm3_surv)

postprob_commit_surv <- tibble("Posterior probability" = post_prob(commit_glmm1_surv, commit_glmm2_surv, commit_glmm3_surv)) %>% 
  mutate("Model" = c("Null", "Reduced", "Full")) %>% 
  relocate("Model",.before = "Posterior probability")

BF_commit_surv <- tibble("Bayes factor" = sort(postprob_commit_surv$"Posterior probability",T)[1]/sort(postprob_commit_surv$"Posterior probability",T)[2])

greaterthan0_commitsurv <- hypothesis(commit_glmm1_surv,
                                      c("Intercept > 0",
                                        "age_group_stringOlder > 0",
                                        "gender_stringMale > 0",
                                        "location_stringMarbles > 0",
                                        "location_stringMOLS > 0"))

#Helping: Cox PH models (Removed order and location random intercepts; included location fixed effect but not order fixed effect because too many levels of order)
helping_glmm1_surv <- brm(helping_latency | cens(censored) ~ age_group_string + gender_string + location_string,
                          data = help_latency_surv,
                          family = cox,
                          prior = prior_survival,
                          iter = 10000,
                          warmup = 2000,
                          control = list(adapt_delta = 0.9999, max_treedepth = 50),
                          save_pars = save_pars(all = TRUE),
                          cores = 4)

helping_glmm2_surv <- brm(helping_latency | cens(censored) ~ condition_string + age_group_string + gender_string + location_string,
                          data = help_latency_surv,
                          family = cox,
                          prior = prior_survival,
                          iter = 10000,
                          warmup = 2000,
                          control = list(adapt_delta = 0.9999, max_treedepth = 50),
                          save_pars = save_pars(all = TRUE),
                          cores = 4)

helping_glmm3_surv <- brm(helping_latency | cens(censored) ~ condition_string * age_group_string + gender_string + location_string,
                           data = help_latency_surv,
                           family = cox,
                           prior = prior_survival,
                           iter = 10000,
                           warmup = 2000,
                           control = list(adapt_delta = 0.9999, max_treedepth = 50),
                           save_pars = save_pars(all = TRUE),
                           cores = 4)

plot(helping_glmm1_surv)
plot(helping_glmm2_surv)
plot(helping_glmm3_surv)

postprob_help_surv <- tibble("Posterior probability" = post_prob(helping_glmm1_surv, helping_glmm2_surv, helping_glmm3_surv)) %>% 
  mutate("Model" = c("Null", "Reduced", "Full")) %>% 
  relocate("Model",.before = "Posterior probability")
BF_help_surv <- tibble("Bayes factor" = sort(postprob_help_surv$"Posterior probability",T)[1]/sort(postprob_help_surv$"Posterior probability",T)[2])

greaterthan0_helpsurv <- hypothesis(helping_glmm1_surv,
                                    c("Intercept > 0",
                                      "age_group_stringOlder > 0",
                                      "gender_stringMale > 0",
                                      "location_stringMarbles > 0",
                                      "location_stringMOLS > 0"))

#Supplementary Material Table 1
commit_glmm1_surv_df <- fixef(commit_glmm1_surv) %>% 
  as.data.frame %>% 
  mutate("exp(Estimate)" = exp(Estimate),
         "exp(Error)" = exp(Est.Error),
         "exp(Q2.5)" = exp(Q2.5),
         "exp(Q97.5)" = exp(Q97.5)) %>% 
  select("exp(Estimate)":"exp(Q97.5)") %>% 
  round(digits = 2)

  as.data.frame %>% 
  mutate("exp(Estimate)" = exp(Estimate),
         "exp(Error)" = exp(Est.Error),
         "exp(Q2.5)" = exp(Q2.5),
         "exp(Q97.5)" = exp(Q97.5)) %>% 
  select("exp(Estimate)":"exp(Q97.5)") %>% 
  round(digits = 2)