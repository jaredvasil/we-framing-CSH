##Packages and seed
library(readxl)
library(tidyverse)
library(janitor)
library(brms)

set.seed(31)

##Data
wesch_data <- read_xlsx("wesch data.xlsx") %>% 
  clean_names() %>% 
  mutate(condition_string = as.factor(condition_string),
         age_group_string = as.factor(age_group_string),
         location_string = as.factor(location_string),
         gender_string = as.factor(gender_string))

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
  set_prior("normal(0,1)", class = "b")) # exp(0) = median hazard ratio of non-intercept predictors; hazard ratio of 1 = no change between levels of predictor; exp(1) # = 1SD of hazard ratio of non-intercept predictors; coefficient estimates must be exponentiated to get the hazard ratio

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

commit_latency_surv <- wesch_data %>% 
  mutate(censored = ifelse(commitment_latency > 60, "right", "none"),
         commitment_latency = ifelse(censored == "right", 60, commitment_latency)) %>% 
  select(censored, commitment_latency, age_stndrd, age_group_string, condition_string, location_string, first, gender_string) %>% 
  relocate(censored, .after = commitment_latency)

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

conditional_effects(commit_glmm1)
conditional_effects(commit_glmm2)
conditional_effects(commit_glmm3)

plot(commit_glmm1)
plot(commit_glmm2)
plot(commit_glmm3)

pp_check(commit_glmm1, nsamples = 100)
pp_check(commit_glmm2, nsamples = 100)
pp_check(commit_glmm3, nsamples = 100)

postprob_commit_tendency <- tibble("Posterior probability" = post_prob(commit_glmm1, commit_glmm2, commit_glmm3)) %>% 
  mutate("Model" = c("Null", "Reduced", "Full")) %>% 
  relocate("Model",.before = "Posterior probability")
BF_commit_tendency <- tibble("Bayes factor" = sort(postprob_commit_tendency$"Posterior probability",T)[1]/sort(postprob_commit_tendency$"Posterior probability",T)[2])

loo_commit1 <- loo(commit_glmm1)
loo_commit2 <- loo(commit_glmm2)
loo_commit3 <- loo(commit_glmm3)
loo_commit_tendency <- loo_compare(loo_commit1, loo_commit2, loo_commit3)

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
pp_check(commit_glmm6, nsamples = 100)

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
commit_glmm8 <- brm(commitment_latency ~ age_group_string + location_string + first + gender_string,
                    data = commit_latency,
                    family = lognormal,
                    prior = prior_latencies,
                    iter = 10000,
                    warmup = 2000,
                    control = list(adapt_delta = 0.9999, max_treedepth = 50),
                    save_pars = save_pars(all = TRUE),
                    cores = 4)

commit_glmm9 <- brm(commitment_latency ~ age_group_string + condition_string + location_string + first + gender_string,
                    data = commit_latency,
                    family = lognormal,
                    prior = prior_latencies,
                    iter = 10000,
                    warmup = 2000,
                    control = list(adapt_delta = 0.9999, max_treedepth = 50),
                    save_pars = save_pars(all = TRUE),
                    cores = 4)

commit_glmm10 <- brm(commitment_latency ~ age_group_string * condition_string + location_string + first + gender_string,
                     data = commit_latency,
                     family = lognormal,
                     prior = prior_latencies,
                     iter = 10000,
                     warmup = 2000,
                     control = list(adapt_delta = 0.9999, max_treedepth = 50),
                     save_pars = save_pars(all = TRUE),
                     cores = 4)

#using truncated distribution for prior predictive checks provides reasonable estimates.
#the bridge sampler (used to compute the posterior probability) tends to give an error
#when computing the posterior probability of truncated distributions, and the posterior
#probabilities seem implausible when the bridge sampler does work. i think the best
#method currently available is to use the truncated model (ub = 60) to generate simulated
#data used in prior predictive checks, and then fit the non-truncated model to the data.

conditional_effects(commit_glmm8)
conditional_effects(commit_glmm9)
conditional_effects(commit_glmm10)

plot(commit_glmm8)
plot(commit_glmm9)
plot(commit_glmm10)

pp_check(commit_glmm8, nsamples = 100)
pp_check(commit_glmm9, nsamples = 100)
pp_check(commit_glmm10, nsamples = 100)

postprob_commit_latency <- tibble("Posterior probability" = post_prob(commit_glmm8, commit_glmm9, commit_glmm10)) %>% 
  mutate("Model" = c("Null", "Reduced", "Full")) %>% 
  relocate("Model",.before = "Posterior probability")
BF_commit_latency <- tibble("Bayes factor" = sort(postprob_commit_latency$"Posterior probability",T)[1]/sort(postprob_commit_latency$"Posterior probability",T)[2])

loo_commit8 <- loo(commit_glmm8)
loo_commit9 <- loo(commit_glmm9)
loo_commit10 <- loo(commit_glmm10)
loo_commit_latency <- loo_compare(loo_commit8, loo_commit9, loo_commit10)

#Leave-taking
commit_glmm11 <- brm(leave_taking_verbal_and_nonverbal ~ age_group_string + location_string + first + gender_string +
                       (1|location_string),
                     data = leave_taking,
                     family = bernoulli,
                     prior = prior_tendencies,
                     iter = 10000,
                     warmup = 2000,
                     control = list(adapt_delta = 0.9999, max_treedepth = 50),
                     save_pars = save_pars(all = TRUE),
                     cores = 4)

commit_glmm12 <- brm(leave_taking_verbal_and_nonverbal ~ age_group_string + condition_string + location_string + first + gender_string +
                       (1|location_string),
                     data = leave_taking,
                     family = bernoulli,
                     prior = prior_tendencies,
                     iter = 10000,
                     warmup = 2000,
                     control = list(adapt_delta = 0.9999, max_treedepth = 50),
                     save_pars = save_pars(all = TRUE),
                     cores = 4)

commit_glmm13 <- brm(leave_taking_verbal_and_nonverbal ~ age_group_string * condition_string + location_string + first + gender_string +
                       (1|location_string),
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

pp_check(commit_glmm11, nsamples = 100)
pp_check(commit_glmm12, nsamples = 100)
pp_check(commit_glmm13, nsamples = 100)

postprob_commit_takeleave <- tibble("Posterior probability" = post_prob(commit_glmm11, commit_glmm12, commit_glmm13)) %>% 
  mutate("Model" = c("Null", "Reduced", "Full")) %>% 
  relocate("Model", .before = "Posterior probability")
BF_commit_takeleave <- tibble("Bayes factor" = sort(postprob_commit_takeleave$"Posterior probability",T)[1]/sort(postprob_commit_takeleave$"Posterior probability",T)[2])

loo_commit11 <- loo(commit_glmm11)
loo_commit12 <- loo(commit_glmm12)
loo_commit13 <- loo(commit_glmm13)
loo_commit_takeleave <- loo_compare(loo_commit11, loo_commit12, loo_commit13)

postprob_greaterthan0_takeleave <- hypothesis(commit_glmm13, "condition_stringWe > 0")

contingency_table_order_takeleave <- table(leave_taking$first, leave_taking$leave_taking_verbal_and_nonverbal)

leave_taking_noorder3 <- leave_taking %>% 
  filter(first != 5)

commit_glmm13_noorder3 <- brm(leave_taking_verbal_and_nonverbal ~ age_group_string * condition_string + location_string + first + gender_string +
                       (1|location_string),
                       data = leave_taking_noorder3,
                       family = bernoulli,
                       prior = prior_tendencies,
                       iter = 10000,
                       warmup = 2000,
                       control = list(adapt_delta = 0.9999, max_treedepth = 50),
                       save_pars = save_pars(all = TRUE),
                       cores = 4)

#Leave-taking tendency: Looking at condition within age [2.5,3.5]
commit_glmm14 <- brm(leave_taking_verbal_and_nonverbal ~ condition_string +
                       (1|location_string),
                     data = commit_leavetake_3yo,
                     family = bernoulli,
                     prior = prior_tendencies,
                     iter = 10000,
                     warmup = 2000,
                     control = list(adapt_delta = 0.9999, max_treedepth = 50),
                     save_pars = save_pars(all = TRUE),
                     cores = 4)

plot(commit_glmm14)
pp_check(commit_glmm14, nsamples = 100)

hypothesis(commit_glmm14, "condition_stringWe > 0")

#Leave-taking tendency: Looking at condition within age [3.5,4.5]
commit_glmm15 <- brm(leave_taking_verbal_and_nonverbal ~ condition_string +
                       (1|location_string),
                     data = commit_leavetake_4yo,
                     family = bernoulli,
                     prior = prior_tendencies,
                     iter = 10000,
                     warmup = 2000,
                     control = list(adapt_delta = 0.9999, max_treedepth = 50),
                     save_pars = save_pars(all = TRUE),
                     cores = 4)

plot(commit_glmm15)
pp_check(commit_glmm15, nsamples = 100)

#Leave-taking tendency: Looking at age within we-framing
commit_glmm16 <- brm(leave_taking_verbal_and_nonverbal ~ age_group_string +
                       (1|location_string),
                    data = commit_leavetake_we,
                    family = bernoulli,
                    prior = prior_tendencies,
                    iter = 10000,
                    warmup = 2000,
                    control = list(adapt_delta = 0.9999, max_treedepth = 50),
                    save_pars = save_pars(all = TRUE),
                    cores = 4)

plot(commit_glmm16)
pp_check(commit_glmm16, nsamples = 100)

#Leave-taking tendency: Looking at age within you-framing
commit_glmm17 <- brm(leave_taking_verbal_and_nonverbal ~ age_group_string +
                       (1|location_string),
                    data = commit_leavetake_you,
                    family = bernoulli,
                    prior = prior_tendencies,
                    iter = 10000,
                    warmup = 2000,
                    control = list(adapt_delta = 0.9999, max_treedepth = 50),
                    save_pars = save_pars(all = TRUE),
                    cores = 4)

plot(commit_glmm17)
pp_check(commit_glmm17, nsamples = 100)

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
sharing_glmm1 <- brm(equal_or_generous_share ~ age_group_string + location_string + first + gender_string +
                       (1|location_string),
                     data = sharing,
                     family = bernoulli,
                     prior = prior_tendencies,
                     iter = 10000,
                     warmup = 2000,
                     control = list(adapt_delta = 0.9999, max_treedepth = 50),
                     save_pars = save_pars(all = TRUE),
                     cores = 4)

sharing_glmm2 <- brm(equal_or_generous_share ~ age_group_string + condition_string + location_string + first + gender_string +
                       (1|location_string),
                     data = sharing,
                     family = bernoulli,
                     prior = prior_tendencies,
                     iter = 10000,
                     warmup = 2000,
                     control = list(adapt_delta = 0.9999, max_treedepth = 50),
                     save_pars = save_pars(all = TRUE),
                     cores = 4)

sharing_glmm3 <- brm(equal_or_generous_share ~ age_group_string * condition_string + location_string + first + gender_string +
                       (1|location_string),
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

pp_check(sharing_glmm1, nsamples = 100)
pp_check(sharing_glmm2, nsamples = 100)
pp_check(sharing_glmm3, nsamples = 100)

postprob_share <- tibble("Posterior probability" = post_prob(sharing_glmm1, sharing_glmm2, sharing_glmm3)) %>% 
  mutate("Model" = c("Null", "Reduced", "Full")) %>% 
  relocate("Model",.before = "Posterior probability")
BF_share <- tibble("Bayes factor" = sort(postprob_share$"Posterior probability",T)[1]/sort(postprob_share$"Posterior probability",T)[2])

loo_share1 <- loo(sharing_glmm1)
loo_share2 <- loo(sharing_glmm2)
loo_share3 <- loo(sharing_glmm3)
loo_share <- loo_compare(loo_share1, loo_share2, loo_share3)

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
  select(censored, age_group_string, helping_latency, age_stndrd, condition_string, location_string, first, gender_string) %>% 
  relocate(censored, .after = helping_latency) %>% 
  filter(helping_latency != 60)

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

conditional_effects(helping_glmm1)
conditional_effects(helping_glmm2)
conditional_effects(helping_glmm3)

plot(helping_glmm1)
plot(helping_glmm2)
plot(helping_glmm3)

pp_check(helping_glmm1, nsamples = 100)
pp_check(helping_glmm2, nsamples = 100)
pp_check(helping_glmm3, nsamples = 100)

postprob_help_tendency <- tibble("Posterior probability" = post_prob(helping_glmm1, helping_glmm2, helping_glmm3)) %>% 
  mutate("Model" = c("Null", "Reduced", "Full")) %>% 
  relocate("Model",.before = "Posterior probability")
BF_help_tendency <- tibble("Bayes factor" = sort(postprob_help_tendency$"Posterior probability",T)[1]/sort(postprob_help_tendency$"Posterior probability",T)[2])

loo_help1 <- loo(helping_glmm1)
loo_help2 <- loo(helping_glmm2)
loo_help3 <- loo(helping_glmm3)
loo_help_tendency <- loo_compare(loo_help1, loo_help2, loo_help3)

help_tendency_age <- table(helping$age_group_string, helping$help_n_y)
help_tendency_condition <- table(helping$condition_string, helping$help_n_y)
help_tendency_gender <- table(helping$gender_string, helping$help_n_y)

#Helping tendency: Looking at age within we-framing
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

#Helping tendency: Looking at age within you-framing
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

#Helping tendency: Looking at condition within age [2.5,3.5]
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

#Helping tendency: Looking at condition within age [3.5,4.5]
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
helping_glmm8 <- brm(helping_latency ~ age_group_string + location_string + first + gender_string,
                     data = help_latency,
                     family = lognormal,
                     prior = prior_latencies,
                     iter = 10000,
                     warmup = 2000,
                     control = list(adapt_delta = 0.9999, max_treedepth = 50),
                     save_pars = save_pars(all = TRUE),
                     cores = 4)

helping_glmm9 <- brm(helping_latency ~ age_group_string + condition_string + location_string + first + gender_string,
                     data = help_latency,
                     family = lognormal,
                     prior = prior_latencies,
                     iter = 10000,
                     warmup = 2000,
                     control = list(adapt_delta = 0.9999, max_treedepth = 50),
                     save_pars = save_pars(all = TRUE),
                     cores = 4)

helping_glmm10 <- brm(helping_latency ~ age_group_string * condition_string + location_string + first + gender_string,
                      data = help_latency,
                      family = lognormal,
                      prior = prior_latencies,
                      iter = 10000,
                      warmup = 2000,
                      control = list(adapt_delta = 0.9999, max_treedepth = 50),
                      save_pars = save_pars(all = TRUE),
                      cores = 4)

conditional_effects(helping_glmm8)
conditional_effects(helping_glmm9)
conditional_effects(helping_glmm10)

plot(helping_glmm8)
plot(helping_glmm9)
plot(helping_glmm10)

pp_check(helping_glmm8, nsamples = 100)
pp_check(helping_glmm9, nsamples = 100)
pp_check(helping_glmm10, nsamples = 100)

postprob_help_latency <- tibble("Posterior probability" = post_prob(helping_glmm8, helping_glmm9, helping_glmm10)) %>% 
  mutate("Model" = c("Null", "Reduced", "Full")) %>% 
  relocate("Model",.before = "Posterior probability")
BF_help_latency <- tibble("Bayes factor" = sort(postprob_help_latency$"Posterior probability",T)[1]/sort(postprob_help_latency$"Posterior probability",T)[2])

loo_help8 <- loo(helping_glmm8)
loo_help9 <- loo(helping_glmm9)
loo_help10 <- loo(helping_glmm10)
loo_help_latency <- loo_compare(loo_help8, loo_help9, loo_help10)

##Main Text: Table 4
table_fixef <- rbind(fixef(commit_glmm1),
                     fixef(commit_glmm8),
                     fixef(commit_glmm13),
                     fixef(sharing_glmm1),
                     fixef(helping_glmm1),
                     fixef(helping_glmm8)) %>% 
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
  geom_text(aes(x = 1.5,  y = 1.03, label = asterisk), size = 10) + 
  theme(text = element_text(size = 30),
        axis.title.y = element_text(vjust = 2),
        axis.title.x = element_text(vjust = 0.25))

ggsave("Figure 2.png", height = 14, width = 14)

##Main Text: Figure 3
fig3_leave_taking <- leave_taking %>% 
  mutate("asterisk" = ifelse(age_group_string == "Younger", "*", NA),
         condition_string = fct_recode(condition_string,
                                       "We-framing" = "We",
                                       "You-framing" = "You"))

ggplot(fig3_leave_taking, aes(x = condition_string, fill = factor(leave_taking_verbal_and_nonverbal, levels = c(0,1)))) +
  geom_bar(position = "fill") +
  facet_grid(cols = vars(age_group_string)) +
  labs(x = "Age group", y = "Proportion of participants", fill = "Participant decision") + 
  scale_fill_manual(values = c("navajowhite2","palegreen4"), labels = c("Spontaneously leave", "Take leave"), limits = c("0", "1")) +
  geom_text(stat = "count", aes(label = ..count..), position = "fill", vjust = 1.5,  size = 8) +
  geom_text(aes(x = 1.5,  y = 1.03, label = asterisk), size = 10) + 
  theme(text = element_text(size = 30),
        axis.title.y = element_text(vjust = 2),
        axis.title.x = element_text(vjust = 0.25))

ggsave("Figure 3.png", height = 14, width = 14)

##Supplementary Material
#Commitment: Cox PH models
commit_glmm1_surv <- brm(commitment_latency | cens(censored) ~ age_group_string + location_string + first + gender_string,
                         data = commit_latency_surv,
                         family = cox,
                         prior = prior_survival,
                         iter = 10000,
                         warmup = 2000,
                         control = list(adapt_delta = 0.9999, max_treedepth = 60),
                         save_pars = save_pars(all = TRUE),
                         cores = 4)

commit_glmm2_surv <- brm(commitment_latency | cens(censored) ~ age_group_string + condition_string + location_string + first + gender_string,   
                         data = commit_latency_surv,
                         family = cox,
                         prior = prior_survival,
                         iter = 10000,
                         warmup = 2000,
                         control = list(adapt_delta = 0.9999, max_treedepth = 60),
                         save_pars = save_pars(all = TRUE),
                         cores = 4)

commit_glmm3_surv <- brm(commitment_latency | cens(censored) ~ age_group_string * condition_string + location_string + first + gender_string,
                          data = commit_latency_surv,
                          family = cox,
                          prior = prior_survival,
                          iter = 10000,
                          warmup = 2000,
                          control = list(adapt_delta = 0.9999, max_treedepth = 60),
                          save_pars = save_pars(all = TRUE),
                          cores = 4)

plot(commit_glmm1_surv)
plot(commit_glmm2_surv)
plot(commit_glmm3_surv)

postprob_commit_surv <- tibble("Posterior probability" = post_prob(commit_glmm1_surv, commit_glmm2_surv, commit_glmm3_surv)) %>% 
  mutate("Model" = c("Null", "Reduced", "Full")) %>% 
  relocate("Model",.before = "Posterior probability")
BF_commit_surv <- tibble("Bayes factor" = sort(postprob_commit_surv$"Posterior probability",T)[1]/sort(postprob_commit_surv$"Posterior probability",T)[2])

loo_commit1_surv <- loo(commit_glmm1_surv)
loo_commit2_surv <- loo(commit_glmm2_surv)
loo_commit3_surv <- loo(commit_glmm3_surv)
loo_commit_surv <- loo_compare(loo_commit1_surv, loo_commit2_surv, loo_commit3_surv)

#Helping: Cox PH models
helping_glmm1_surv <- brm(helping_latency | cens(censored) ~ age_group_string + location_string + first + gender_string,
                          data = help_latency_surv,
                          family = cox,
                          prior = prior_survival,
                          iter = 10000,
                          warmup = 2000,
                          control = list(adapt_delta = 0.9999, max_treedepth = 50),
                          save_pars = save_pars(all = TRUE),
                          cores = 4)

helping_glmm2_surv <- brm(helping_latency | cens(censored) ~ age_group_string + condition_string + location_string + first + gender_string,   
                          data = help_latency_surv,
                          family = cox,
                          prior = prior_survival,
                          iter = 10000,
                          warmup = 2000,
                          control = list(adapt_delta = 0.9999, max_treedepth = 50),
                          save_pars = save_pars(all = TRUE),
                          cores = 4)

helping_glmm3_surv <- brm(helping_latency | cens(censored) ~ age_group_string * condition_string + location_string + first + gender_string,
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

loo_help1_surv <- loo(helping_glmm1_surv)
loo_help2_surv <- loo(helping_glmm2_surv)
loo_help3_surv <- loo(helping_glmm3_surv)
loo_help_surv <- loo_compare(loo_help1_surv, loo_help2_surv, loo_help3_surv)

#Supplementary Material Table 1
commit_glmm1_surv_df <- fixef(commit_glmm1_surv) %>% 
  as.data.frame %>% 
  mutate("exp(Estimate)" = exp(Estimate),
         "exp(Error)" = exp(Est.Error),
         "exp(Q2.5)" = exp(Q2.5),
         "exp(Q97.5)" = exp(Q97.5)) %>% 
  select("exp(Estimate)":"exp(Q97.5)") %>% 
  round(digits = 2)

helping_glmm1_surv_df <- fixef(helping_glmm1_surv) %>% 
  as.data.frame %>% 
  mutate("exp(Estimate)" = exp(Estimate),
         "exp(Error)" = exp(Est.Error),
         "exp(Q2.5)" = exp(Q2.5),
         "exp(Q97.5)" = exp(Q97.5)) %>% 
  select("exp(Estimate)":"exp(Q97.5)") %>% 
  round(digits = 2)
