##Packages and seed
library(readxl)
library(tidyverse)
library(janitor)
library(brms)

set.seed(31)

##Data
wesch_data <- read_xlsx("WESCH Data.xlsx") %>% 
  clean_names() %>% 
  mutate(condition_string = as.factor(condition_string),
         age_group_string = as.factor(age_group_string),
         location_string = as.factor(location_string),
         gender_string = as.factor(gender_string))

wesch_data$age_group_string <- relevel(wesch_data$age_group_string, ref = "Younger")
wesch_data$condition_string <- relevel(wesch_data$condition_string, ref = "You")

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
  filter(approach_n_y != 0)

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
commit_glmm8 <- brm(commitment_latency ~ age_group_string + location_string + first + gender_string +
                      (1|location_string),
                    data = commit_latency,
                    family = lognormal,
                    prior = prior_latencies,
                    iter = 10000,
                    warmup = 2000,
                    control = list(adapt_delta = 0.9999, max_treedepth = 70),
                    save_pars = save_pars(all = TRUE),
                    cores = 4)

commit_glmm9 <- brm(commitment_latency ~ age_group_string + condition_string + location_string + first + gender_string +
                      (1|location_string),
                    data = commit_latency,
                    family = lognormal,
                    prior = prior_latencies,
                    iter = 10000,
                    warmup = 2000,
                    control = list(adapt_delta = 0.9999, max_treedepth = 50),
                    save_pars = save_pars(all = TRUE),
                    cores = 4)

commit_glmm10 <- brm(commitment_latency ~ age_group_string * condition_string + location_string + first + gender_string +
                       (1|location_string),
                     data = commit_latency,
                     family = lognormal,
                     prior = prior_latencies,
                     iter = 10000,
                     warmup = 2000,
                     control = list(adapt_delta = 0.9999, max_treedepth = 70),
                     save_pars = save_pars(all = TRUE),
                     cores = 4)

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

commit_glmm8_surv <- brm(commitment_latency | cens(censored) ~ age_group_string + location_string + first + gender_string,
                         data = commit_latency_surv,
                         family = cox,
                         prior = prior_survival,
                         iter = 10000,
                         warmup = 2000,
                         control = list(adapt_delta = 0.9999, max_treedepth = 60),
                         save_pars = save_pars(all = TRUE),
                         cores = 4)

commit_glmm9_surv <- brm(commitment_latency | cens(censored) ~ age_group_string + condition_string + location_string + first + gender_string,   
                         data = commit_latency_surv,
                         family = cox,
                         prior = prior_survival,
                         iter = 10000,
                         warmup = 2000,
                         control = list(adapt_delta = 0.9999, max_treedepth = 60),
                         save_pars = save_pars(all = TRUE),
                         cores = 4)

commit_glmm10_surv <- brm(commitment_latency | cens(censored) ~ age_group_string * condition_string + location_string + first + gender_string,
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

postprob_commit_surv <- tibble("Posterior probability" = post_prob(commit_glmm8_surv, commit_glmm9_surv, commit_glmm10_surv)) %>% 
  mutate("Model" = c("Null", "Reduced", "Full")) %>% 
  relocate("Model",.before = "Posterior probability")
BF_commit_surv <- tibble("Bayes factor" = sort(postprob_commit_surv$"Posterior probability",T)[1]/sort(postprob_commit_surv$"Posterior probability",T)[2])

loo_commit8_surv <- loo(commit_glmm8_surv)
loo_commit9_surv <- loo(commit_glmm9_surv)
loo_commit10_surv <- loo(commit_glmm10_surv)
loo_commit_surv <- loo_compare(loo_commit8_surv, loo_commit9_surv, loo_commit10_surv)

#Leave-taking
commit_glmm11 <- brm(leave_taking_type_recoded ~ age_group_string + location_string + first + gender_string +
                       (1|location_string),
                     data = leave_taking,
                     family = bernoulli,
                     prior = prior_tendencies,
                     iter = 10000,
                     warmup = 2000,
                     control = list(adapt_delta = 0.9999, max_treedepth = 50),
                     save_pars = save_pars(all = TRUE),
                     cores = 4)

commit_glmm12 <- brm(leave_taking_type_recoded ~ age_group_string + condition_string + location_string + first + gender_string +
                       (1|location_string),
                     data = leave_taking,
                     family = bernoulli,
                     prior = prior_tendencies,
                     iter = 10000,
                     warmup = 2000,
                     control = list(adapt_delta = 0.9999, max_treedepth = 50),
                     save_pars = save_pars(all = TRUE),
                     cores = 4)

commit_glmm13 <- brm(leave_taking_type_recoded ~ age_group_string * condition_string + location_string + first + gender_string +
                       (1|location_string),
                     data = leave_taking,
                     family = bernoulli,
                     prior = prior_tendencies,
                     iter = 10000,
                     warmup = 2000,
                     control = list(adapt_delta = 0.9999, max_treedepth = 50),
                     save_pars = save_pars(all = TRUE),
                     cores = 4)

plot(commit_glmm11)
plot(commit_glmm12)
plot(commit_glmm13)

pp_check(commit_glmm11, nsamples = 500)
pp_check(commit_glmm12, nsamples = 500)
pp_check(commit_glmm13, nsamples = 500)

postprob_commit_takeleave <- tibble("Posterior probability" = post_prob(commit_glmm11, commit_glmm12, commit_glmm13)) %>% 
  mutate("Model" = c("Null", "Reduced", "Full")) %>% 
  relocate("Model",.before = "Posterior probability")
BF_commit_takeleave <- tibble("Bayes factor" = sort(postprob_commit_takeleave$"Posterior probability",T)[1]/sort(postprob_commit_takeleave$"Posterior probability",T)[2])

loo_commit11 <- loo(commit_glmm11)
loo_commit12 <- loo(commit_glmm12)
loo_commit13 <- loo(commit_glmm13)
loo_commit_takeleave <- loo_compare(loo_commit11, loo_commit12, loo_commit13)

hypothesis(commit_glmm13, "condition_stringWe > 0")

posterior_samples_takeleave <- posterior_samples(commit_glmm13)
HPD_takeleave_diff <- mean(posterior_samples_takeleave$b_condition_stringWe < 0) - 
  mean(posterior_samples_takeleave$b_condition_stringWe < -0.08)
HPD_takeleave_ratio <- (HPD_takeleave_diff/.95)

table(leave_taking$first, leave_taking$leave_taking_type_recoded)

leave_taking_noorder3 <- leave_taking %>% 
  filter(first != 5)

commit_glmm13_noorder3 <- brm(leave_taking_type_recoded ~ age_group_string * condition_string + location_string + first + gender_string +
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
commit_glmm14 <- brm(leave_taking_type_recoded ~ age_group_string +
                       (1|location_string),
                    data = commit_leavetake_we,
                    family = bernoulli,
                    prior = prior_tendencies,
                    iter = 10000,
                    warmup = 2000,
                    control = list(adapt_delta = 0.9999, max_treedepth = 50),
                    save_pars = save_pars(all = TRUE),
                    cores = 4)

plot(commit_glmm14)
pp_check(commit_glmm14, nsamples = 100)

#Leave-taking tendency: Looking at age within you-framing
commit_glmm15 <- brm(leave_taking_type_recoded ~ age_group_string +
                       (1|location_string),
                    data = commit_leavetake_you,
                    family = bernoulli,
                    prior = prior_tendencies,
                    iter = 10000,
                    warmup = 2000,
                    control = list(adapt_delta = 0.9999, max_treedepth = 50),
                    save_pars = save_pars(all = TRUE),
                    cores = 4)

plot(commit_glmm15)
pp_check(commit_glmm15, nsamples = 100)

#Leave-taking tendency: Looking at condition within age [2.5,3.5]
commit_glmm16 <- brm(leave_taking_type_recoded ~ condition_string +
                       (1|location_string),
                    data = commit_leavetake_3yo,
                    family = bernoulli,
                    prior = prior_tendencies,
                    iter = 10000,
                    warmup = 2000,
                    control = list(adapt_delta = 0.9999, max_treedepth = 50),
                    save_pars = save_pars(all = TRUE),
                    cores = 4)

plot(commit_glmm16)
pp_check(commit_glmm16, nsamples = 100)

hypothesis(commit_glmm16, "condition_stringWe > 0")

#Leave-taking tendency: Looking at condition within age [3.5,4.5]
commit_glmm17 <- brm(leave_taking_type_recoded ~ condition_string +
                       (1|location_string),
                    data = commit_leavetake_4yo,
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
  filter(number_of_erasers_shared < 7)

#Resource distribution
sharing_glmm1 <- brm(number_of_erasers_shared_recoded ~ age_group_string + location_string + first + gender_string +
                       (1|location_string),
                     data = sharing,
                     family = bernoulli,
                     prior = prior_tendencies,
                     iter = 10000,
                     warmup = 2000,
                     control = list(adapt_delta = 0.9999, max_treedepth = 50),
                     save_pars = save_pars(all = TRUE),
                     cores = 4)

sharing_glmm2 <- brm(number_of_erasers_shared_recoded ~ age_group_string + condition_string + location_string + first + gender_string +
                       (1|location_string),
                     data = sharing,
                     family = bernoulli,
                     prior = prior_tendencies,
                     iter = 10000,
                     warmup = 2000,
                     control = list(adapt_delta = 0.9999, max_treedepth = 50),
                     save_pars = save_pars(all = TRUE),
                     cores = 4)

sharing_glmm3 <- brm(number_of_erasers_shared_recoded ~ age_group_string * condition_string + location_string + first + gender_string +
                       (1|location_string),
                     data = sharing,
                     family = bernoulli,
                     prior = prior_tendencies,
                     iter = 10000,
                     warmup = 2000,
                     control = list(adapt_delta = 0.9999, max_treedepth = 50),
                     save_pars = save_pars(all = TRUE),
                     cores = 4)

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
                     control = list(adapt_delta = 0.9999, max_treedepth = 50),
                     save_pars = save_pars(all = TRUE),
                     cores = 4)

helping_glmm9 <- brm(helping_latency ~ age_group_string + condition_string + location_string + first + gender_string +
                       (1|location_string),
                     data = help_latency,
                     family = lognormal,
                     prior = prior_latencies,
                     iter = 10000,
                     warmup = 2000,
                     control = list(adapt_delta = 0.9999, max_treedepth = 50),
                     save_pars = save_pars(all = TRUE),
                     cores = 4)

helping_glmm10 <- brm(helping_latency ~ age_group_string * condition_string + location_string + first + gender_string +
                        (1|location_string),
                      data = help_latency,
                      family = lognormal,
                      prior = prior_latencies,
                      iter = 10000,
                      warmup = 2000,
                      control = list(adapt_delta = 0.9999, max_treedepth = 50),
                      save_pars = save_pars(all = TRUE),
                      cores = 4)

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

helping_glmm8_surv <- brm(helping_latency | cens(censored) ~ age_group_string + location_string + first + gender_string,
                          data = help_latency_surv,
                          family = cox,
                          prior = prior_survival,
                          iter = 10000,
                          warmup = 2000,
                          control = list(adapt_delta = 0.9999, max_treedepth = 60),
                          save_pars = save_pars(all = TRUE),
                          cores = 4)

helping_glmm9_surv <- brm(helping_latency | cens(censored) ~ age_group_string + condition_string + location_string + first + gender_string,   
                          data = help_latency_surv,
                          family = cox,
                          prior = prior_survival,
                          iter = 10000,
                          warmup = 2000,
                          control = list(adapt_delta = 0.9999, max_treedepth = 60),
                          save_pars = save_pars(all = TRUE),
                          cores = 4)

helping_glmm10_surv <- brm(helping_latency | cens(censored) ~ age_group_string * condition_string + location_string + first + gender_string,
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

postprob_help_surv <- tibble("Posterior probability" = post_prob(helping_glmm8_surv, helping_glmm9_surv, helping_glmm10_surv)) %>% 
  mutate("Model" = c("Null", "Reduced", "Full")) %>% 
  relocate("Model",.before = "Posterior probability")
BF_help_surv <- tibble("Bayes factor" = sort(postprob_help_surv$"Posterior probability",T)[1]/sort(postprob_help_surv$"Posterior probability",T)[2])

loo_help8_surv <- loo(helping_glmm8_surv)
loo_help9_surv <- loo(helping_glmm9_surv)
loo_help10_surv <- loo(helping_glmm10_surv)
loo_help_surv <- loo_compare(loo_help8_surv, loo_help9_surv, loo_help10_surv)

##Table 4
table_fixef <- rbind(fixef(commit_glmm1),
                     fixef(commit_glmm8),
                     fixef(commit_glmm13),
                     fixef(sharing_glmm1),
                     fixef(helping_glmm1),
                     fixef(helping_glmm8)) %>% 
  round(digits = 2)
write.csv(table_fixef, "Table 4.csv")

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
         "exp(Q97.5)" = exp(Q97.5)) %>% 
  select("exp(Estimate)":"exp(Q97.5)") %>% 
  round(digits = 2)

helping_glmm8_surv_df <- fixef(helping_glmm8_surv) %>% 
  as.data.frame %>% 
  mutate("exp(Estimate)" = exp(Estimate),
         "exp(Error)" = exp(Est.Error),
         "exp(Q2.5)" = exp(Q2.5),
         "exp(Q97.5)" = exp(Q97.5)) %>% 
  select("exp(Estimate)":"exp(Q97.5)") %>% 
  round(digits = 2)