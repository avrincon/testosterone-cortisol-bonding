# Title: Testosterone and cortisol are negatively associated with ritualized bonding behavior in male macaques
# Authors: Alan V. Rincon, Michael Heistermann, Oliver Schülke, Julia Ostner
# Analysis script for testosterone and cortisol models

# options for Stan to run faster
rstan::rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

library(tidyverse)
library(brms)
library(bayesplot)

# read data
d6h <- read_tsv("data-6h-excretion-window.txt")
d13h <- read_tsv("data-13.5h-excretion-window.txt")

# put data in a list for easier data wrangling
d <- list(d6h = d6h,
          d13h = d13h)

# exclude urine samples collected with salivettes (for testosterone analysis)
# and append to list of data frames
d.excl.salivette <- 
  d %>% 
  map(~.x %>% filter(collection_method == "pipette"))

d <- 
  d %>% 
  append(d.excl.salivette) %>% 
  set_names(c("d6h", "d13h", "d6h.excl.sl", "d13h.excl.sl"))


# prep data ---------------------------------------------------------------

# log transform testosterone and cortisol
# keep only samples with > 2hrs of focal observation prior to sample collection
d2 <- 
  d %>% 
  map(~.x %>% 
        filter(sample_window_duration_h > 2) %>%
        mutate(log_cortisol = log(cortisol_ng_mg_crea),
               log_testosterone = log(testosterone_ng_mg_crea),
               infant_care_m_per_h = infant_care_s_per_h/60,
               grooming_m_per_h = grooming_s_per_h/60) %>% 
        mutate_at(vars(triadic_interactions_per_h, infant_care_m_per_h,
                       grooming_m_per_h, aggression_mm_per_h),
                  list(log = ~log(. + 1))))

# within subjects centering
# gather behaviours into one column to calculate mean per subject
d3 <- 
  d2 %>% 
  map(~.x %>% 
        pivot_longer(cols = c("triadic_interactions_per_h_log",
                              "infant_care_m_per_h_log"),
                     names_to = "test_predictors", 
                     values_to = "value") %>% 
        group_by(male_id, test_predictors) %>% 
        mutate(mean = mean(value)) %>% 
        ungroup() %>% 
        # center by subtracting mean from value
        mutate(within = value - mean) %>% 
        pivot_longer(cols = c("value", "mean", "within"), 
                     names_to = "name") %>% 
        # spread behaviours back into columns
        unite(test_predictors_2, test_predictors, name, sep = "_") %>% 
        pivot_wider(names_from = test_predictors_2, values_from = value))

# re-sclae predictors to mean 0, SD 1
d4 <- 
  d3 %>% 
  map(~.x %>% 
        mutate_at(
          vars(triadic_interactions_per_h_log_mean, 
               triadic_interactions_per_h_log_within,
               infant_care_m_per_h_log_mean,
               infant_care_m_per_h_log_within,
               sample_collection_time_s, david_score_norm,
               aggression_mm_per_h_log, grooming_m_per_h_log), 
          ~as.vector(scale(.))))

# fit models --------------------------------------------------------------


# testosterone model, using a 6 hour excretion window
# includes between-subjects effect for infant care but not triadic interactions
# excluding samples collected with salivettes
# Table 2 (model 1)
m1 <-
  brm(data = d4$d6h.excl.sl, family = gaussian,
      log_testosterone ~ 1 + triadic_interactions_per_h_log_mean + 
        triadic_interactions_per_h_log_within + 
        infant_care_m_per_h_log_within +
        grooming_m_per_h_log + aggression_mm_per_h_log + david_score_norm +
        sample_collection_time_s +
        (1 + triadic_interactions_per_h_log_within + 
           infant_care_m_per_h_log_within +
           grooming_m_per_h_log +
           aggression_mm_per_h_log + sample_collection_time_s | male_id),
      prior = c(set_prior("normal(0, 1)", class = "Intercept"),
                set_prior("normal(0, 1)", class = "b"),
                set_prior("cauchy(0, 1)", class = "sd"),
                set_prior("cauchy(0, 1)", class = "sigma"),
                set_prior("lkj(2)", class = "cor")),
      iter = 4000, warmup = 1000, chains = 4, cores = 4,
      control = list(adapt_delta = 0.99),
      sample_prior = "yes", seed = 20)

# cortisol model, using 13.5 hour excretion window
# includes between-subjects effect for infant care but not triadic interactions
# includes all samples
# Table 3 (model 2)
m2 <- 
  update(m1,
         log_cortisol ~ ., 
         newdata = d4$d13h)

# complementary testosterone model, using 6 hour hormone excretion window
# includes between-subjects effect for triadic interactions but not infant care
# Table S3
mS3 <- 
  update(m1, 
         ~ 1 + triadic_interactions_per_h_log_within + 
           infant_care_m_per_h_log_within + infant_care_m_per_h_log_mean +
           grooming_m_per_h_log + aggression_mm_per_h_log + david_score_norm +
           sample_collection_time_s +
           (1 + triadic_interactions_per_h_log_within + 
              infant_care_m_per_h_log_within +
              grooming_m_per_h_log +
              aggression_mm_per_h_log + sample_collection_time_s | male_id),
         newdata = d4$d6h.excl.sl)

# complementary cortisol model, using 13.5 hour hormone excretion window
# includes between-subjects effect for triadic interactions but not infant care
# Table S4
mS4 <- 
  update(mS3,
         log_cortisol ~ ., 
         newdata = d4$d13h)


# models for testosterone using 13.5 hour excretion window (excluding salivette samples)
# Table S1a and Table S1b
mS1a <- update(m1, newdata = d4$d13h.excl.sl)
mS1b <- update(mS3, newdata = d4$d13h.excl.sl)

# models for cortisol excluding salivette samples
# Table S2a and Table S2b
mS2a <- update(m2, newdata = d4$d13h.excl.sl)
mS2b <- update(mS4, newdata = d4$d13h.excl.sl)

# put models in a list for easier processing
m <- 
  list(
    m1 = m1,
    mS3 = mS3,
    m2 = m2,
    mS4 = mS4,
    mS1a = mS1a,
    mS1b = mS1b,
    mS2a = mS2a,
    mS2b = mS2b
  )

# model summaries ---------------------------------------------------------


# veiw coefficients
m %>% map(~fixef(.x) %>% round(digits = 2))

# get posterior samples
post <- 
  m %>% 
  map(~ .x %>% posterior_samples(add_chain = T))

# get posterior beta coefficients
post.betas <- 
  post %>% 
  map(~ .x %>% 
        as_tibble %>% 
        # select beta coefs
        select(starts_with("b_")) %>% 
        gather(key = "predictor", value = "beta_coef") %>%
        mutate(predictor = str_remove(predictor, "b_"))
        )

# calculate proportion of posterior sammples that fall on same side of 0 as the mean (Pr)
post.pr <- 
  post.betas %>% 
  map(~ .x %>% 
        group_by(predictor) %>% 
        summarise(mean = mean(beta_coef),
                  sd = sd(beta_coef),
                  lower_95_CI = quantile(beta_coef, probs = 0.025),
                  upper_95_CI = quantile(beta_coef, probs = 0.975),
                  N_post_samples = n(),
                  N_above_0 = sum(beta_coef > 0),
                  perc_above_0 = N_above_0/N_post_samples*100,
                  perc_below_0 = 100-perc_above_0) %>% 
        mutate(Pr = ifelse(mean > 0, perc_above_0/100, perc_below_0/100)))


# model diagnostics -------------------------------------------------------

# compute number of effective samples
ratios_cp <- 
  m %>% 
  map(~ .x %>% 
        neff_ratio() %>% 
        mcmc_neff(.x, size = 2))


# plot of chains
post.chains <- 
  post %>% 
  map(~.x %>% 
        # select variables from which you want the chains
        select(starts_with("b_"), starts_with("sd_"), sigma, chain) %>%
        # select(starts_with("r_"), chain) %>%
        mcmc_trace(facet_args = list(ncol = 4), 
                   size = .15) +
        labs(title = "Chains trace: beta coefs") +
        theme(legend.position = c(.95, .2))
  )

# check autocorrelation (plots)
post.autocorr <- 
  post %>% 
  map(~.x %>% 
        select(starts_with("b_"), starts_with("sd_"), sigma) %>% 
        mcmc_acf(lags = 5))
