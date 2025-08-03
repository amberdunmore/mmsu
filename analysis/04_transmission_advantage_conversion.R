library(tidyverse)

# Step 1: Get ratio increase in gametocyte density between subpatent and patent -------------

## Read in data
dat <- lapply(1:7, function(x) {
  readxl::read_xlsx(path = "analysis/data-raw/slater_natcomms_2019_fig_6.xlsx", sheet = x)
  }) %>% do.call(rbind, .) %>%
  janitor::clean_names()

## Get study groups
dat <- dat %>% separate_wider_delim(1, "_", names = c("G", "study", "id"))

# get summary statistics as subpatent and patent gams not shared in data
dat <- dat %>% group_by(study) %>%
  mutate(N = n()) %>%
  group_by(study, microscoscopy_status) %>%
  summarise(n = n(),
    perc = n()/unique(N),
            mean = mean(gametocyte_density),
            q025 = quantile(gametocyte_density, 0.025),
            q25 = quantile(gametocyte_density, 0.25),
            q50 = median(gametocyte_density),
            q75 = quantile(gametocyte_density, 0.75),
            q975 = quantile(gametocyte_density, 0.975)
            ) %>% ungroup() %>%
  rename(microscopy_status = microscoscopy_status)

## Patent data from paper
patdf <- data.frame(
  study = dat$study,
  microscopy_status = dat$microscopy_status,
  pcr_gams_x = c(105, 142, 57, 232, 92, 102, 9, 4,
                     65, 41, 143, 90, 21, 3),
  pcr_gams_n = c(242, 154, 123, 289, 165, 109, 876, 6,
                  357, 56, 1966, 117, 213, 13),
  pfpr = c(66.5, 66.5, 82.5, 82.5, 90.6, 90.6, 9.7, 9.7, 38, 38, 18.5, 18.5, 38.1, 38.1)
)

patdf <- patdf %>%
  mutate(pcr_gams = pcr_gams_x / pcr_gams_n) %>%
  mutate(gams_n = pcr_gams_n)

## Make data frame
df <- left_join(dat, patdf, by = c("study", "microscopy_status"))
df <- df %>% group_by(study) %>%
  summarise(across(mean:q975, .fns = ~ log(.x[.data$microscopy_status == 1]/.x[.data$microscopy_status == 0])),
           pcr_gams_x = sum(pcr_gams_x),
           gams_n = sum(pcr_gams_n),
           n = min(n))

# make plot to
gg1 <- df %>%
  ggplot(aes(qlogis(pcr_gams_x/gams_n), (q50))) +
  geom_point(aes(size = n)) +
  geom_smooth(aes(weight = gams_n), method = "lm") +
  theme_bw(base_size = 14, base_family = "Helvetica") +
  xlab("Log Odds Gametocyte Prevalence\n") +
  ylab("\nLog Gametocyte Density (p/microL)") +
  scale_size_binned("Sample Size\n", guide = guide_bins(show.limits = TRUE)) +
  theme(legend.key.size = unit(2, "line"))
gg1
fit <- lm(q50 ~ qlogis(pcr_gams_x/gams_n), data = df, weights = gams_n)
grad <- as.numeric(fit$coefficients[2])
grad

# how we go from gam prevalence to fold change in gametocyte density
exp(grad*(qlogis(0.45) - qlogis(0.15)))


# Step 2: Translate this for our model

# Our data comes from symptomatic patients who have cD probability of infecting mosquiotes:
cD <- ICDMM::model_param_list_create()$cD

# cD is 0.067, which based on various literature corresponds to the following gam density:
#
# 260 (100 - NA) (https://elifesciences.org/articles/00626#s2 - Churcher at al elife 2013) - central only as higher is off plot
# 1.5 (0.2-10) (https://journals.plos.org/ploscompbiol/article/figure?id=10.1371/journal.pcbi.1003025.g003 - Johnston at al Plos Comp Bio 2013)
# 100 (15 - 1100) (Fig3a https://malariajournal.biomedcentral.com/articles/10.1186/s12936-016-1584-z - Pett et al Malaria Jounral 2016)
# 95 (12 - 2800) (Fig3b https://malariajournal.biomedcentral.com/articles/10.1186/s12936-016-1584-z - Pett et al Malaria Jounral 2016)
# 100 (20 - 990) (Fig3c https://malariajournal.biomedcentral.com/articles/10.1186/s12936-016-1584-z - Pett et al Malaria Jounral 2016)
# 30 (20 - 100) (Fig1 https://academic.oup.com/trstmh/article/115/12/1462/6295667 - Ahmad et al Transaction RSTMH 2021)

# Taking just a simple central mean gives ~100
mean_gams <- mean(c(260,1.5, 100, 95, 100, 30))

# however, they all have very different relationships between density and infection
# and our data will be anywhere between 1 and 3x game density based on gam prevalence data
# realistically. Based on this we will ignore Churcher et al as 3x 250 goes off scale
# (though suggests flattening to max 0.2 - so 3x max increase in cD)
#
# Johnstan et al are around 10-20 probability - 1.5 - 3x max
# Pett et al are all around 9-12.5-18 (except for microscopy which is very different and likely innaccurate)
# Ahmad et al are all around 9-13-20 as well
# so likely a consistent increase observed for the ranges we are seeing

# fit logistic to this general trend:
logistic <- function(x, L, k, x0) {
  L / (1 + exp(-k * (x - x0)))
}


# Observed data
x <- c(0.01, 100, 200, 300)
y <- c(0.04, 0.067, 0.11, 0.15)

# Fit the model using nls
fit <- nls(
  y ~ logistic(x, L, k, x0),
  start = list(L = 0.3, k = 0.02, x0 = 150),
  algorithm = "port",
  lower = c(0.001, 0.0001, -Inf),
  upper = c(1, 1, Inf)
)

# View estimated parameters
summary(fit)

# Predict and plot
xnew <- exp(seq(-6, 6, length.out = 100))
plot(xnew, predict(fit, newdata = data.frame(x = xnew)), type = "l", ylim = c(0, 0.3))
points(x, y, col = "blue")

# Great let's create this function and save it be used to convert with
c_multiplier <- function(transmission_multiplier) {

  # gradient from out gam fitting
  # and 0.1 gam prevalence at 0% mutation
  gam_fold <- exp(0.2342576*(qlogis(0.1*transmission_multiplier) - qlogis(0.1)))

  # new gam density
  new_gam <- 100*gam_fold

  # Custom logistic that starts at zero model fit parameters
  L <- 2.304505e-01
  k <- 7.412884e-03
  x0 <- 2.151131e+02

  # what is out multiplier
  mult <- (L / (1 + exp(-k * (new_gam - x0))))

  # and normalise against 1
  gam_fold <- exp(0.2342576*(qlogis(0.1) - qlogis(0.1)))
  new_gam <- 100*gam_fold
  norm <-  (L / (1 + exp(-k * (new_gam - x0))))

  return(mult/norm)

}
saveRDS(c_multiplier, "analysis/data-derived/c_mult_func.rds")
usethis::use_data(c_multiplier, overwrite = TRUE)
