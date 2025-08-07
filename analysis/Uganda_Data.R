
library(tidyverse)

# -------------------------------------------------------- o
# 1. First lets work out suitable EIR and ft ranges --------
# -------------------------------------------------------- o

install.packages("RCurl")

library(RCurl)

url <- "https://raw.githubusercontent.com/OJWatson/hrpup/main/analysis/data_derived/chaokun_uganda.csv"

# simple plots to show what this is:
uganda_dat %>% filter(year>2015) %>% ggplot(aes(year, pfpr210, color = name_1)) + geom_line() +
  theme_bw() + xlab("Year") + ylab("Malaria Prevalence") +
  scale_color_discrete(name = "Admin")

uganda_dat %>% filter(year>2015) %>% ggplot(aes(year, ft, color = name_1)) + geom_line() +
  theme_bw() + xlab("Year") + ylab("Treatment Coverage (ft)") +
  scale_color_discrete(name = "Admin")

# Next we need to get an ft range and an EIR range for these malaraia prevalence
pfpr_to_eir_heuristic <- function(ft, pfpr, sdd = 0.03){

  if(sdd < 0.03) stop("sdd must not be less than 0.03 or likelihood may fail")

  ll_function <- function(params, data) {

    # get eq
    start <- phi_eir_rel(exp(params[["eir"]]), data$ft)

    # what is the prev
    prev <- start$D + start$T + start$A

    # calculate the log likelihood using a truncated normal
    ll <- log(truncnorm::dtruncnorm(prev, a = 0, b = 1, mean = data$pfpr, sd = data$sdd))
    return(ll)
  }

  # starting conditions
  # we use a log to make it easier for the solver
  start <- c("eir" = log(20))

  # create the data list
  data <- list("ft" = ft, "pfpr" = pfpr, "sdd" = sdd)

  # create bounds
  lower <- log(c(0.01))
  upper <- log(c(500))

  # fit our best model
  fit <- optim(
    par = start,
    fn = ll_function,
    data = data,
    method = "L-BFGS-B",
    lower = lower,
    upper = upper,
    control = list(
      trace = TRUE,
      fnscale = -1,
      maxit = 10000
    ),
    hessian = TRUE
  )

  return(exp(as.numeric(fit$par)))

}

# what are the ranges in pfpr and ft for each region
uga_ranges <- uganda_dat %>%
  filter(year > 2015) %>%
  group_by(name_1) %>%
  summarise(pfpr_low = min(pfpr210),
            pfpr_high = max(pfpr210),
            ft_low = min(ft),
            ft_high = max(ft))

install.packages("truncnorm")

# now work out what the eir should be
uga_eir_ranges <- uga_ranges %>%
  split(.$name_1) %>%
  map(.f = function(x){
    x$eir_high <- pfpr_to_eir_heuristic(x$ft_high, x$pfpr_high)
    x$eir_low <- pfpr_to_eir_heuristic(x$ft_low, x$pfpr_low)
    return(x)
  }) %>%
  do.call(rbind, .) %>%
  rename(district = name_1)

# -------------------------------------------------------- o
# 2. Create data for range function --------
# -------------------------------------------------------- o

# Downloading the Uganda allele frequency data using RCurl
url2 <- "https://raw.githubusercontent.com/bailey-lab/selmar/main/analysis/data/data-derived/Uganda_allele_frequency.txt"
uga_data <- read.csv(text = getURL(url2, ssl.verifypeer = FALSE), sep = " ")

# Then continuing with filtering
uga_data <- uga_data %>%
  filter(Locus == "K13") %>%
  select(district = District, year, n, freq)

# work out suitable starting res
start_freq <- function(x){
  if(x[1]>0) {
    return(x[1])
  } else {
    return(mean(x[seq_len(which(x>0)[1])]))
  }
}

# find out starting freq and join with eir, ft data
uga_start_res <- uga_data %>% group_by(district) %>%
  summarise(f1 = start_freq(freq),
            years = max(year)-min(year),
            year0 = year[1])

# not to just fix a few spelling and admin name differences
uga_eir_ranges <- uga_eir_ranges %>%
  mutate(district = replace(district, district == "Kabale", "Rukiga"))
uga_eir_ranges <- uga_eir_ranges %>%
  mutate(district = replace(district, district == "Amolatar", "Amoleta"))

# and merge this data together
uga_param_data <- left_join(uga_eir_ranges %>% mutate(district = replace(district, district == "Amolatar", "Amoleta")), uga_start_res)

