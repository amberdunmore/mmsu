make_dat <- function(eir = 2, ft = 0.5, simulation_years = 30, day0_res=0.1, treatment_failure_rate = 0.43, baseline_ratio = 1.265, treated_ratio = 1.265, asymptomatic_ratio = 1, print = FALSE) {

  model <- malaria_model(
    EIR = eir,
    ft = ft,
    ton = 365,
    toff = 365 + (simulation_years * 365),
    day0_res = day0_res,
    treatment_failure_rate = treatment_failure_rate,
    rT_r_cleared = 0.1,
    rT_r_failed = 0.1,
    resistance_baseline_ratio = baseline_ratio,
    resistance_cleared_ratio = treated_ratio,
    resistance_failed_ratio = treated_ratio,
    resistance_asymptomatic_ratio = asymptomatic_ratio
  )

  times <- seq(0, 365 + (simulation_years * 365), by = 30)
  output <- model$run(times)

  dat <- data.frame(
    time_years = (output[, "t"] - 365) / 365,  # Adjusted to start from resistance introduction
    resistance_prevalence = output[, "prevalence_res"] * 100,
    total_prevalence = output[, "prevalence"] * 100
  )

  if(print) {
    print(ggplot(dat, aes(time_years, resistance_prevalence)) + geom_line())
  }
  return(dat)

}

selfn <- function(out) {

  # Calculate selection coefficient using the actual years of observation
  p0_idx <- which.min(abs(out[, 1] - (365 + 365)))  # ~2 years from start
  p1_idx <- which.min(abs(out[, 1] - (365 + (2*365))))  # ~3 years from start

  p0 <- out$resistance_prevalence[20]/100
  p1 <- out$resistance_prevalence[32]/100

  sel_coeff <- log(p1/(1-p1)) - log(p0/(1-p0))
  return(sel_coeff)
}


make_dat(ft = 0.5, eir = 2, asymptomatic_ratio = 1.265) %>% selfn()
uga_param_data <- readRDS("analysis/data-derived/uga_param_data.rds")

clinical_kb <- uga_param_data %>%
  split(.$district) %>%
  lapply(function(x){
    grid <- expand.grid(
      eir = seq(x$eir_low[1], x$eir_high[1], length.out = 10),
      ft = seq(x$ft_low[1], x$ft_high[1], length.out = 10),
      day0_res = x$f1[1])
    for(i in seq_along(grid$eir)) {
      grid$s <- selfn(make_dat(eir = grid$eir[i], ft = grid$ft[i], day0_res = grid$day0_res[i], asymptomatic_ratio = 1))
    }
    grid$district <- x$district[1]
    grid$kb <- 1
    return(grid)
  })

all_kb <- uga_param_data %>%
  split(.$district) %>%
  lapply(function(x){
    grid <- expand.grid(
      eir = seq(x$eir_low[1], x$eir_high[1], length.out = 10),
      ft = seq(x$ft_low[1], x$ft_high[1], length.out = 10),
      baseline_ratio = c(0.903, 1.265, 1.759),
      asymptomatic_ratio = c(1, 1.265),
      day0_res = x$f1[1])
    for(i in seq_along(grid$eir)) {
      grid$s[i] <- selfn(make_dat(eir = grid$eir[i], ft = grid$ft[i], day0_res = grid$day0_res[i], asymptomatic_ratio = grid$asymptomatic_ratio[i], baseline_ratio = grid$baseline_ratio[i], treated_ratio = grid$baseline_ratio[i], simulation_years = 5))
    }
    grid$district <- x$district[1]
    return(grid)
  })

# now get real data
# 3. Comparison with Meier-Scherling estimates
individual_estimates <- data.frame(
  mutation = c("Pro441Leu", "Cys469Phe", "Cys469Tyr", "Ala675Val"),
  estimate = c(0.494, 0.324, 0.383, 0.237),
  lower = c(-0.462, -0.629, 0.207, 0.087),
  upper = c(1.410, 1.150, 0.591, 0.403),
  type = "Individual",
  order = 4:7
)

uganda_combined <- data.frame(
  mutation = "Uganda Combined",
  estimate = 0.381,
  lower = 0.298,
  upper = 0.472,
  type = "Uganda Combined",
  order = 3
)

model_summary <- do.call(rbind, all_kb)  %>%
  group_by(asymptomatic_ratio) %>%
  summarise(
    mean_sel = mean(s, na.rm = TRUE),
    lower_ci = quantile(s, 0.025, na.rm = TRUE),
    upper_ci = quantile(s, 0.975, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  mutate(
    mutation = paste("Model Combined\n",
                     ifelse(asymptomatic_ratio == 1, "Clinical Only", "All Resistant")),
    type = "Model",
    order = ifelse(asymptomatic_ratio == 1, 1, 2)
  )

meier_data <- rbind(individual_estimates, uganda_combined)

p3 <- ggplot() +
  geom_hline(yintercept = 0) +
  geom_errorbar(data = model_summary,
                aes(x = reorder(mutation, order), ymin = lower_ci, ymax = upper_ci, color = as.factor(asymptomatic_ratio)),
                width = 0.2, size = 1) +
  geom_point(data = model_summary,
             aes(x = reorder(mutation, order), y = mean_sel, color = as.factor(asymptomatic_ratio)),
             size = 3) +
  geom_errorbar(data = meier_data,
                aes(x = reorder(mutation, order), ymin = lower, ymax = upper,
                    color = type),
                width = 0.2, size = 1) +
  geom_point(data = meier_data,
             aes(x = reorder(mutation, order), y = estimate, color = type),
             size = 3) +
  labs(
    x = "Estimate Type",
    y = "Selection coefficient (per year)",
    color = "Estimate Source"
  ) +
  theme_kb_16() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
    legend.position = "bottom",
    axis.line = element_line(), plot.background = element_rect(fill = "white", color = "white")
  ) +
  scale_color_manual(
    values = c(
      "1" = "#E31A1C",
      "1.265" = "#1F78B4",
      "Uganda Combined" = "#FF7F00",
      "Individual" = "#228B22"
    ),
    labels = c(
      "clinical_only" = "Model: Clinical States Only",
      "all_resistant" = "Model: All Resistant States",
      "Uganda Combined" = "Uganda: Combined Estimate",
      "Individual" = "Uganda: Individual Mutations"
    ),
    breaks = c("clinical_only", "all_resistant", "Uganda Combined", "Individual")
  )
p3
ggsave("uganda_meier_scherling_comparison_fixed.png", p3, width = 16, height = 8, dpi = 300)


# and now for the Agago comparison
dat1 <- make_dat(eir = uga_param_data$eir_low[1] + ((uga_param_data$eir_high[1]-uga_param_data$eir_low[1])/2),
                 ft = uga_param_data$ft_low[1] + ((uga_param_data$ft_high[1]-uga_param_data$ft_low[1])/2),
                  day0_res = uga_param_data$f1[1], asymptomatic_ratio = 1, simulation_years = 50)

dat2 <- make_dat(eir = uga_param_data$eir_low[1] + ((uga_param_data$eir_high[1]-uga_param_data$eir_low[1])/2),
                 ft = uga_param_data$ft_low[1] + ((uga_param_data$ft_high[1]-uga_param_data$ft_low[1])/2),
                 day0_res = uga_param_data$f1[1], asymptomatic_ratio = 1.265, simulation_years = 50)

p4 <- rbind(dat1 %>% mutate(scen = "kB Clinical Only"),
      dat2 %>% mutate(scen = "kB All States")) %>%
  ggplot(aes(time_years, resistance_prevalence/100, color = scen)) +
  geom_line() +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  scale_x_continuous(limits = c(0, 50)) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "gray50", alpha = 0.7) +
  labs(
    x = "Time (years)",
    y = "Resistance Prevalence"
  ) +
  theme_kb_16() +
  theme(plot.title = element_text(face = "bold", size = 16, color = "#E31A1C"),
        axis.line = element_line(), plot.background = element_rect(fill = "white", color = "white")) +
  scale_color_discrete(name = "kB Infetious States")
ggsave("kb_comparison_fixed.png", p4, width = 10, height = 6, dpi = 300)


# 36.5 years
(dat1 %>% filter(resistance_prevalence > 95))$time_years[1]
# 19.7 years
(dat2 %>% filter(resistance_prevalence > 95))$time_years[1]
