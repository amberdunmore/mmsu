# analysis/resistance_sweep.R - Explore impact of transmission/duration advantages on resistance dynamics

# Setup environment and libraries
devtools::load_all()  # load the mmsu package (assumes package is installed with above changes)
library(ggplot2)      # for plotting
library(future.apply) # for parallel simulation

# Define parameter grids for exploration
EIR_values <- c(10, 50)                         # low and high transmission intensities (annual EIR)
ft_values  <- c(0.3, 0.7)                       # low and high treatment coverages
resistance_trans_mult_values <- c(1, 1.5, 2, 2.5, 3) # multipliers for transmission advantage (1 = no advantage)
resistance_dur_mult_values  <- c(1, 1.05, 1.1, 1.15)  # multipliers for duration advantage (1 = no advantage)

# set when res is turned on and off
# Run simulation up to 1 year after toff (just to check it stops) days
ton <- 365
toff <- 365 + (20*365)
times <- seq(0, toff+365, by = 1)

# Pre-compute baseline rT_s (treatment clearance rate for sensitive infections) for each EIR-ft combination
base_rT_s <- list()
for (EIR in EIR_values) {
  for (ft in ft_values) {
    params_eq <- phi_eir_rel(EIR, ft)
    key <- paste(EIR, ft)
    base_rT_s[[key]] <- params_eq$rT_s[1]   # use equilibrium value of rT_s for baseline
  }
}

# Prepare a grid of all combinations
param_grid <- expand.grid(EIR = EIR_values,
                          ft  = ft_values,
                          resistance_trans_mult = resistance_trans_mult_values,
                          resistance_dur_mult  = resistance_dur_mult_values)

# Parallelize simulations over the grid
future::plan(future::multisession)  # use multiple cores for parallel loop
results_list <- future_lapply(1:nrow(param_grid), function(i) {

  # Extract combination parameters
  combo <- param_grid[i, ]
  EIR  <- combo$EIR
  ft   <- combo$ft
  res_trans <- combo$resistance_trans_mult
  res_dur   <- combo$resistance_dur_mult

  # Ensure baseline rT_r_true equals the sensitive clearance rate for a fair comparison
  rT_baseline <- base_rT_s[[paste(EIR, ft)]]

  # Initialize model with given parameters
  model <- malaria_model(EIR = EIR, ft = ft, rT_r_true = rT_baseline,
                         resistance_trans_mult = res_trans,
                         resistance_dur_mult  = res_dur,
                         ton = ton, toff = toff)
  output <- model$run(times)

  # Compute final resistance prevalence (fraction of infections that are resistant at end)
  final_prev_res <- output[nrow(output), "prevalence_res"]

  # Estimate selection coefficient from initial growth of resistance frequency
  # (use logit difference between 1 day post-introduction and one year later)
  p0 <- output[output[ , "t"] == ton+1, "prevalence_res"]   # resistant frequency just after introduction (at day 3001)
  p1 <- output[output[ , "t"] == toff, "prevalence_res"]   # resistant frequency one year later (day 3366)
  # Guard against extreme values for logit calculation
  p0 <- max(min(p0, 0.99999999), 1e-8)
  p1 <- max(min(p1, 0.99999999), 1e-8)
  # Calculate selection coefficient (per day) and scale to per year
  sel_rate_day <- (log(p1/(1 - p1)) - log(p0/(1 - p0))) / (toff - ton+1)
  sel_coeff_year <- sel_rate_day * 365  # per-year growth rate of resistant fraction

  # Return results for this combination
  return(data.frame(EIR = EIR, ft = ft,
                    resistance_trans_mult = res_trans,
                    resistance_dur_mult = res_dur,
                    final_resistance_prevalence = final_prev_res,
                    selection_coefficient = sel_coeff_year))
})
# Combine all results into one data frame
results_df <- do.call(rbind, results_list)

# Run illustrative time-series simulations for resistance dynamics (with vs without advantages)
EIR_ts <- 50  # example scenario: high transmission
ft_ts  <- 0.5 # example scenario: moderate treatment coverage

# Ensure baseline rT_r_true equals sensitive clearance rate for this scenario
params_base <- phi_eir_rel(EIR_ts, ft_ts)
rT_base <- params_base$rT_s[1]

# Define scenarios: no advantage vs transmission/duration advantages
scenario_params <- list(
  "No advantage"           = list(res_trans_mult = 1,   res_dur_mult = 1),
  "Transmission advantage" = list(res_trans_mult = 2, res_dur_mult = 1),
  "Duration advantage"     = list(res_trans_mult = 1,   res_dur_mult = 1.1),
  "Both advantages"        = list(res_trans_mult = 2, res_dur_mult = 1.1)
)

# Simulate each scenario
ts_results <- list()
for (scenario in names(scenario_params)) {
  pars <- scenario_params[[scenario]]
  model <- malaria_model(EIR = EIR_ts, ft = ft_ts, rT_r_true = rT_base,
                         resistance_trans_mult = pars$res_trans_mult,
                         resistance_dur_mult  = pars$res_dur_mult,
                         ton = ton, toff = toff)
  out <- model$run(times)
  ts_results[[scenario]] <- data.frame(time = out[ , "t"],
                                       prevalence_res = out[ , "prevalence_res"],
                                       scenario = scenario)
}
ts_df <- do.call(rbind, ts_results)

# Create time-series plot comparing resistance dynamics
p_ts <- ggplot(ts_df, aes(x = time, y = prevalence_res, color = scenario)) +
  geom_line(size = 1) +
  geom_vline(xintercept = params_base$res_time[1], linetype = "dashed", color = "black") +
  labs(x = "Time (days)", y = "Resistant infection prevalence", color = "Scenario",
       title = "Resistance Dynamics with and without Advantages",
       subtitle = paste("EIR =", EIR_ts, ", ft =", ft_ts, "- resistance introduced at day", model$contents()$ton)) +
  theme_minimal()
p_ts

# Prepare data for heatmaps (convert EIR and ft to factors for faceting)
results_df$EIR <- factor(results_df$EIR)
results_df$ft  <- factor(results_df$ft)

# Heatmap of final resistance prevalence
p_prev <- ggplot(results_df, aes(x = resistance_trans_mult, y = resistance_dur_mult, fill = final_resistance_prevalence)) +
  geom_tile() +
  facet_grid(EIR ~ ft, labeller = label_both) +
  scale_fill_viridis_c(name = "Final resistant\nprevalence", limits = c(0, 1)) +
  labs(x = "Transmission multiplier", y = "Duration multiplier",
       title = "Final resistance prevalence",
       subtitle = "Proportion of infections with resistant parasites at 10,000 days") +
  theme_minimal()
p_prev

# Heatmap of selection coefficients
p_sel <- ggplot(results_df, aes(x = resistance_trans_mult, y = resistance_dur_mult, fill = selection_coefficient)) +
  geom_tile() +
  facet_grid(EIR ~ ft, labeller = label_both) +
  scale_fill_viridis_c(name = "Selection coeff.\n(per year)") +
  labs(x = "Transmission multiplier", y = "Duration multiplier",
       title = "Selection coefficient",
       subtitle = "Annual growth rate of resistant parasite frequency after introduction") +
  theme_minimal()
p_sel

# Save plots to files in analysis/figures
if (!dir.exists("analysis/figures")) {
  dir.create("analysis/figures", recursive = TRUE)
}
ggsave("analysis/figures/resistance_dynamics_timeseries.png", p_ts, width = 8, height = 6)
ggsave("analysis/figures/final_resistance_prevalence_heatmap.png", p_prev, width = 8, height = 6)
ggsave("analysis/figures/selection_coefficient_heatmap.png", p_sel, width = 8, height = 6)
