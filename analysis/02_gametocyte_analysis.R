
devtools::load_all()

library(mmsu)
library(ggplot2)
library(dplyr)
library(tidyr)
library(viridis)
library(patchwork)



# EXPERIMENTAL PARAMETER RANGES -------------------------------------------

define_experimental_ranges <- function() {
  estimates <- list(
    in_vivo_baseline = list(
      conservative = 2.08,
      realistic = 2.9,
      optimistic = 3.98
    ),
    in_vivo_treated = list(
      conservative = 0.39,
      realistic = 0.89,
      optimistic = 1.98
    ),
    in_vitro_baseline = list(
      conservative = 0.8, # Just set this myself, not provided by study
      realistic = 1.0,   # Proellochs 2024
      optimistic = 1.0
    ),
    in_vitro_treated = list(
      conservative = 7.0,
      realistic = 7.5, # Rajapandi 2019
      optimistic = 8.0
    )
  )
  return(estimates)
}


# CLASSIFICATION FUNCTION -------------------------------------------------


classify_parameter_combination <- function(baseline_r, treated_r, estimates) {
  tol <- 0.05

  # In vivo scenarios
  if (abs(baseline_r - estimates$in_vivo_baseline$realistic) < tol &&
      abs(treated_r - estimates$in_vivo_treated$realistic) < tol) {
    return("In vivo (central)")
  }

  if (abs(baseline_r - estimates$in_vivo_baseline$conservative) < tol &&
      abs(treated_r - estimates$in_vivo_treated$conservative) < tol) {
    return("In vivo (conservative)")
  }

  if (abs(baseline_r - estimates$in_vivo_baseline$optimistic) < tol &&
      abs(treated_r - estimates$in_vivo_treated$optimistic) < tol) {
    return("In vivo (optimistic)")
  }

  # In vitro scenarios
  if (abs(baseline_r - estimates$in_vitro_baseline$realistic) < tol &&
      abs(treated_r - estimates$in_vitro_treated$realistic) < tol) {
    return("In vitro (central)")
  }

  if (abs(baseline_r - estimates$in_vitro_baseline$conservative) < tol &&
      abs(treated_r - estimates$in_vitro_treated$conservative) < tol) {
    return("In vitro (conservative)")
  }

  if (abs(baseline_r - estimates$in_vitro_baseline$optimistic) < tol &&
      abs(treated_r - estimates$in_vitro_treated$optimistic) < tol) {
    return("In vitro (optimistic)")
  }

  # Check if within experimental ranges
  in_vivo_baseline_range <- baseline_r >= estimates$in_vivo_baseline$conservative &&
    baseline_r <= estimates$in_vivo_baseline$optimistic

  in_vivo_treated_range <- treated_r >= estimates$in_vivo_treated$conservative &&
    treated_r <= estimates$in_vivo_treated$optimistic

  in_vitro_baseline_range <- baseline_r >= estimates$in_vitro_baseline$conservative &&
    baseline_r <= estimates$in_vitro_baseline$optimistic

  in_vitro_treated_range <- treated_r >= estimates$in_vitro_treated$conservative &&
    treated_r <= estimates$in_vitro_treated$optimistic

  if (in_vivo_baseline_range || in_vivo_treated_range ||
      in_vitro_baseline_range || in_vitro_treated_range) {
    return("Within experimental range")
  }

  return("Other")
}


# MODEL DEMONSTRATION FUNCTIONS -------------------------------------------


# Basic model demonstration data
create_model_demonstration_data <- function() {

  scenarios <- list(
    "Wild-type (no advantage)" = list(
      baseline_ratio = 1.0,
      treated_ratio = 1.0,
      color = "#2C3E50",
      description = "Normal parasites - no transmission advantage"
    ),
    "Moderate advantage" = list(
      baseline_ratio = 3.0,
      treated_ratio = 3.0,
      color = "#F39C12",
      description = "3x more infectious than wild-type"
    ),
    "Strong advantage" = list(
      baseline_ratio = 6.0,
      treated_ratio = 6.0,
      color = "#E74C3C",
      description = "6x more infectious than wild-type"
    )
  )

  EIR <- 50
  ft <- 0.6
  simulation_years <- 3

  results_list <- list()

  for (scenario_name in names(scenarios)) {
    params <- scenarios[[scenario_name]]

    model <- malaria_model(
      EIR = EIR,
      ft = ft,
      ton = 365,
      toff = 365 + (simulation_years * 365),
      day0_res = 0.01,
      treatment_failure_rate = 0.43,
      rT_r_cleared = 0.1,
      rT_r_failed = 0.1,
      resistance_baseline_ratio = params$baseline_ratio,
      resistance_cleared_ratio = params$treated_ratio,
      resistance_failed_ratio = params$treated_ratio
    )

    times <- seq(0, 365 + (simulation_years * 365), by = 30)
    output <- model$run(times)

    results_list[[scenario_name]] <- data.frame(
      time_years = output[, "t"] / 365,
      resistance_prevalence = output[, "prevalence_res"] * 100,
      total_prevalence = output[, "prevalence"] * 100,
      scenario = scenario_name,
      advantage_level = params$baseline_ratio,
      color = params$color
    )
  }

  combined_results <- do.call(rbind, results_list)
  return(combined_results)
}

# Transmission mechanism data
create_transmission_mechanism_data <- function() {

  model <- malaria_model(
    EIR = 50, ft = 0.6, ton = 365, toff = 365 + (2*365),
    day0_res = 0.05, treatment_failure_rate = 0.43,
    rT_r_cleared = 0.1, rT_r_failed = 0.1,
    resistance_baseline_ratio = 2.9,
    resistance_cleared_ratio = 0.89,
    resistance_failed_ratio = 0.89
  )

  times <- seq(0, 365 + (2*365), by = 14)
  output <- model$run(times)

  mechanism_data <- data.frame(
    time_years = output[, "t"] / 365,
    total_prevalence = output[, "prevalence"] * 100,
    resistant_prevalence = output[, "prevalence_res"] * 100,
    sensitive_prevalence = output[, "prevalence_sensitive"] * 100,
    treated_resistant = output[, "prevalence_resistant_treated"] * 100,
    EIR_sensitive = output[, "EIR_s"],
    EIR_resistant = output[, "EIR_r"],
    EIR_total = output[, "EIR_global"]
  )

  return(mechanism_data)
}

# Parameter sweep data
create_parameter_sweep_data <- function() {

  infectiousness_ratios <- c(1, 2, 4, 6, 8, 10)
  results_list <- list()

  for (ratio in infectiousness_ratios) {
    model <- malaria_model(
      EIR = 50, ft = 0.6, ton = 365, toff = 365 + (2*365),
      day0_res = 0.01, treatment_failure_rate = 0.43,
      rT_r_cleared = 0.1, rT_r_failed = 0.1,
      resistance_baseline_ratio = ratio,
      resistance_cleared_ratio = ratio,
      resistance_failed_ratio = ratio
    )

    times <- seq(0, 365 + (2*365), by = 30)
    output <- model$run(times)

    results_list[[paste0("ratio_", ratio)]] <- data.frame(
      time_years = output[, "t"] / 365,
      resistance_prevalence = output[, "prevalence_res"] * 100,
      infectiousness_ratio = ratio
    )
  }

  combined_results <- do.call(rbind, results_list)
  return(combined_results)
}

# Duration analysis data
create_duration_analysis_data <- function() {

  infectiousness_durations <- c(5, 10, 20, 50, 100)
  results_list <- list()

  for (duration in infectiousness_durations) {
    model <- malaria_model(
      EIR = 50, ft = 0.6, ton = 365,
      toff = 365 + duration,
      day0_res = 0.01, treatment_failure_rate = 0.43,
      rT_r_cleared = 0.1, rT_r_failed = 0.1,
      resistance_baseline_ratio = 2.9,
      resistance_cleared_ratio = 0.89,
      resistance_failed_ratio = 0.89
    )

    times <- seq(0, 365 + (3*365), by = 14)
    output <- model$run(times)

    results_list[[paste0("duration_", duration)]] <- data.frame(
      time_years = output[, "t"] / 365,
      prevalence_res = output[, "prevalence_res"] * 100,
      duration = duration,
      duration_label = paste(duration, "days")
    )
  }

  combined_results <- do.call(rbind, results_list)
  return(combined_results)
}


# SIMPLE TIME-SERIES COMPARISON DATA --------------------------------------


create_simple_timeseries_data <- function() {
  scenarios <- data.frame(
    study_type = c("In vitro", "In vitro", "In vivo", "In vivo", "In vivo"),
    confidence = c("Central", "Optimistic", "Conservative", "Central", "Optimistic"),
    baseline_ratio = c(1.0, 1.0, 2.08, 2.9, 3.98),
    treated_ratio = c(7.5, 8.0, 0.39, 0.89, 1.98),
    color = c("blue", "lightblue", "orange", "red", "darkred")
  )

  results_list <- list()

  for (i in 1:nrow(scenarios)) {
    scenario <- scenarios[i, ]

    model <- malaria_model(
      EIR = 50, ft = 0.6, ton = 365, toff = 365 + (4*365),
      day0_res = 0.01, treatment_failure_rate = 0.43,
      rT_r_cleared = 0.1, rT_r_failed = 0.1,
      resistance_baseline_ratio = scenario$baseline_ratio,
      resistance_cleared_ratio = scenario$treated_ratio,
      resistance_failed_ratio = scenario$treated_ratio
    )

    times <- seq(0, 365 + (4*365), by = 30)
    output <- model$run(times)

    results_list[[i]] <- data.frame(
      time_years = output[, "t"] / 365,
      resistance_prevalence = output[, "prevalence_res"] * 100,
      study_type = scenario$study_type,
      confidence = scenario$confidence,
      scenario = paste(scenario$study_type, scenario$confidence),
      baseline_ratio = scenario$baseline_ratio,
      treated_ratio = scenario$treated_ratio,
      color = scenario$color
    )
  }

  combined_data <- do.call(rbind, results_list)
  return(combined_data)
}


# COMPREHENSIVE HEATMAP DATA ----------------------------------------------


create_comprehensive_heatmap_data <- function() {

  EIR_values <- c(1, 5, 10, 20, 50, 100, 200, 400)
  ft_values <- seq(0.1, 0.9, 0.1)

  scenarios <- list(
    "Wild-type" = list(baseline = 1.0, treated = 1.0, color = "#2C3E50"),
    "In vivo" = list(baseline = 2.9, treated = 0.89, color = "#E74C3C"),
    "In vitro" = list(baseline = 1.0, treated = 7.5, color = "#3498DB"),
    "Combined" = list(baseline = 3.73, treated = 7.5, color = "#9B59B6")
  )

  param_grid <- expand.grid(
    EIR = EIR_values,
    ft = ft_values,
    scenario = names(scenarios)
  )

  results_list <- list()

  for (i in 1:nrow(param_grid)) {
    eir <- param_grid$EIR[i]
    ft <- param_grid$ft[i]
    scenario <- param_grid$scenario[i]
    params <- scenarios[[scenario]]

    if (i %% 50 == 0) cat(paste("Heatmap progress:", i, "/", nrow(param_grid), "\n"))

    model <- malaria_model(
      EIR = eir, ft = ft, ton = 365, toff = 365 + (5*365),
      day0_res = 0.01, treatment_failure_rate = 0.43,
      rT_r_cleared = 0.1, rT_r_failed = 0.1,
      resistance_baseline_ratio = params$baseline,
      resistance_cleared_ratio = params$treated,
      resistance_failed_ratio = params$treated
    )

    times <- seq(0, 365 + (5*365), by = 30)
    output <- model$run(times)

    final_resistance <- output[nrow(output), "prevalence_res"]

    results_list[[i]] <- data.frame(
      EIR = eir,
      ft = ft,
      scenario = scenario,
      final_resistance = final_resistance,
      baseline_ratio = params$baseline,
      treated_ratio = params$treated
    )
  }

  combined_results <- do.call(rbind, results_list)
  return(combined_results)
}


# COMPREHENSIVE PARAMETER ANALYSIS ----------------------------------------


comprehensive_parameter_analysis <- function() {
  estimates <- define_experimental_ranges()

  baseline_ratios <- c(
    0.8, 1.0, 1.5, 2.0, 3.0,
    estimates$in_vivo_baseline$conservative,
    4.5, 5.0,
    estimates$in_vivo_baseline$realistic,
    6.0, 7.0,
    estimates$in_vivo_baseline$optimistic,
    8.0, 9.0, 10.0
  )

  treated_ratios <- c(
    1.0, 1.5, 2.0,
    estimates$in_vivo_treated$conservative,
    3.0, 3.5,
    estimates$in_vivo_treated$realistic,
    4.5, 5.0, 6.0,
    estimates$in_vivo_treated$optimistic,
    estimates$in_vitro_treated$conservative,
    estimates$in_vitro_treated$realistic,
    estimates$in_vitro_treated$optimistic,
    9.0, 10.0
  )

  EIR_values <- c(10, 25, 50, 100)
  ft_values <- c(0.3, 0.5, 0.7)

  param_grid <- expand.grid(
    EIR = EIR_values,
    ft = ft_values,
    baseline_ratio = baseline_ratios,
    treated_ratio = treated_ratios
  )

  sample_size <- min(500, nrow(param_grid))
  sampled_indices <- sample(1:nrow(param_grid), sample_size)
  param_grid_sample <- param_grid[sampled_indices, ]

  results_list <- list()

  for (i in 1:nrow(param_grid_sample)) {
    combo <- param_grid_sample[i, ]

    if (i %% 50 == 0) cat(paste("Comprehensive progress:", i, "/", nrow(param_grid_sample), "\n"))

    tryCatch({
      model <- malaria_model(
        EIR = combo$EIR,
        ft = combo$ft,
        ton = 365,
        toff = 365 + (5*365),
        day0_res = 0.01,
        treatment_failure_rate = 0.43,
        rT_r_cleared = 0.1,
        rT_r_failed = 0.1,
        resistance_baseline_ratio = combo$baseline_ratio,
        resistance_cleared_ratio = combo$treated_ratio,
        resistance_failed_ratio = combo$treated_ratio
      )

      times <- seq(0, 365 + (5*365), by = 60)
      output <- model$run(times)

      final_resistance <- output[nrow(output), "prevalence_res"]

      p0_idx <- which.min(abs(output[, "t"] - 395))
      p1_idx <- which.min(abs(output[, "t"] - 730))

      p0 <- max(min(output[p0_idx, "prevalence_res"], 0.999), 0.001)
      p1 <- max(min(output[p1_idx, "prevalence_res"], 0.999), 0.001)

      sel_coeff <- log(p1*(1-p0)/(p0*(1-p1)))

      param_type <- classify_parameter_combination(combo$baseline_ratio, combo$treated_ratio, estimates)

      results_list[[i]] <- data.frame(
        EIR = combo$EIR,
        ft = combo$ft,
        baseline_ratio = combo$baseline_ratio,
        treated_ratio = combo$treated_ratio,
        final_resistance = final_resistance,
        selection_coefficient = sel_coeff,
        parameter_type = param_type,
        is_experimental_estimate = param_type != "Other"
      )

    }, error = function(e) {
      # Skip failed runs
    })
  }

  results_df <- do.call(rbind, results_list)
  results_df <- results_df[!is.na(results_df$final_resistance), ]

  return(results_df)
}


# ALL INDIVIDUAL PLOT FUNCTIONS -------------------------------------------


# Plot 1: Model demonstration
plot_model_demonstration <- function(demo_data) {
  p <- ggplot(demo_data, aes(x = time_years, y = resistance_prevalence,
                             color = scenario, linetype = scenario)) +
    geom_line(size = 2) +
    geom_vline(xintercept = 1, linetype = "dashed", alpha = 0.7, color = "black", size = 1) +
    scale_color_manual(values = setNames(unique(demo_data$color), unique(demo_data$scenario))) +
    scale_linetype_manual(values = c("solid", "longdash", "dotted")) +
    scale_x_continuous(breaks = seq(0, 4, 0.5)) +
    scale_y_continuous(limits = c(0, max(demo_data$resistance_prevalence) * 1.1),
                       breaks = seq(0, 100, 20)) +
    labs(
      title = "Model Demonstration: Impact of Transmission Advantages on Resistance Spread",
      subtitle = "EIR = 50 | Treatment coverage = 60 % | Resistance turned on at 1 year",
      x = "Time (years)",
      y = "Resistant infection prevalence (%)",
      color = "Parasite type",
      linetype = "Parasite type",
      caption = "Dashed vertical line shows when transmission advantages become active"
    ) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      plot.title = element_text(face = "bold", size = 16),
      plot.subtitle = element_text(size = 12),
      axis.title = element_text(face = "bold", size = 12),
      axis.text = element_text(size = 11),
      legend.title = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    ) +
    guides(
      color = guide_legend(override.aes = list(size = 3)),
      linetype = guide_legend(override.aes = list(size = 3))
    ) +
    annotate("text", x = 0.5, y = max(demo_data$resistance_prevalence) * 0.9,
             label = "Baseline period\n(no advantage)",
             hjust = 0.5, vjust = 1, size = 4, fontface = "italic") +
    annotate("text", x = 2.5, y = max(demo_data$resistance_prevalence) * 0.9,
             label = "Transmission advantage\nperiod",
             hjust = 0.5, vjust = 1, size = 4, fontface = "italic")

  return(p)
}

# Plot 2: Infection prevalence mechanism
plot_infection_prevalence <- function(mechanism_data) {
  prevalence_data <- mechanism_data %>%
    select(time_years, resistant_prevalence, sensitive_prevalence) %>%
    pivot_longer(cols = -time_years, names_to = "type", values_to = "prevalence") %>%
    mutate(
      type = case_when(
        type == "resistant_prevalence" ~ "Resistant infections",
        type == "sensitive_prevalence" ~ "Sensitive infections"
      )
    )

  p <- ggplot(prevalence_data, aes(x = time_years, y = prevalence, color = type, fill = type)) +
    geom_area(alpha = 0.6, position = "identity") +
    geom_line(size = 1.2) +
    geom_vline(xintercept = 1, linetype = "dashed", alpha = 0.7) +
    scale_color_manual(values = c("Resistant infections" = "#E74C3C", "Sensitive infections" = "#3498DB")) +
    scale_fill_manual(values = c("Resistant infections" = "#E74C3C", "Sensitive infections" = "#3498DB")) +
    labs(
      title = "Infection Prevalence Over Time",
      x = "Time (years)",
      y = "Prevalence (%)",
      color = "Infection type",
      fill = "Infection type"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")

  return(p)
}

# Plot 3: Transmission intensity by strain
plot_transmission_intensity <- function(mechanism_data) {
  EIR_data <- mechanism_data %>%
    select(time_years, EIR_sensitive, EIR_resistant) %>%
    pivot_longer(cols = -time_years, names_to = "strain", values_to = "EIR") %>%
    mutate(
      strain = case_when(
        strain == "EIR_sensitive" ~ "Sensitive strain transmission",
        strain == "EIR_resistant" ~ "Resistant strain transmission"
      )
    )

  p <- ggplot(EIR_data, aes(x = time_years, y = EIR, color = strain)) +
    geom_line(size = 1.2) +
    geom_vline(xintercept = 1, linetype = "dashed", alpha = 0.7) +
    scale_color_manual(values = c("Sensitive strain transmission" = "#3498DB", "Resistant strain transmission" = "#E74C3C")) +
    labs(
      title = "Transmission Intensity by Strain",
      x = "Time (years)",
      y = "Annual EIR",
      color = "Transmission type"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")

  return(p)
}

# Plot 4: Parameter effect sweep
plot_parameter_sweep <- function(sweep_data) {
  p <- ggplot(sweep_data, aes(x = time_years, y = resistance_prevalence,
                              color = factor(infectiousness_ratio))) +
    geom_line(size = 1.5) +
    geom_vline(xintercept = 1, linetype = "dashed", alpha = 0.7) +
    scale_color_viridis_d(name = "Infectiousness\nratio", option = "plasma") +
    labs(
      title = "Parameter Effect: Infectiousness Ratio Impact",
      x = "Time (years)",
      y = "Resistant infection prevalence (%)",
    ) +
    theme_minimal() +
    theme(legend.position = "right")

  return(p)
}

# Plot 5: Duration analysis
plot_duration_analysis <- function(duration_data) {
  p <- ggplot(duration_data, aes(x = time_years, y = prevalence_res, color = factor(duration))) +
    geom_line(size = 1.2) +
    geom_vline(xintercept = 1, linetype = "dashed", alpha = 0.7) +
    scale_color_viridis_d(name = "Infectiousness\nadvantage\nduration (days)", option = "plasma") +
    labs(
      title = "Impact of Infectiousness Advantage Duration",
      x = "Time (years)",
      y = "Resistant infection prevalence (%)",
    ) +
    theme_minimal() +
    theme(legend.position = "right")

  return(p)
}

# Plot 6: Simple timeseries comparison
plot_simple_timeseries <- function(timeseries_data) {
  p <- ggplot(timeseries_data, aes(x = time_years, y = resistance_prevalence,
                                   color = scenario, linetype = study_type)) +
    geom_line(size = 1.5) +
    geom_vline(xintercept = 1, linetype = "dashed", alpha = 0.7, color = "black") +
    scale_color_manual(values = setNames(unique(timeseries_data$color), unique(timeseries_data$scenario))) +
    scale_linetype_manual(values = c("In vitro" = "solid", "In vivo" = "longdash")) +
    labs(
      title = "Resistance Prevalence: In Vivo vs In Vitro Estimates",
      x = "Time (years)",
      y = "Resistant infection prevalence (%)",
      color = "Experimental scenario",
      linetype = "Study type",
      caption = "Dashed vertical line: when resistance advantages turn on"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      legend.position = "bottom"
    )

  return(p)
}

# Plot 7: Comprehensive heatmap
plot_comprehensive_heatmap <- function(heatmap_data) {
  p <- ggplot(heatmap_data, aes(x = factor(ft), y = factor(EIR), fill = final_resistance)) +
    geom_tile(color = "white", size = 0.2) +
    geom_text(aes(label = round(final_resistance, 2)), size = 2.5, color = "white") +
    facet_wrap(~scenario, nrow = 2) +
    scale_fill_viridis_c(name = "Final\nresistance\nprevalence", option = "plasma") +
    scale_x_discrete(breaks = seq(0.1, 0.9, 0.2)) +
    labs(
      title = "Gametocytogenesis Impact on Resistance Spread Across Epidemiological Settings",
      subtitle = "Final resistance prevalence after 5 years",
      x = "Treatment coverage (ft)",
      y = "Annual EIR",
      caption = "Numbers show final resistance prevalence"
    ) +
    theme_minimal() +
    theme(
      strip.text = element_text(face = "bold"),
      plot.title = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )

  return(p)
}

# Plot 8: Detailed parameter heatmap
plot_detailed_heatmap <- function(results_df) {
  p <- ggplot(results_df, aes(x = baseline_ratio, y = treated_ratio,
                              fill = final_resistance)) +
    geom_point(aes(size = ifelse(is_experimental_estimate, 3, 1),
                   color = parameter_type), alpha = 0.8) +
    facet_grid(EIR ~ ft, labeller = labeller(
      EIR = function(x) paste("EIR =", x),
      ft = function(x) paste("Treatment =", round(as.numeric(x)*100), "%")
    )) +
    scale_fill_viridis_c(name = "Final\nresistance\nprevalence", limits = c(0, 1)) +
    scale_size_identity() +
    scale_color_manual(
      name = "Parameter type",
      values = c(
        "In vivo (central)" = "red",
        "In vivo (conservative)" = "orange",
        "In vivo (optimistic)" = "darkred",
        "In vitro (central)" = "blue",
        "In vitro (conservative)" = "lightblue",
        "In vitro (optimistic)" = "darkblue",
        "Within experimental range" = "purple",
        "Other" = "gray"
      )
    ) +
    labs(
      title = "Resistance Spread Across Parameter Space",
      subtitle = "Colored points show experimental estimates and ranges",
      x = "Baseline infectiousness ratio",
      y = "Post-treatment infectiousness ratio",
      caption = "Larger points = exact experimental estimates"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "bottom"
    ) +
    guides(color = guide_legend(override.aes = list(size = 3)))

  return(p)
}

# Plot 9: Selection coefficient landscape
plot_selection_landscape <- function(results_df) {
  p <- ggplot(results_df, aes(x = baseline_ratio, y = treated_ratio,
                              fill = selection_coefficient)) +
    geom_point(aes(size = ifelse(is_experimental_estimate, 3, 1),
                   color = parameter_type), alpha = 0.8) +
    facet_grid(EIR ~ ft, labeller = labeller(
      EIR = function(x) paste("EIR =", x),
      ft = function(x) paste("Treatment =", round(as.numeric(x)*100), "%")
    )) +
    scale_fill_gradient2(name = "Selection\ncoefficient\n(per year)",
                         low = "blue", mid = "white", high = "red", midpoint = 0) +
    scale_size_identity() +
    scale_color_manual(
      name = "Parameter type",
      values = c(
        "In vivo (central)" = "red",
        "In vivo (conservative)" = "orange",
        "In vivo (optimistic)" = "darkred",
        "In vitro (central)" = "blue",
        "In vitro (conservative)" = "lightblue",
        "In vitro (optimistic)" = "darkblue",
        "Within experimental range" = "purple",
        "Other" = "gray"
      )
    ) +
    labs(
      title = "Selection Coefficient Landscape",
      subtitle = "Quantifying resistance advantages across parameter space",
      x = "Baseline infectiousness ratio",
      y = "Post-treatment infectiousness ratio"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "bottom"
    ) +
    guides(color = guide_legend(override.aes = list(size = 3)))

  return(p)
}

# Plot 10: Experimental summary
plot_experimental_summary <- function(results_df) {
  exp_data <- results_df[results_df$is_experimental_estimate, ]

  if (nrow(exp_data) > 0) {
    p <- ggplot(exp_data, aes(x = parameter_type, y = final_resistance,
                              fill = parameter_type)) +
      geom_boxplot(alpha = 0.7) +
      geom_jitter(width = 0.2, alpha = 0.6) +
      scale_fill_manual(
        values = c(
          "In vivo (central)" = "red",
          "In vivo (conservative)" = "orange",
          "In vivo (optimistic)" = "darkred",
          "In vitro (central)" = "blue",
          "In vitro (conservative)" = "lightblue",
          "In vitro (optimistic)" = "darkblue",
          "Within experimental range" = "purple"
        )
      ) +
      labs(
        title = "Final Resistance for Experimental Estimates",
        subtitle = "Comparing in vivo vs in vitro scenarios across confidence intervals",
        x = "Parameter scenario",
        y = "Final resistance prevalence",
        fill = "Scenario"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none"
      )

    return(p)
  }
  return(NULL)
}

# Plot 11: Selection coefficient comparison
plot_selection_comparison <- function() {
  cat("Creating selection coefficient comparison...\n")

  estimates <- data.frame(
    scenario = c("In vitro\n(central)", "In vitro\n(optimistic)",
                 "In vivo\n(conservative)", "In vivo\n(central)", "In vivo\n(optimistic)"),
    baseline_ratio = c(1.0, 1.0, 2.08, 2.9, 3.98),
    treated_ratio = c(7.5, 8.0, 0.39, 0.89, 1.98),
    study_type = c("In vitro", "In vitro", "In vivo", "In vivo", "In vivo"),
    confidence = c("Central", "Optimistic", "Conservative", "Central", "Optimistic")
  )

  selection_results <- list()

  for (i in 1:nrow(estimates)) {
    est <- estimates[i, ]

    model <- malaria_model(
      EIR = 50, ft = 0.6, ton = 365, toff = 365 + (2*365),
      day0_res = 0.01, treatment_failure_rate = 0.43,
      rT_r_cleared = 0.1, rT_r_failed = 0.1,
      resistance_baseline_ratio = est$baseline_ratio,
      resistance_cleared_ratio = est$treated_ratio,
      resistance_failed_ratio = est$treated_ratio
    )

    times <- seq(0, 365 + (2*365), by = 14)
    output <- model$run(times)

    p0_idx <- which.min(abs(output[, "t"] - 395))
    p1_idx <- which.min(abs(output[, "t"] - 730))

    p0 <- max(min(output[p0_idx, "prevalence_res"], 0.999), 0.001)
    p1 <- max(min(output[p1_idx, "prevalence_res"], 0.999), 0.001)

    sel_coeff <- log(p1*(1-p0)/(p0*(1-p1)))

    selection_results[[i]] <- data.frame(
      scenario = est$scenario,
      study_type = est$study_type,
      confidence = est$confidence,
      baseline_ratio = est$baseline_ratio,
      treated_ratio = est$treated_ratio,
      selection_coefficient = sel_coeff
    )
  }

  selection_data <- do.call(rbind, selection_results)

  colors <- c("In vitro" = "#3498DB", "In vivo" = "#E74C3C")

  p <- ggplot(selection_data, aes(x = scenario, y = selection_coefficient,
                                  fill = study_type, alpha = confidence)) +
    geom_col(color = "black", size = 0.5) +
    scale_fill_manual(values = colors) +
    scale_alpha_manual(values = c("Conservative" = 0.6, "Central" = 1.0, "Optimistic" = 0.8)) +
    labs(
      title = "Selection Coefficients: Quantifying Resistance Advantages",
      x = "Experimental scenario",
      y = "Selection coefficient (per year)",
      fill = "Study type",
      alpha = "Confidence level",
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )

  return(p)
}

# Plot 12: Key findings summary
plot_key_findings <- function(timeseries_data, selection_data) {
  # Extract key numbers
  in_vitro_final <- timeseries_data %>%
    filter(study_type == "In vitro", confidence == "Central", time_years == max(time_years)) %>%
    pull(resistance_prevalence)

  in_vivo_final <- timeseries_data %>%
    filter(study_type == "In vivo", confidence == "Central", time_years == max(time_years)) %>%
    pull(resistance_prevalence)

  in_vitro_sel <- selection_data %>%
    filter(study_type == "In vitro", confidence == "Central") %>%
    pull(selection_coefficient)

  in_vivo_sel <- selection_data %>%
    filter(study_type == "In vivo", confidence == "Central") %>%
    pull(selection_coefficient)

  summary_data <- data.frame(
    metric = rep(c("Final resistance (%)", "Selection coefficient"), each = 2),
    study_type = rep(c("In vitro", "In vivo"), 2),
    value = c(in_vitro_final, in_vivo_final, in_vitro_sel, in_vivo_sel),
    interpretation = c("Minimal spread", "Complete fixation", "Weak selection", "Strong selection")
  )

  p <- ggplot(summary_data, aes(x = study_type, y = value, fill = study_type)) +
    geom_col(alpha = 0.8, color = "black") +
    geom_text(aes(label = paste0(round(value, 1), "\n(", interpretation, ")")),
              vjust = 0.5, size = 3.5, fontface = "bold") +
    facet_wrap(~metric, scales = "free_y") +
    scale_fill_manual(values = c("In vitro" = "#3498DB", "In vivo" = "#E74C3C")) +
    labs(
      title = "Key Finding: In Vivo vs In Vitro Estimates",
      x = "Study type",
      y = "Value",
      fill = "Study type"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "none",
      strip.text = element_text(face = "bold")
    )

  return(p)
}


# MAIN ANALYSIS FUNCTION --------------------------------------------------

run_complete_gametocyte_analysis <- function() {
  demo_data <- create_model_demonstration_data()
  mechanism_data <- create_transmission_mechanism_data()
  sweep_data <- create_parameter_sweep_data()
  duration_data <- create_duration_analysis_data()
  timeseries_data <- create_simple_timeseries_data()
  heatmap_data <- create_comprehensive_heatmap_data()
  comprehensive_data <- comprehensive_parameter_analysis()
  cat("Creating and saving all individual plots...\n")
  # Plot 1: Model demonstration
  p1 <- plot_model_demonstration(demo_data)
  print(p1)
  ggsave("01_model_demonstration.png", p1, width = 12, height = 8, dpi = 300)
  # Plot 2: Infection prevalence mechanism
  p2 <- plot_infection_prevalence(mechanism_data)
  print(p2)
  ggsave("02_infection_prevalence.png", p2, width = 10, height = 6, dpi = 300)
  # Plot 3: Transmission intensity
  p3 <- plot_transmission_intensity(mechanism_data)
  print(p3)
  ggsave("03_transmission_intensity.png", p3, width = 10, height = 6, dpi = 300)
  # Plot 4: Parameter sweep
  p4 <- plot_parameter_sweep(sweep_data)
  print(p4)
  ggsave("04_parameter_sweep.png", p4, width = 10, height = 8, dpi = 300)
  # Plot 5: Duration analysis
  p5 <- plot_duration_analysis(duration_data)
  print(p5)
  ggsave("05_duration_analysis.png", p5, width = 10, height = 8, dpi = 300)
  # Plot 6: Simple timeseries comparison
  p6 <- plot_simple_timeseries(timeseries_data)
  print(p6)
  ggsave("06_simple_timeseries.png", p6, width = 12, height = 8, dpi = 300)
  # Plot 7: Comprehensive heatmap
  p7 <- plot_comprehensive_heatmap(heatmap_data)
  print(p7)
  ggsave("07_comprehensive_heatmap.png", p7, width = 12, height = 10, dpi = 300)
  # Plot 8: Detailed parameter heatmap
  p8 <- plot_detailed_heatmap(comprehensive_data)
  print(p8)
  ggsave("08_detailed_heatmap.png", p8, width = 14, height = 10, dpi = 300)
  # Plot 9: Selection landscape
  p9 <- plot_selection_landscape(comprehensive_data)
  print(p9)
  ggsave("09_selection_landscape.png", p9, width = 14, height = 10, dpi = 300)
  # Plot 10: Experimental summary
  p10 <- plot_experimental_summary(comprehensive_data)
  if (!is.null(p10)) {
    print(p10)
    ggsave("10_experimental_summary.png", p10, width = 10, height = 8, dpi = 300)
  }
  # Plot 11: Selection comparison
  p11 <- plot_selection_comparison()
  print(p11)
  ggsave("11_selection_comparison.png", p11, width = 10, height = 8, dpi = 300)
  # Plot 12: Key findings (need to create selection data first)
  selection_estimates <- data.frame(
    study_type = c("In vitro", "In vivo"),
    confidence = c("Central", "Central"),
    selection_coefficient = c(0.3, 10.4)  # Approximate values from your images
  )
  p12 <- plot_key_findings(timeseries_data, selection_estimates)
  print(p12)
  ggsave("12_key_findings.png", p12, width = 10, height = 8, dpi = 300)
}


# RUNNING THE COMPLETE ANALYSIS -----------------------------------------------

complete_results <- run_complete_gametocyte_analysis()

