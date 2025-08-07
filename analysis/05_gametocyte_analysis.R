
devtools::load_all()

library(mmsu)
library(ggplot2)
library(dplyr)
library(tidyr)
library(viridis)
library(patchwork)
library(scales)


# Parameter ranges --------------------------------------------------------

define_experimental_ranges <- function() {
  estimates <- list(
    in_vivo_baseline = list(
      conservative = 0.903,
      central = 1.265,
      optimistic = 1.759
    ),
    in_vivo_treated = list(
      conservative = 0.903,
      central = 1.265,
      optimistic = 1.759
    ),
    in_vitro_baseline = list(
      conservative = 0.8,
      central = 1.0,
      optimistic = 1.0
    ),
    in_vitro_treated = list(
      conservative = 7.0,
      central = 7.5,
      optimistic = 8.0
    )
  )
  return(estimates)
}



# Data generation functions -----------------------------------------------

create_model_demonstration_data <- function(EIR = 40, ft = 0.34, simulation_years = 10,
                                            treatment_failure_rate = 0.43, day0_res = 0.01) {
  scenarios <- list(
    "No advantage" = list(
      baseline_ratio = 1.0,
      treated_ratio = 1.0,
      color = "#2C3E50",
      description = "Normal parasites - no transmission advantage"
    ),
    "Moderate advantage (2x)" = list(
      baseline_ratio = 2.0,
      treated_ratio = 2.0,
      color = "#F39C12",
      description = "2x more infectious than normal"
    ),
    "Strong advantage (4x)" = list(
      baseline_ratio = 4.0,
      treated_ratio = 4.0,
      color = "#E74C3C",
      description = "4x more infectious than normal"
    )
  )

  results_list <- list()
  for (scenario_name in names(scenarios)) {
    params <- scenarios[[scenario_name]]

    model <- malaria_model(
      EIR = EIR,
      ft = ft,
      ton = 365,
      toff = 365 + (simulation_years * 365),
      day0_res = day0_res,
      treatment_failure_rate = treatment_failure_rate,
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

create_simple_timeseries_data <- function() {
  # Updated with consistent parameter values
  scenarios <- data.frame(
    study_type = c("In vitro", "In vitro", "In vitro", "In vivo", "In vivo", "In vivo"),
    confidence = c("Conservative", "Central", "Optimistic", "Conservative", "Central", "Optimistic"),
    baseline_ratio = c(0.8, 1.0, 1.0, 0.903, 1.265, 1.759),
    treated_ratio = c(7.0, 7.5, 8.0, 0.903, 1.265, 1.759),
    color = c("darkblue", "blue", "lightblue", "orange", "red", "darkred")
  )

  results_list <- list()
  for (i in 1:nrow(scenarios)) {
    scenario <- scenarios[i, ]
    model <- malaria_model(
      EIR = 40, ft = 0.34,
      ton = 365, toff = 365 + (4*365),
      day0_res = 0.01,
      treatment_failure_rate = 0.43,
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


# Classification function for heatmaps ------------------------------------

classify_parameter_combination <- function(baseline_r, treated_r, estimates) {
  tol <- 0.05

  # In vivo scenarios
  if (abs(baseline_r - estimates$in_vivo_baseline$central) < tol &&
      abs(treated_r - estimates$in_vivo_treated$central) < tol) {
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
  if (abs(baseline_r - estimates$in_vitro_baseline$central) < tol &&
      abs(treated_r - estimates$in_vitro_treated$central) < tol) {
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

  return("Other")
}


# Parameter analysis ------------------------------------------------------

comprehensive_parameter_analysis <- function() {
  estimates <- define_experimental_ranges()

  # Use the consistent parameter values
  baseline_ratios <- c(
    0.8, 1.0,
    estimates$in_vivo_baseline$conservative,   # 0.903
    estimates$in_vivo_baseline$central,        # 1.265
    estimates$in_vivo_baseline$optimistic,     # 1.759
    2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0
  )

  treated_ratios <- c(
    estimates$in_vivo_treated$conservative,    # 0.903
    1.0,
    estimates$in_vivo_treated$central,         # 1.265
    estimates$in_vivo_treated$optimistic,      # 1.759
    2.0, 3.0, 4.0, 5.0, 6.0,
    estimates$in_vitro_treated$conservative,   # 7.0
    estimates$in_vitro_treated$central,        # 7.5
    estimates$in_vitro_treated$optimistic,     # 8.0
    9.0, 10.0
  )

  # Updated broader ranges
  EIR_values <- c(1, 10, 25, 50, 100, 200, 400, 600)  # 0-600 range
  ft_values <- c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)  # 20-90% range

  param_grid <- expand.grid(
    EIR = EIR_values,
    ft = ft_values,
    baseline_ratio = baseline_ratios,
    treated_ratio = treated_ratios
  )

  # ENSURE experimental estimates are included for ALL EIR/ft combinations
  experimental_combos <- data.frame(
    baseline_ratio = c(1.265, 1.0, 0.903, 1.759, 0.8, 1.0),  # All experimental baselines
    treated_ratio = c(1.265, 7.5, 0.903, 1.759, 7.0, 8.0)    # All experimental treated
  )

  # Create explicit experimental combinations for each EIR/ft
  exp_grid <- expand.grid(
    EIR = EIR_values,
    ft = ft_values,
    baseline_ratio = experimental_combos$baseline_ratio,
    treated_ratio = experimental_combos$treated_ratio
  )

  # Filter to only valid experimental combinations
  exp_grid <- exp_grid %>%
    filter(
      (baseline_ratio == 1.265 & treated_ratio == 1.265) |  # In vivo central
        (baseline_ratio == 0.903 & treated_ratio == 0.903) |  # In vivo conservative
        (baseline_ratio == 1.759 & treated_ratio == 1.759) |  # In vivo optimistic
        (baseline_ratio == 1.0 & treated_ratio == 7.5) |      # In vitro central
        (baseline_ratio == 0.8 & treated_ratio == 7.0) |      # In vitro conservative
        (baseline_ratio == 1.0 & treated_ratio == 8.0)        # In vitro optimistic
    )

  # Combine with regular grid and remove duplicates
  param_grid <- rbind(param_grid, exp_grid)
  param_grid <- unique(param_grid)

  # Don't sample - use more points to ensure experimental estimates appear
  sample_size <- min(1000, nrow(param_grid))  # Increased sample size further
  if (nrow(param_grid) > sample_size) {
    # Always include ALL experimental estimates
    exp_indices <- which(
      (param_grid$baseline_ratio == 1.265 & param_grid$treated_ratio == 1.265) |  # In vivo central
        (param_grid$baseline_ratio == 0.903 & param_grid$treated_ratio == 0.903) |  # In vivo conservative
        (param_grid$baseline_ratio == 1.759 & param_grid$treated_ratio == 1.759) |  # In vivo optimistic
        (param_grid$baseline_ratio == 1.0 & param_grid$treated_ratio == 7.5) |      # In vitro central
        (param_grid$baseline_ratio == 0.8 & param_grid$treated_ratio == 7.0) |      # In vitro conservative
        (param_grid$baseline_ratio == 1.0 & param_grid$treated_ratio == 8.0)        # In vitro optimistic
    )

    other_indices <- setdiff(1:nrow(param_grid), exp_indices)
    sampled_other <- sample(other_indices, sample_size - length(exp_indices))

    param_grid_sample <- param_grid[c(exp_indices, sampled_other), ]
  } else {
    param_grid_sample <- param_grid
  }

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

      # Calculate selection coefficient
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

  # Verify experimental estimates are present
  exp_summary <- results_df %>%
    filter(is_experimental_estimate == TRUE) %>%
    group_by(parameter_type, EIR, ft) %>%
    summarise(count = n(), .groups = 'drop')
  print(exp_summary)

  return(results_df)
}

# Additional  -------------------------------------------------------------

create_transmission_mechanism_data <- function() {
  model <- malaria_model(
    EIR = 40, ft = 0.34,
    ton = 365, toff = 365 + (2*365),
    day0_res = 0.05,
    treatment_failure_rate = 0.43,
    rT_r_cleared = 0.1, rT_r_failed = 0.1,
    resistance_baseline_ratio = 1.265,  # Central value
    resistance_cleared_ratio = 1.265,
    resistance_failed_ratio = 1.265
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

create_parameter_sweep_data <- function() {
  infectiousness_ratios <- c(1, 2, 4, 6, 8, 10)
  results_list <- list()

  for (ratio in infectiousness_ratios) {
    model <- malaria_model(
      EIR = 40, ft = 0.34,
      ton = 365, toff = 365 + (2*365),
      day0_res = 0.01,
      treatment_failure_rate = 0.43,
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

create_duration_analysis_data <- function() {
  infectiousness_durations <- c(5, 10, 20, 50, 100)
  results_list <- list()

  for (duration in infectiousness_durations) {
    model <- malaria_model(
      EIR = 40, ft = 0.34,
      ton = 365,
      toff = 365 + duration,
      day0_res = 0.01,
      treatment_failure_rate = 0.43,
      rT_r_cleared = 0.1, rT_r_failed = 0.1,
      resistance_baseline_ratio = 1.265,  # Central value
      resistance_cleared_ratio = 1.265,
      resistance_failed_ratio = 1.265
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



# Heatmap data ------------------------------------------------------------

create_comprehensive_heatmap_data <- function() {
  EIR_values <- c(1, 5, 10, 20, 50, 100, 200, 400)
  ft_values <- seq(0.1, 0.9, 0.1)

  scenarios <- list(
    "Wild-type (no advantage)" = list(
      baseline = 1.0,
      treated = 1.0,
      color = "#2C3E50",
      treatment_failure_rate = 0.0  # No treatment failure for wild-type
    ),
    "In vivo" = list(
      baseline = 1.265,  # Central value
      treated = 1.265,   # Central value
      color = "#E74C3C",
      treatment_failure_rate = 0.43
    ),
    "In vitro" = list(
      baseline = 1.0,    # No baseline advantage
      treated = 7.5,     # Large post-treatment advantage
      color = "#3498DB",
      treatment_failure_rate = 0.43
    ),
    "Combined" = list(
      baseline = 1.265,  # In vivo baseline advantage
      treated = 7.5,     # In vitro post-treatment advantage
      color = "#9B59B6",
      treatment_failure_rate = 0.43
    )
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

    if (i %% 50 == 0) cat(paste("Grid heatmap progress:", i, "/", nrow(param_grid), "\n"))

    model <- malaria_model(
      EIR = eir,
      ft = ft,
      ton = 365,
      toff = 365 + (5*365),
      day0_res = 0.01,
      treatment_failure_rate = params$treatment_failure_rate,
      rT_r_cleared = 0.1,
      rT_r_failed = 0.1,
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
      treated_ratio = params$treated,
      treatment_failure_rate = params$treatment_failure_rate
    )
  }

  combined_results <- do.call(rbind, results_list)
  return(combined_results)
}


# Additional --------------------------------------------------------------

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
    geom_vline(xintercept = 1, linetype = "dashed", alpha = 0.7, color = "black") +
    scale_color_manual(values = c("Resistant infections" = "#E74C3C", "Sensitive infections" = "#3498DB")) +
    scale_fill_manual(values = c("Resistant infections" = "#E74C3C", "Sensitive infections" = "#3498DB")) +
    labs(
      title = "Infection Prevalence Over Time: Sensitive vs Resistant Strains",
      subtitle = "EIR = 40 | Treatment coverage = 34% | Resistance advantage active after year 1",
      x = "Time (years)",
      y = "Prevalence (%)",
      color = "Infection type",
      fill = "Infection type",
      caption = "Dashed line shows when transmission advantages become active"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      legend.position = "bottom"
    )

  return(p)
}

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
    geom_line(size = 1.5) +
    geom_vline(xintercept = 1, linetype = "dashed", alpha = 0.7, color = "black") +
    scale_color_manual(values = c("Sensitive strain transmission" = "#3498DB",
                                  "Resistant strain transmission" = "#E74C3C")) +
    labs(
      title = "Transmission Intensity by Strain Over Time",
      subtitle = "Annual EIR contribution from each parasite strain",
      x = "Time (years)",
      y = "Annual EIR",
      color = "Transmission type",
      caption = "Dashed line shows when resistance advantages become active"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      legend.position = "bottom"
    )

  return(p)
}

plot_parameter_sweep <- function(sweep_data) {
  p <- ggplot(sweep_data, aes(x = time_years, y = resistance_prevalence,
                              color = factor(infectiousness_ratio))) +
    geom_line(size = 1.5) +
    geom_vline(xintercept = 1, linetype = "dashed", alpha = 0.7, color = "black") +
    scale_color_viridis_d(name = "Infectiousness\nratio", option = "plasma") +
    labs(
      title = "Impact of Infectiousness Ratio on Resistance Spread",
      subtitle = "Parameter sweep showing effect of different transmission advantages",
      x = "Time (years)",
      y = "Resistant infection prevalence (%)",
      caption = "Dashed line shows when resistance advantages become active"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      legend.position = "right"
    )

  return(p)
}

plot_duration_analysis <- function(duration_data) {
  p <- ggplot(duration_data, aes(x = time_years, y = prevalence_res, color = factor(duration))) +
    geom_line(size = 1.2) +
    geom_vline(xintercept = 1, linetype = "dashed", alpha = 0.7, color = "black") +
    scale_color_viridis_d(name = "Advantage\nduration\n(days)", option = "plasma") +
    labs(
      title = "Impact of Infectiousness Advantage Duration",
      subtitle = "How long the transmission advantage lasts affects resistance spread",
      x = "Time (years)",
      y = "Resistant infection prevalence (%)",
      caption = "Dashed line shows when resistance advantages start (duration varies by scenario)"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      legend.position = "right"
    )

  return(p)
}


# Heatmap plot functions --------------------------------------------------

plot_comprehensive_heatmap_grid <- function(heatmap_data) {
  p <- ggplot(heatmap_data, aes(x = factor(ft), y = factor(EIR))) +
    geom_tile(aes(fill = final_resistance), color = "white", size = 0.5) +
    geom_text(aes(label = sprintf("%.2f", final_resistance)),
              color = "white", fontface = "bold", size = 3) +

    facet_wrap(~scenario, ncol = 2, scales = "fixed") +

    scale_fill_viridis_c(
      name = "Final\nresistance\nprevalence",
      option = "plasma",
      labels = scales::percent,
      limits = c(0, 1)
    ) +

    scale_x_discrete(
      name = "Treatment coverage (%)",
      labels = function(x) paste0(as.numeric(x) * 100, "%")
    ) +

    scale_y_discrete(
      name = "Annual EIR",
      labels = function(x) as.character(x)
    ) +

    labs(
      title = "Final Resistance Prevalence: Parameter Grid",
      subtitle = "5-year simulation across transmission intensities and treatment coverages",
      caption = "Values show final resistance prevalence after 5 years"
    ) +

    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 16),
      plot.subtitle = element_text(size = 12),
      strip.text = element_text(face = "bold", size = 12),
      axis.title = element_text(face = "bold", size = 11),
      axis.text = element_text(size = 10),
      legend.position = "right",
      panel.spacing = unit(1, "lines"),
      strip.background = element_rect(fill = "grey90", color = "white")
    ) +

    guides(
      fill = guide_colorbar(
        barwidth = 1.5,
        barheight = 10,
        title.position = "top",
        title.hjust = 0.5
      )
    )

  return(p)
}

plot_detailed_heatmap_clean <- function(results_df) {
  # Filter to only central estimates and "Other" to reduce overcrowding
  central_data <- results_df %>%
    filter(parameter_type %in% c("In vivo (central)", "In vitro (central)", "Other")) %>%
    mutate(
      parameter_category = case_when(
        parameter_type == "In vivo (central)" ~ "In vivo central",
        parameter_type == "In vitro (central)" ~ "In vitro central",
        TRUE ~ "Other combinations"
      ),
      point_size = ifelse(is_experimental_estimate, 3.5, 2)
    )

  p <- ggplot(central_data, aes(x = baseline_ratio, y = treated_ratio)) +
    geom_point(aes(color = final_resistance,
                   shape = parameter_category,
                   size = point_size),
               alpha = 0.8, stroke = 0.8) +

    scale_color_viridis_c(name = "Final resistance\nprevalence",
                          limits = c(0, 1),
                          option = "plasma",
                          labels = scales::percent) +

    scale_shape_manual(
      name = "Parameter type",
      values = c("In vivo central" = 18,      # Diamond
                 "In vitro central" = 15,     # Square
                 "Other combinations" = 16)   # Circle
    ) +

    scale_size_identity() +

    facet_grid(EIR ~ ft,
               labeller = labeller(
                 EIR = function(x) paste("EIR =", x),
                 ft = function(x) paste("Treatment =", scales::percent(as.numeric(x)))
               )) +

    labs(
      title = "Final Resistance Prevalence Across Parameter Space",
      subtitle = "Central experimental estimates highlighted (larger points)",
      x = "Baseline infectiousness ratio (κB)",
      y = "Post-treatment infectiousness ratio",
      caption = "Diamond = in vivo central; Square = in vitro central; Circle = other parameter combinations"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 12),
      legend.position = "right",
      strip.text = element_text(size = 10, face = "bold"),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "grey80", fill = NA, size = 0.5)
    ) +

    guides(
      color = guide_colorbar(barwidth = 1.2, barheight = 8),
      shape = guide_legend(override.aes = list(size = 4),
                           title.position = "top")
    )

  return(p)
}

plot_selection_landscape_clean <- function(results_df) {
  central_data <- results_df %>%
    filter(parameter_type %in% c("In vivo (central)", "In vitro (central)", "Other")) %>%
    mutate(
      selection_category = case_when(
        selection_coefficient < 0 ~ "Negative selection",
        selection_coefficient < 1 ~ "Weak (0-1)",
        selection_coefficient < 3 ~ "Moderate (1-3)",
        selection_coefficient < 7 ~ "Strong (3-7)",
        selection_coefficient < 15 ~ "Very strong (7-15)",
        TRUE ~ "Extreme (>15)"
      ),
      selection_category = factor(selection_category, levels = c(
        "Negative selection", "Weak (0-1)", "Moderate (1-3)",
        "Strong (3-7)", "Very strong (7-15)", "Extreme (>15)"
      )),
      parameter_category = case_when(
        parameter_type == "In vivo (central)" ~ "In vivo central",
        parameter_type == "In vitro (central)" ~ "In vitro central",
        TRUE ~ "Other combinations"
      ),
      point_size = ifelse(is_experimental_estimate, 3.5, 2.5)
    )

  p <- ggplot(central_data, aes(x = baseline_ratio, y = treated_ratio)) +
    geom_point(aes(fill = selection_category,
                   shape = parameter_category,
                   size = point_size),
               alpha = 0.8, stroke = 0.8, color = "white") +

    facet_grid(EIR ~ ft,
               labeller = labeller(
                 EIR = function(x) paste("EIR =", x),
                 ft = function(x) paste("Treatment =", scales::percent(as.numeric(x)))
               )) +

    scale_fill_manual(
      name = "Selection coefficient\n(per year)",
      values = c(
        "Negative selection" = "#3498DB",
        "Weak (0-1)" = "#95A5A6",
        "Moderate (1-3)" = "#F39C12",
        "Strong (3-7)" = "#E67E22",
        "Very strong (7-15)" = "#E74C3C",
        "Extreme (>15)" = "#8E44AD"
      )
    ) +

    scale_shape_manual(
      name = "Parameter type",
      values = c("In vivo central" = 23,      # Diamond with fill
                 "In vitro central" = 22,     # Square with fill
                 "Other combinations" = 21)   # Circle with fill
    ) +

    scale_size_identity() +

    labs(
      title = "Selection Coefficient Landscape",
      subtitle = "Discrete categories with experimental estimates highlighted",
      x = "Baseline infectiousness ratio (κB)",
      y = "Post-treatment infectiousness ratio"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 12),
      legend.position = "right",
      strip.text = element_text(size = 10, face = "bold"),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "grey80", fill = NA, size = 0.5)
    ) +

    guides(
      fill = guide_legend(override.aes = list(size = 4, shape = 21, color = "white"),
                          title.position = "top"),
      shape = guide_legend(override.aes = list(size = 4, fill = "grey", color = "white"),
                           title.position = "top")
    )

  return(p)
}

plot_experimental_summary_clean <- function(results_df) {
  # Get all experimental estimates (not just central ones)
  exp_data <- results_df %>%
    filter(is_experimental_estimate == TRUE) %>%
    mutate(
      # Create cleaner labels
      study_label = case_when(
        parameter_type == "In vivo (conservative)" ~ "In vivo conservative",
        parameter_type == "In vivo (central)" ~ "In vivo central",
        parameter_type == "In vivo (optimistic)" ~ "In vivo optimistic",
        parameter_type == "In vitro (conservative)" ~ "In vitro conservative",
        parameter_type == "In vitro (central)" ~ "In vitro central",
        parameter_type == "In vitro (optimistic)" ~ "In vitro optimistic",
        TRUE ~ parameter_type
      ),
      # Group by study type for coloring
      study_type = case_when(
        grepl("In vivo", parameter_type) ~ "In vivo",
        grepl("In vitro", parameter_type) ~ "In vitro",
        TRUE ~ "Other"
      )
    ) %>%
    # Get average values across EIR/ft combinations for each parameter set
    group_by(study_label, study_type, baseline_ratio, treated_ratio) %>%
    summarise(
      final_resistance = mean(final_resistance, na.rm = TRUE),
      selection_coefficient = mean(selection_coefficient, na.rm = TRUE),
      .groups = 'drop'
    )

  if (nrow(exp_data) == 0) {
    return(ggplot() +
             annotate("text", x = 0.5, y = 0.5, label = "No experimental estimates found", size = 6) +
             theme_void())
  }

  # Create a simple table-style summary plot
  p <- exp_data %>%
    select(study_label, final_resistance, selection_coefficient) %>%
    pivot_longer(cols = c(final_resistance, selection_coefficient),
                 names_to = "metric", values_to = "value") %>%
    mutate(
      metric = case_when(
        metric == "final_resistance" ~ "Final Resistance (%)",
        metric == "selection_coefficient" ~ "Selection Coefficient"
      ),
      # Format values appropriately
      formatted_value = case_when(
        metric == "Final Resistance (%)" ~ paste0(round(value * 100, 1), "%"),
        metric == "Selection Coefficient" ~ as.character(round(value, 2))
      )
    ) %>%
    ggplot(aes(x = study_label, y = metric)) +
    geom_tile(aes(fill = value), color = "white", size = 1) +
    geom_text(aes(label = formatted_value), color = "white", fontface = "bold", size = 4) +
    scale_fill_viridis_c(option = "plasma", name = "Value") +
    labs(
      title = "Experimental Estimates Summary",
      subtitle = "Average outcomes across transmission settings",
      x = "Experimental scenario",
      y = "Outcome metric"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 12),
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title.x = element_blank(),
      panel.grid = element_blank(),
      legend.position = "none"
    )

  return(p)
}



# Continued plot functions ------------------------------------------------

plot_model_demonstration <- function(demo_data) {
  demo_data$scenario <- factor(demo_data$scenario,
                               levels = c("No advantage",
                                          "Moderate advantage (2x)",
                                          "Strong advantage (4x)"))

  ggplot(demo_data, aes(x = time_years, y = resistance_prevalence,
                        color = scenario, linetype = scenario)) +
    geom_line(size = 2) +
    geom_vline(xintercept = 1, linetype = "dashed", alpha = 0.7, color = "black", size = 1) +
    scale_color_manual(values = setNames(unique(demo_data$color), unique(demo_data$scenario))) +
    scale_linetype_manual(values = c("solid", "longdash", "dotted")) +
    labs(
      title = "Model Demonstration: Impact of Transmission Advantages on Resistance Spread",
      subtitle = "EIR = 40 | Treatment coverage = 34% | Resistance advantages active after year 1",
      x = "Time (years)",
      y = "Resistant infection prevalence (%)",
      color = "Infectiousness advantage",
      linetype = "Infectiousness advantage"
    ) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      plot.title = element_text(face = "bold", size = 14)
    )
}

plot_simple_timeseries <- function(timeseries_data) {
  ggplot(timeseries_data, aes(x = time_years, y = resistance_prevalence,
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
      linetype = "Study type"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      legend.position = "bottom"
    )
}



run_integrated_gametocyte_analysis <- function() {

  # Results storage
  all_results <- list()

  # Generate all datasets
  demo_data <- create_model_demonstration_data()

  # Generating timeseries comparison data
  timeseries_data <- create_simple_timeseries_data()

  # Generating transmission mechanism data
  mechanism_data <- create_transmission_mechanism_data()

  # Generating parameter sweep data
  sweep_data <- create_parameter_sweep_data()

  #Generating duration analysis data
  duration_data <- create_duration_analysis_data()

  # Generating comprehensive heatmap data
  heatmap_data <- create_comprehensive_heatmap_data()

  # Running comprehensive parameter analysis
  comprehensive_data <- comprehensive_parameter_analysis()

  # Creating and saving all plots
  p1 <- plot_model_demonstration(demo_data)
  ggsave("01_model_demonstration_integrated.png", p1, width = 12, height = 8, dpi = 300)

  p2 <- plot_simple_timeseries(timeseries_data)
  ggsave("02_timeseries_comparison_integrated.png", p2, width = 12, height = 8, dpi = 300)

  p3 <- plot_infection_prevalence(mechanism_data)
  ggsave("03_infection_prevalence_mechanism.png", p3, width = 12, height = 8, dpi = 300)

  p4 <- plot_transmission_intensity(mechanism_data)
  ggsave("04_transmission_intensity_strain.png", p4, width = 12, height = 8, dpi = 300)

  p5 <- plot_parameter_sweep(sweep_data)
  ggsave("05_parameter_sweep_ratios.png", p5, width = 12, height = 8, dpi = 300)

  p6 <- plot_duration_analysis(duration_data)
  ggsave("06_duration_analysis.png", p6, width = 12, height = 8, dpi = 300)

  p7 <- plot_comprehensive_heatmap_grid(heatmap_data)
  ggsave("07_comprehensive_heatmap_grid.png", p7, width = 14, height = 10, dpi = 300, bg = "white")

  p8 <- plot_detailed_heatmap_clean(comprehensive_data)
  ggsave("08_detailed_heatmap_clean.png", p8, width = 16, height = 12, dpi = 300, bg = "white")

  p9 <- plot_selection_landscape_clean(comprehensive_data)
  ggsave("09_selection_landscape_clean.png", p9, width = 16, height = 12, dpi = 300, bg = "white")

  p10 <- plot_experimental_summary_clean(comprehensive_data)
  ggsave("10_experimental_summary_clean.png", p10, width = 12, height = 8, dpi = 300, bg = "white")

  # Saving all datasets
  write.csv(demo_data, "demo_data_integrated.csv", row.names = FALSE)
  write.csv(timeseries_data, "timeseries_data_integrated.csv", row.names = FALSE)
  write.csv(mechanism_data, "mechanism_data_integrated.csv", row.names = FALSE)
  write.csv(sweep_data, "sweep_data_integrated.csv", row.names = FALSE)
  write.csv(duration_data, "duration_data_integrated.csv", row.names = FALSE)
  write.csv(heatmap_data, "heatmap_data_integrated.csv", row.names = FALSE)
  write.csv(comprehensive_data, "comprehensive_data_integrated.csv", row.names = FALSE)

  # Generating summary statistics
  exp_summary <- comprehensive_data %>%
    filter(is_experimental_estimate == TRUE) %>%
    group_by(parameter_type) %>%
    summarise(
      n_points = n(),
      mean_final_resistance = mean(final_resistance, na.rm = TRUE),
      mean_selection_coeff = mean(selection_coefficient, na.rm = TRUE),
      median_final_resistance = median(final_resistance, na.rm = TRUE),
      median_selection_coeff = median(selection_coefficient, na.rm = TRUE),
      .groups = 'drop'
    )

  print(exp_summary)

  # Overall parameter space coverage
  param_summary <- comprehensive_data %>%
    summarise(
      total_combinations = n(),
      experimental_estimates = sum(is_experimental_estimate),
      mean_final_resistance = mean(final_resistance, na.rm = TRUE),
      mean_selection_coeff = mean(selection_coefficient, na.rm = TRUE),
      resistance_range = paste(round(range(final_resistance, na.rm = TRUE), 3), collapse = " - "),
      selection_range = paste(round(range(selection_coefficient, na.rm = TRUE), 2), collapse = " - ")
    )

  print(param_summary)

  # Store results
  all_results$plots <- list(
    model_demo = p1,
    timeseries = p2,
    infection_prevalence = p3,
    transmission_intensity = p4,
    parameter_sweep = p5,
    duration_analysis = p6,
    heatmap_grid = p7,
    heatmap_clean = p8,
    selection_clean = p9,
    summary_clean = p10
  )

  all_results$data <- list(
    demo = demo_data,
    timeseries = timeseries_data,
    mechanism = mechanism_data,
    sweep = sweep_data,
    duration = duration_data,
    heatmap_grid = heatmap_data,
    comprehensive = comprehensive_data
  )

  all_results$summaries <- list(
    experimental = exp_summary,
    overall = param_summary
  )

  return(all_results)
}



# Run analysis! -----------------------------------------------------------

results <- run_integrated_gametocyte_analysis()



# Getting selection coefficients for baseline  ----------------------------


get_baseline_experimental_results <- function() {
  estimates <- define_experimental_ranges()

  scenarios <- list(
    "In vivo conservative" = list(baseline = 0.903, treated = 0.903),
    "In vivo central" = list(baseline = 1.265, treated = 1.265),
    "In vivo optimistic" = list(baseline = 1.759, treated = 1.759),
    "In vitro conservative" = list(baseline = 0.8, treated = 7.0),
    "In vitro central" = list(baseline = 1.0, treated = 7.5),
    "In vitro optimistic" = list(baseline = 1.0, treated = 8.0)
  )

  results_list <- list()

  for (scenario_name in names(scenarios)) {
    params <- scenarios[[scenario_name]]

    model <- malaria_model(
      EIR = 40, ft = 0.34,
      ton = 365, toff = 365 + (5*365),
      day0_res = 0.01,
      treatment_failure_rate = 0.43,
      rT_r_cleared = 0.1, rT_r_failed = 0.1,
      resistance_baseline_ratio = params$baseline,
      resistance_cleared_ratio = params$treated,
      resistance_failed_ratio = params$treated
    )

    times <- seq(0, 365 + (5*365), by = 60)
    output <- model$run(times)

    final_resistance <- output[nrow(output), "prevalence_res"]

    # Calculate selection coefficient
    p0_idx <- which.min(abs(output[, "t"] - 395))
    p1_idx <- which.min(abs(output[, "t"] - 730))

    p0 <- max(min(output[p0_idx, "prevalence_res"], 0.999), 0.001)
    p1 <- max(min(output[p1_idx, "prevalence_res"], 0.999), 0.001)
    sel_coeff <- log(p1*(1-p0)/(p0*(1-p1)))

    results_list[[scenario_name]] <- data.frame(
      scenario = scenario_name,
      selection_coefficient = sel_coeff,
      final_resistance = final_resistance
    )
  }

  return(do.call(rbind, results_list))
}

baseline_results <- get_baseline_experimental_results()
print(baseline_results)



# Broader ranges ----------------------------------------------------------


comprehensive_data_broad <- comprehensive_parameter_analysis()

# Updated summary plot
p10_broad <- plot_experimental_summary_clean(comprehensive_data_broad)
ggsave("10_experimental_summary_broad_ranges.png", p10_broad,
       width = 12, height = 8, dpi = 300, bg = "white")


# Kb Comparison Plots (CHECK AS DOESN'T SEEM RIGHT) ---------------------------------------------

generate_single_kb_scenario <- function(scenario_name = "current_model") {
  cat(sprintf("=== GENERATING %s SCENARIO ===\n", toupper(scenario_name)))

  # Parameters for analysis
  EIR_values <- c(10, 25, 50, 100)
  ft_values <- c(0.3, 0.5, 0.7)

  # In vivo experimental estimates only
  estimates <- data.frame(
    scenario = c("In vivo conservative", "In vivo central", "In vivo optimistic"),
    baseline_ratio = c(0.903, 1.265, 1.759),
    treated_ratio = c(0.903, 1.265, 1.759)
  )

  results_list <- list()

  # Running analysis for each combination
  for (eir in EIR_values) {
    for (ft in ft_values) {
      for (i in 1:nrow(estimates)) {
        est <- estimates[i, ]


        # Creating model with current parameters
        model <- malaria_model(
          EIR = eir, ft = ft,
          ton = 365, toff = 365 + (5*365),
          day0_res = 0.01,
          treatment_failure_rate = 0.43,
          rT_r_cleared = 0.1, rT_r_failed = 0.1,
          resistance_baseline_ratio = est$baseline_ratio,
          resistance_cleared_ratio = est$treated_ratio,
          resistance_failed_ratio = est$treated_ratio
        )

        # Run model
        times <- seq(0, 365 + (5*365), by = 60)
        output <- model$run(times)

        # Calculating selection coefficient
        p0_idx <- which.min(abs(output[, "t"] - 395))  # ~1 year
        p1_idx <- which.min(abs(output[, "t"] - 730))  # ~2 years

        p0 <- max(min(output[p0_idx, "prevalence_res"], 0.999), 0.001)
        p1 <- max(min(output[p1_idx, "prevalence_res"], 0.999), 0.001)

        sel_coeff <- log(p1*(1-p0)/(p0*(1-p1)))
        final_res <- output[nrow(output), "prevalence_res"]

        # Store results
        results_list[[length(results_list) + 1]] <- data.frame(
          EIR = eir,
          ft = ft,
          scenario = est$scenario,
          baseline_ratio = est$baseline_ratio,
          treated_ratio = est$treated_ratio,
          selection_coefficient = sel_coeff,
          final_resistance = final_res,
          model_version = scenario_name
        )
      }
    }
  }

  results_df <- do.call(rbind, results_list)

  # Print summary statistics
  summary_stats <- results_df %>%
    summarise(
      mean_selection = mean(selection_coefficient),
      sd_selection = sd(selection_coefficient),
      min_selection = min(selection_coefficient),
      max_selection = max(selection_coefficient),
      mean_final_resistance = mean(final_resistance),
      .groups = 'drop'
    )
  print(summary_stats)

  return(results_df)
}


# PLOTTING FUNCTIONS

create_histogram_plot <- function(results_df, scenario_name) {
  p <- ggplot(results_df, aes(x = selection_coefficient)) +
    geom_histogram(fill = "#2E86AB", alpha = 0.7, bins = 15, color = "white", size = 0.5) +
    labs(
      title = paste("Selection Coefficient Distribution:", scenario_name),
      subtitle = "In vivo estimates across all EIR/treatment combinations",
      x = "Selection coefficient (per year)",
      y = "Frequency"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 12)
    )

  return(p)
}

create_boxplot_plot <- function(results_df, scenario_name) {
  p <- ggplot(results_df, aes(x = scenario, y = selection_coefficient)) +
    geom_boxplot(fill = "#2E86AB", alpha = 0.7, color = "black") +
    geom_jitter(width = 0.2, alpha = 0.6, size = 2) +
    labs(
      title = paste("Selection Coefficients by Estimate:", scenario_name),
      subtitle = "Distribution across all EIR/treatment combinations",
      x = "In vivo estimate",
      y = "Selection coefficient (per year)"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 12),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )

  return(p)
}

create_scatter_plot <- function(results_df, scenario_name) {
  summary_data <- results_df %>%
    group_by(scenario) %>%
    summarise(
      mean_selection = mean(selection_coefficient),
      mean_final_resistance = mean(final_resistance),
      .groups = 'drop'
    )

  p <- ggplot(summary_data, aes(x = mean_selection, y = mean_final_resistance)) +
    geom_point(size = 4, alpha = 0.8, color = "#2E86AB") +
    geom_text(aes(label = gsub("In vivo ", "", scenario)),
              vjust = -0.7, hjust = 0.5, size = 3.5, fontface = "bold") +
    scale_y_continuous(labels = scales::percent) +
    labs(
      title = paste("Selection vs Final Resistance:", scenario_name),
      subtitle = "Average across all EIR/treatment combinations",
      x = "Selection coefficient (per year)",
      y = "Final resistance prevalence"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 12)
    )

  return(p)
}


# MAIN EXECUTION FUNCTIONS

# Run this BEFORE modifying model.R, so file should still say cA_r <- cA * if (t > ton && t < toff) 1 else 1
# This means κB only applies to clinical resistant states (Dr, Tr_cleared, Tr_failed)

run_clinical_only_analysis <- function() {

  # Generate results
  results <- generate_single_kb_scenario("Clinical_States_Only")

  # Create plots
  p1 <- create_histogram_plot(results, "Clinical States Only")
  p2 <- create_boxplot_plot(results, "Clinical States Only")
  p3 <- create_scatter_plot(results, "Clinical States Only")

  # Save plots
  ggsave("kb_clinical_only_histogram.png", p1, width = 10, height = 6, dpi = 300)
  ggsave("kb_clinical_only_boxplot.png", p2, width = 10, height = 6, dpi = 300)
  ggsave("kb_clinical_only_scatter.png", p3, width = 10, height = 6, dpi = 300)

  # Save results
  write.csv(results, "kb_clinical_only_results.csv", row.names = FALSE)

  return(results)
}

# Run clinical analysis
run_clinical_only_analysis()


# Now modify code, should be cA_r <- cA * if (t > ton && t < toff) resistance_baseline_ratio else 1
# Save file and recompile
# This means κB applies to ALL resistant states (Ar, Dr, Tr_cleared, Tr_failed)

devtools::load_all()

run_all_resistant_analysis <- function() {

  # Generate results
  results <- generate_single_kb_scenario("All_Resistant_States")

  # Create plots
  p1 <- create_histogram_plot(results, "All Resistant States")
  p2 <- create_boxplot_plot(results, "All Resistant States")
  p3 <- create_scatter_plot(results, "All Resistant States")

  # Save plots
  ggsave("kb_all_resistant_histogram.png", p1, width = 10, height = 6, dpi = 300)
  ggsave("kb_all_resistant_boxplot.png", p2, width = 10, height = 6, dpi = 300)
  ggsave("kb_all_resistant_scatter.png", p3, width = 10, height = 6, dpi = 300)

  # Save results
  write.csv(results, "kb_all_resistant_results.csv", row.names = FALSE)

  # Compare with clinical only results
  if (file.exists("kb_clinical_only_results.csv")) {
    clinical_results <- read.csv("kb_clinical_only_results.csv")

    cat("=== COMPARISON SUMMARY ===\n")
    cat("Clinical states only - mean selection coefficient:", round(mean(clinical_results$selection_coefficient), 3), "\n")
    cat("All resistant states - mean selection coefficient:", round(mean(results$selection_coefficient), 3), "\n")

    difference <- mean(results$selection_coefficient) - mean(clinical_results$selection_coefficient)
    cat("Difference:", round(difference, 3), "\n")

  }

  return(results)
}

# Run all resistant analysis
run_all_resistant_analysis()
