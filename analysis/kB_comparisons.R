devtools::clean_dll()
devtools::document()
odin::odin_package(".")
devtools::load_all(".")


# kB comparison plot generation -------------------------------------------

library(tidyverse)
library(ggplot2)
library(patchwork)

# Creating analysis function for Uganda district data
run_uganda_district_analysis <- function(kb_scenario = "clinical_only") {

  cat(sprintf("=== RUNNING %s SCENARIO FOR UGANDA DISTRICTS ===\n", toupper(kb_scenario)))

  results_list <- list()

  # Loop through each district
  for(i in 1:nrow(uga_param_data)) {
    district_data <- uga_param_data[i, ]
    district_name <- district_data$district

    cat(sprintf("Processing district: %s\n", district_name))

    # Test both high and low parameter combinations for each district
    scenarios <- list(
      list(eir = district_data$eir_high, ft = district_data$ft_high, label = "high"),
      list(eir = district_data$eir_low, ft = district_data$ft_low, label = "low")
    )

    for(scenario in scenarios) {
      # In vivo experimental estimates (from your original code)
      estimates <- data.frame(
        scenario = c("In vivo conservative", "In vivo central", "In vivo optimistic"),
        baseline_ratio = c(0.903, 1.265, 1.759),
        treated_ratio = c(0.903, 1.265, 1.759)
      )

      for(j in 1:nrow(estimates)) {
        est <- estimates[j, ]

        # Create model with district-specific parameters
        model <- malaria_model(
          EIR = scenario$eir,
          ft = scenario$ft,
          ton = 365,
          toff = 365 + (district_data$years * 365),
          day0_res = district_data$f1,  # Use district-specific starting resistance
          treatment_failure_rate = 0.43,
          rT_r_cleared = 0.1,
          rT_r_failed = 0.1,
          resistance_baseline_ratio = est$baseline_ratio,
          resistance_cleared_ratio = est$treated_ratio,
          resistance_failed_ratio = est$treated_ratio
        )

        # Run model for the duration observed in that district + some extra years
        total_years <- district_data$years + 5
        times <- seq(0, 365 + (total_years * 365), by = 30)
        output <- model$run(times)

        # Calculate selection coefficient using the actual years of observation
        # Use 1 year after start and 2 years after start for consistency
        p0_idx <- which.min(abs(output[, "t"] - (365 + 365)))  # ~2 years from start
        p1_idx <- which.min(abs(output[, "t"] - (365 + (2*365))))  # ~3 years from start

        p0 <- max(min(output[p0_idx, "prevalence_res"], 0.999), 0.001)
        p1 <- max(min(output[p1_idx, "prevalence_res"], 0.999), 0.001)

        sel_coeff <- log(p1*(1-p0)/(p0*(1-p1)))
        final_res <- output[nrow(output), "prevalence_res"]

        # Store results
        results_list[[length(results_list) + 1]] <- data.frame(
          district = district_name,
          EIR = scenario$eir,
          ft = scenario$ft,
          scenario_type = scenario$label,
          vivo_scenario = est$scenario,
          baseline_ratio = est$baseline_ratio,
          treated_ratio = est$treated_ratio,
          selection_coefficient = sel_coeff,
          final_resistance = final_res,
          years_observed = district_data$years,
          starting_freq = district_data$f1,
          kb_scenario = kb_scenario
        )
      }
    }
  }

  results_df <- do.call(rbind, results_list)

  # Print summary statistics
  summary_stats <- results_df %>%
    summarise(
      mean_selection = mean(selection_coefficient, na.rm = TRUE),
      sd_selection = sd(selection_coefficient, na.rm = TRUE),
      min_selection = min(selection_coefficient, na.rm = TRUE),
      max_selection = max(selection_coefficient, na.rm = TRUE),
      mean_final_resistance = mean(final_resistance, na.rm = TRUE),
      .groups = 'drop'
    )

  cat("Summary Statistics:\n")
  print(summary_stats)

  return(results_df)
}

# Function to test which district shows greatest difference
test_district_differences <- function() {

  cat("=== TESTING DISTRICTS FOR GREATEST κB DIFFERENCE ===\n")

  # Test top 5 districts with highest EIR
  test_districts <- c("Agago", "Katakwi", "Kole", "Lamwo", "Kaabong")

  differences <- data.frame(
    district = character(),
    clinical_fixation = numeric(),
    all_resistant_fixation = numeric(),
    difference_years = numeric(),
    clinical_final = numeric(),
    all_resistant_final = numeric()
  )

  baseline_ratio <- 1.265
  treated_ratio <- 1.265

  for(district_name in test_districts) {
    district_data <- uga_param_data %>% filter(district == district_name)

    if(nrow(district_data) == 0) next

    cat(sprintf("Testing %s (EIR: %.1f, ft: %.3f, f1: %.4f)...\n",
                district_name, district_data$eir_high, district_data$ft_high, district_data$f1))

    tryCatch({
      # Clinical only scenario
      model_clinical <- malaria_model(
        EIR = district_data$eir_high,
        ft = district_data$ft_high,
        ton = 365,
        toff = 365 + (30 * 365),  # Extended to 30 years
        day0_res = district_data$f1,
        treatment_failure_rate = 0.43,
        rT_r_cleared = 0.1,
        rT_r_failed = 0.1,
        resistance_baseline_ratio = 1.0,  # No effect on Ar
        resistance_cleared_ratio = treated_ratio,
        resistance_failed_ratio = treated_ratio
      )

      # All resistant scenario
      model_all <- malaria_model(
        EIR = district_data$eir_high,
        ft = district_data$ft_high,
        ton = 365,
        toff = 365 + (30 * 365),  # Extended to 30 years
        day0_res = district_data$f1,
        treatment_failure_rate = 0.43,
        rT_r_cleared = 0.1,
        rT_r_failed = 0.1,
        resistance_baseline_ratio = baseline_ratio,  # Effect on Ar too
        resistance_cleared_ratio = treated_ratio,
        resistance_failed_ratio = treated_ratio
      )

      # Extended time horizon - test up to 200 years
      times <- seq(0, 365 * 200, by = 365 * 2)  # Every 2 years to speed up

      # Run models
      output_clinical <- model_clinical$run(times)
      output_all <- model_all$run(times)

      # Get final resistance levels
      clinical_final <- tail(output_clinical[, "prevalence_res"], 1)
      all_final <- tail(output_all[, "prevalence_res"], 1)

      cat(sprintf("  Final resistance - Clinical: %.3f, All resistant: %.3f\n",
                  clinical_final, all_final))

      # Find fixation times (95%) - if they exist
      clinical_fixation_idx <- which(output_clinical[, "prevalence_res"] >= 0.95)[1]
      all_fixation_idx <- which(output_all[, "prevalence_res"] >= 0.95)[1]

      clinical_years <- if(!is.na(clinical_fixation_idx)) times[clinical_fixation_idx] / 365 else NA
      all_years <- if(!is.na(all_fixation_idx)) times[all_fixation_idx] / 365 else NA

      # Calculate difference based on what we have
      if(!is.na(clinical_years) && !is.na(all_years)) {
        diff_years <- clinical_years - all_years
      } else {
        # Use final resistance levels as proxy for speed
        diff_years <- all_final - clinical_final  # Higher = faster resistance growth
      }

      differences <- rbind(differences, data.frame(
        district = district_name,
        clinical_fixation = ifelse(is.na(clinical_years), 999, clinical_years),
        all_resistant_fixation = ifelse(is.na(all_years), 999, all_years),
        difference_years = diff_years,
        clinical_final = clinical_final,
        all_resistant_final = all_final
      ))

    }, error = function(e) {
      cat(sprintf("  Error with %s: %s\n", district_name, e$message))
    })
  }

  if(nrow(differences) == 0) {
    cat("No districts successfully tested. Using Agago as default.\n")
    return("Agago")
  }

  # Sort by greatest difference (positive = clinical slower, negative = difference in final resistance)
  differences <- differences[order(-differences$difference_years), ]

  cat("\n=== DISTRICT COMPARISON RESULTS ===\n")
  print(differences)

  best_district <- differences$district[1]
  cat(sprintf("\nBest district for showing difference: %s\n", best_district))

  return(best_district)
}

# Function to create individual fixation plots (updated to use best district)
create_individual_fixation_plots <- function(district_name = NULL) {

  cat("=== CREATING INDIVIDUAL FIXATION PLOTS ===\n")

  # Use provided district or find the best one
  if(is.null(district_name)) {
    district_name <- test_district_differences()
  }

  # Handle case where no district was found
  if(is.na(district_name) || length(district_name) == 0) {
    cat("No suitable district found. Using Agago as default.\n")
    district_name <- "Agago"
  }

  # Select the specified district
  district_data <- uga_param_data %>% filter(district == district_name)

  if(nrow(district_data) == 0) {
    cat("District not found. Using Agago as fallback.\n")
    district_name <- "Agago"
    district_data <- uga_param_data %>% filter(district == district_name)
  }

  cat(sprintf("Using district: %s\n", district_name))

  # Use central estimate
  baseline_ratio <- 1.265
  treated_ratio <- 1.265

  # PLOT 1: Clinical States Only (κB doesn't affect Ar)

  tryCatch({
    model_clinical <- malaria_model(
      EIR = district_data$eir_high,
      ft = district_data$ft_high,
      ton = 365,
      toff = 365 + (50 * 365),  # 50 years of intervention
      day0_res = district_data$f1,
      treatment_failure_rate = 0.43,
      rT_r_cleared = 0.1,
      rT_r_failed = 0.1,
      resistance_baseline_ratio = 1.0,  # No effect on asymptomatic (Ar)
      resistance_cleared_ratio = treated_ratio,  # Only clinical get boost
      resistance_failed_ratio = treated_ratio
    )

    # Run for 200 years but plot first 100
    times <- seq(0, 365 * 200, by = 365)
    output_clinical <- model_clinical$run(times)

    df_clinical <- data.frame(
      time_years = output_clinical[, "t"] / 365,
      resistance_prevalence = output_clinical[, "prevalence_res"]
    )

    # Limit to first 100 years for plotting
    df_clinical <- df_clinical[df_clinical$time_years <= 100, ]

    # Find fixation time for clinical only
    fixation_clinical <- df_clinical %>%
      filter(resistance_prevalence >= 0.95) %>%
      slice(1) %>%
      pull(time_years)

    if(length(fixation_clinical) == 0) fixation_clinical <- NA

    plot_clinical <- ggplot(df_clinical, aes(x = time_years, y = resistance_prevalence)) +
      geom_line(color = "#E31A1C", size = 1.2) +
      scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
      scale_x_continuous(limits = c(0, 100)) +
      geom_hline(yintercept = 0.95, linetype = "dashed", color = "gray50", alpha = 0.7) +
      labs(
        title = "κB: Clinical States Only",
        subtitle = ifelse(!is.na(fixation_clinical),
                          paste0("Fixation (95%) at ~", round(fixation_clinical, 1), " years"),
                          paste0("Final resistance at 100y: ", round(tail(df_clinical$resistance_prevalence, 1) * 100, 1), "%")),
        x = "Time (years)",
        y = "Resistance Prevalence"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(face = "bold", size = 14, color = "#E31A1C"),
        plot.subtitle = element_text(size = 12)
      )

  }, error = function(e) {
    cat("Error creating clinical plot:", e$message, "\n")
    plot_clinical <- ggplot() + labs(title = "Error in clinical model")
  })

  # PLOT 2: All Resistant States (κB affects Ar + clinical)

  tryCatch({
    model_all <- malaria_model(
      EIR = district_data$eir_high,
      ft = district_data$ft_high,
      ton = 365,
      toff = 365 + (50 * 365),
      day0_res = district_data$f1,
      treatment_failure_rate = 0.43,
      rT_r_cleared = 0.1,
      rT_r_failed = 0.1,
      resistance_baseline_ratio = baseline_ratio,  # Effect on asymptomatic (Ar) too
      resistance_cleared_ratio = treated_ratio,
      resistance_failed_ratio = treated_ratio
    )

    times <- seq(0, 365 * 200, by = 365)
    output_all <- model_all$run(times)

    df_all <- data.frame(
      time_years = output_all[, "t"] / 365,
      resistance_prevalence = output_all[, "prevalence_res"]
    )

    # Limit to first 100 years for plotting
    df_all <- df_all[df_all$time_years <= 100, ]

    # Find fixation time for all resistant
    fixation_all <- df_all %>%
      filter(resistance_prevalence >= 0.95) %>%
      slice(1) %>%
      pull(time_years)

    if(length(fixation_all) == 0) fixation_all <- NA

    plot_all <- ggplot(df_all, aes(x = time_years, y = resistance_prevalence)) +
      geom_line(color = "#1F78B4", size = 1.2) +
      scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
      scale_x_continuous(limits = c(0, 100)) +
      geom_hline(yintercept = 0.95, linetype = "dashed", color = "gray50", alpha = 0.7) +
      labs(
        title = "κB: All Resistant States",
        subtitle = ifelse(!is.na(fixation_all),
                          paste0("Fixation (95%) at ~", round(fixation_all, 1), " years"),
                          paste0("Final resistance at 100y: ", round(tail(df_all$resistance_prevalence, 1) * 100, 1), "%")),
        x = "Time (years)",
        y = "Resistance Prevalence"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(face = "bold", size = 14, color = "#1F78B4"),
        plot.subtitle = element_text(size = 12)
      )

  }, error = function(e) {
    cat("Error creating all resistant plot:", e$message, "\n")
    plot_all <- ggplot() + labs(title = "Error in all resistant model")
  })

  # Print comparison
  cat("District:", district_name, "\n")
  cat("EIR:", district_data$eir_high, "ft:", district_data$ft_high, "Starting freq:", district_data$f1, "\n")

  return(list(clinical_plot = plot_clinical, all_resistant_plot = plot_all))
}

# Function to create comparison plots
create_uganda_comparison_plots <- function(clinical_results, all_resistant_results) {

  # Combine results for comparison
  combined_results <- rbind(clinical_results, all_resistant_results)

  # 1. Selection coefficient comparison by district
  p1 <- ggplot(combined_results, aes(x = district, y = selection_coefficient, fill = kb_scenario)) +
    geom_boxplot(alpha = 0.7) +
    labs(
      title = "Selection Coefficients by District: κB Scenario Comparison",
      subtitle = "Across all EIR/ft combinations and in vivo estimates",
      x = "District",
      y = "Selection coefficient (per year)",
      fill = "κB Scenario"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom"
    ) +
    scale_fill_manual(values = c("clinical_only" = "#E31A1C", "all_resistant" = "#1F78B4"),
                      labels = c("Clinical States Only", "All Resistant States"))

  # 2. Overall selection coefficient comparison
  p2 <- ggplot(combined_results, aes(x = kb_scenario, y = selection_coefficient, fill = kb_scenario)) +
    geom_violin(alpha = 0.7, width = 0.8) +
    geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
    stat_summary(fun = mean, geom = "point", size = 3, color = "red") +
    labs(
      title = "Overall Selection Coefficient Distribution",
      subtitle = "Comparison between κB scenarios across all districts",
      x = "κB Scenario",
      y = "Selection coefficient (per year)"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      legend.position = "none"
    ) +
    scale_x_discrete(labels = c("Clinical States\nOnly", "All Resistant\nStates")) +
    scale_fill_manual(values = c("clinical_only" = "#E31A1C", "all_resistant" = "#1F78B4"))

  # 3. Comparison with Meier-Scherling estimates
  # Individual mutation estimates
  individual_estimates <- data.frame(
    mutation = c("Pro441Leu", "Cys469Phe", "Cys469Tyr", "Ala675Val"),
    estimate = c(0.494, 0.324, 0.383, 0.237),
    lower = c(-0.462, -0.629, 0.207, 0.087),
    upper = c(1.410, 1.150, 0.591, 0.403),
    type = "Individual",
    order = 4:7
  )

  # Uganda Combined estimate
  uganda_combined <- data.frame(
    mutation = "Uganda Combined",
    estimate = 0.381,
    lower = 0.298,
    upper = 0.472,
    type = "Uganda Combined",
    order = 3
  )

  # Summary statistics for our model results
  model_summary <- combined_results %>%
    group_by(kb_scenario) %>%
    summarise(
      mean_sel = mean(selection_coefficient, na.rm = TRUE),
      lower_ci = quantile(selection_coefficient, 0.025, na.rm = TRUE),
      upper_ci = quantile(selection_coefficient, 0.975, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    mutate(
      mutation = paste("Model Combined\n",
                       ifelse(kb_scenario == "clinical_only", "Clinical Only", "All Resistant")),
      type = "Model",
      order = ifelse(kb_scenario == "clinical_only", 1, 2)
    )

  # Combine Meier-Scherling data
  meier_data <- rbind(individual_estimates, uganda_combined)

  p3 <- ggplot() +
    # Model estimates (colored by scenario)
    geom_errorbar(data = model_summary,
                  aes(x = reorder(mutation, order), ymin = lower_ci, ymax = upper_ci, color = kb_scenario),
                  width = 0.2, size = 1) +
    geom_point(data = model_summary,
               aes(x = reorder(mutation, order), y = mean_sel, color = kb_scenario),
               size = 3) +
    # Meier-Scherling estimates
    geom_errorbar(data = meier_data,
                  aes(x = reorder(mutation, order), ymin = lower, ymax = upper,
                      color = type),
                  width = 0.2, size = 1) +
    geom_point(data = meier_data,
               aes(x = reorder(mutation, order), y = estimate, color = type),
               size = 3) +
    labs(
      title = "Model vs. Meier-Scherling Selection Coefficients",
      subtitle = "Comparison of model predictions with Uganda field estimates",
      x = "Estimate Type",
      y = "Selection coefficient (per year)",
      color = "Estimate Source"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom"
    ) +
    scale_color_manual(
      values = c(
        "clinical_only" = "#E31A1C",
        "all_resistant" = "#1F78B4",
        "Uganda Combined" = "#FF7F00",  # Orange for Uganda combined
        "Individual" = "#228B22"       # Dark green for individual mutations
      ),
      labels = c(
        "clinical_only" = "Model: Clinical States Only",
        "all_resistant" = "Model: All Resistant States",
        "Uganda Combined" = "Uganda: Combined Estimate",
        "Individual" = "Uganda: Individual Mutations"
      ),
      breaks = c("clinical_only", "all_resistant", "Uganda Combined", "Individual")
    )

  return(list(district_plot = p1, overall_plot = p2, comparison_plot = p3))
}

# MAIN EXECUTION

# STEP 1: Run analysis for clinical states only scenario
# Make sure model.R file has: cA_r <- cA * if (t > ton && t < toff) 1 else 1
clinical_results <- run_uganda_district_analysis("clinical_only")

# STEP 2: Modify your model and recompile
# Change: cA_r <- cA * if (t > ton && t < toff) 1 else 1
# To: cA_r <- cA * if (t > ton && t < toff) resistance_baseline_ratio else 1

devtools::clean_dll()
devtools::document()
odin::odin_package(".")
devtools::load_all(".")

# STEP 3: Run analysis for all resistant states scenario
all_resistant_results <- run_uganda_district_analysis("all_resistant")

# STEP 4: Create comparison plots
plots <- create_uganda_comparison_plots(clinical_results, all_resistant_results)

# STEP 5: Create fixation comparison plots (individual)
fixation_plots <- create_individual_fixation_plots()

# Save individual fixation plots
ggsave("uganda_fixation_clinical_only.png", fixation_plots$clinical_plot, width = 10, height = 8, dpi = 300)
ggsave("uganda_fixation_all_resistant.png", fixation_plots$all_resistant_plot, width = 10, height = 8, dpi = 300)

# Side-by-side comparison using patchwork
if(require(patchwork, quietly = TRUE)) {
  combined_fixation <- fixation_plots$clinical_plot + fixation_plots$all_resistant_plot +
    plot_annotation(
      title = "Resistance Fixation Comparison: κB Scenarios",
      subtitle = "Representative district (Tororo) - Shows impact of κB on asymptomatic infections"
    )

  ggsave("uganda_fixation_side_by_side.png", combined_fixation, width = 16, height = 8, dpi = 300)
}

# Save all plots
ggsave("uganda_selection_by_district.png", plots$district_plot, width = 14, height = 8, dpi = 300)
ggsave("uganda_overall_selection_comparison.png", plots$overall_plot, width = 10, height = 8, dpi = 300)
ggsave("uganda_meier_scherling_comparison.png", plots$comparison_plot, width = 12, height = 8, dpi = 300)

# Save results
write.csv(clinical_results, "uganda_clinical_only_results.csv", row.names = FALSE)
write.csv(all_resistant_results, "uganda_all_resistant_results.csv", row.names = FALSE)

# Print final comparison summary
cat("Clinical states only - mean selection coefficient:", round(mean(clinical_results$selection_coefficient, na.rm = TRUE), 3), "\n")
cat("All resistant states - mean selection coefficient:", round(mean(all_resistant_results$selection_coefficient, na.rm = TRUE), 3), "\n")

difference <- mean(all_resistant_results$selection_coefficient, na.rm = TRUE) - mean(clinical_results$selection_coefficient, na.rm = TRUE)
cat("Difference (All resistant - Clinical only):", round(difference, 3), "\n")

# Compare with Meier-Scherling combined estimate (0.381)
meier_combined <- 0.381
cat("Meier-Scherling combined estimate:", meier_combined, "\n")
cat("Clinical only vs Meier-Scherling difference:", round(mean(clinical_results$selection_coefficient, na.rm = TRUE) - meier_combined, 3), "\n")
cat("All resistant vs Meier-Scherling difference:", round(mean(all_resistant_results$selection_coefficient, na.rm = TRUE) - meier_combined, 3), "\n")
