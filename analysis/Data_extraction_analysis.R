# ========================
# K13 RESISTANCE ANALYSIS
# ========================

# Required libraries
library(ggplot2)
library(dplyr)
library(broom)
library(weights)
library(patchwork)
library(readxl)
library(scales)

# Read in data from Excel file
malaria_data <- read_excel("Plot_Data_01.xlsx", sheet = "Sheet1")

#If not, data frame below
malaria_data <- data.frame(
  Study = c(
    "Ladeia-Andrade 2016", "Olivera 2019", "Teklemariam 2015", "Tesfaye 2024",
    "Ippolito 2020", "Getnet 2013", "Atroosh 2014", "Chang 2016", "Chang 2016",
    "Djimde 2016", "Kyaw Myo Tun 2015", "Kyaw Myo Tun 2015", "Vantaux 2020",
    "Vantaux 2020", "Thriemer 2014", "Rovira-Vallbona 2020", "Rovira-Vallbona 2020",
    "Ashley 2014", "Ashley 2014", "Ashley 2014", "Ashley 2014", "Ashley 2014",
    "Ashley 2014", "Ashley 2014", "Carrara 2009", "Lek 2022"
  ),
  k13_mutation_pct = c(
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 46, 46, 63, 63, 80.7, 87, 87, 44.7, 73.6,
    27.5, 26, 0, 1.7, 0, 21.9, 57.9
  ),
  `baseline_gam_prevalence` = c(
    32.5, 12.5, 7.6, 6.3, 13, 10, 40.7, 37, 50, 22.8, 22, 21, 49, 37.5, 18.9,
    NA, 76, 11.8, 6.7, 5.9, 11.3, 0, 30.8, 2.9, 6.3, 46.9
  ),
  baseline_asexual_density = c(
    4700.1, 3527, 27798, 10627, 13000, NA, 8199, NA, NA, NA, 3084, 2140, NA, NA,
    8233, 10166, 15832, 55623, 36529, 49738, 66066, 27130, 60037, 52250, 6982, 11873
  ),
  `baseline_gam_density` = c(
    NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
    NA, NA, NA, NA, NA, NA, NA
  ),
  `day3_gam_prevalence` = c(
    23.6, 18.6, 1, 2.5, 3, 37.5, NA, NA, NA, 17.4, 8, 13, NA, NA, NA, NA, NA,
    NA, NA, NA, NA, NA, NA, NA, NA, NA
  ),
  `day7_gam_prevalence` = c(
    NA, 14.3, NA, 2.5, 2, 0, 34.9, 7.84, 0, 23.9, 0, 1.4, NA, NA, 20.5, NA,
    NA, NA, NA, NA, NA, NA, NA, NA, NA, NA
  ),
  `sample_size` = c(
    162, 84, 92, 73, 94, 80, 89, 61, 9, 92, 78, 76, 53, 56, 95, 33, 27, 460,
    185, 120, 80, 56, 120, 40, 3264, 211
  ),
  stringsAsFactors = FALSE
)

# Cleaning column names
colnames(malaria_data) <- trimws(colnames(malaria_data))

# Renaming
malaria_data <- malaria_data %>%
  rename(
    study_n = `sample_size`,
    baseline_gam_prevalence = `baseline_gam_prevalence`,
    day3_gam_prevalence = `day3_gam_prevalence`,
    day7_gam_prevalence = `day7_gam_prevalence`
  )

# Converting "NA" strings to actual NA values
malaria_data[malaria_data == "NA"] <- NA

numeric_cols <- c("k13_mutation_pct", "baseline_gam_prevalence", "baseline_asexual_density",
                  "baseline_gam_density", "day3_gam_prevalence", "day7_gam_prevalence", "study_n")

for(col in numeric_cols) {
  if(col %in% colnames(malaria_data)) {
    malaria_data[[col]] <- as.numeric(malaria_data[[col]])
  }
}

# Dataset summary
cat("Total study populations:", nrow(malaria_data), "\n")
cat("Columns available:", paste(colnames(malaria_data), collapse = ", "), "\n\n")

# First few rows
print(head(malaria_data))


# Core functions  ---------------------------------------------------------

# Weighted binomial regression function
perform_weighted_binomial <- function(data, outcome_var, predictor_var = "k13_mutation_pct", weight_var = "study_n") {

  # Filter complete cases
  complete_data <- data[!is.na(data[[outcome_var]]) &
                          !is.na(data[[predictor_var]]) &
                          !is.na(data[[weight_var]]), ]

  if(nrow(complete_data) < 3) {
    return(list(error = "Insufficient data"))
  }

  # Converting percentages to proportions and creating binomial counts
  complete_data$successes <- round((complete_data[[outcome_var]]/100) * complete_data[[weight_var]])
  complete_data$failures <- complete_data[[weight_var]] - complete_data$successes

  # Ensuring valid counts
  complete_data$successes <- pmax(0, pmin(complete_data$successes, complete_data[[weight_var]]))
  complete_data$failures <- complete_data[[weight_var]] - complete_data$successes

  cat("\nStudies included in", outcome_var, "analysis:\n")
  for(i in 1:nrow(complete_data)) {
    cat(sprintf("  %s: %.1f%% (n=%d)\n",
                complete_data$Study[i],
                complete_data[[outcome_var]][i],
                complete_data[[weight_var]][i]))
  }


  formula_str <- paste("cbind(successes, failures) ~", predictor_var)
  model <- glm(as.formula(formula_str), family = binomial, data = complete_data)

  # Extracting results
  model_summary <- summary(model)
  coefficients <- model_summary$coefficients

  if(!predictor_var %in% rownames(coefficients)) {
    return(list(error = "Model fitting failed"))
  }

  # Calculate odds ratio and confidence intervals
  beta <- coefficients[predictor_var, "Estimate"]
  se <- coefficients[predictor_var, "Std. Error"]
  odds_ratio <- exp(beta)
  or_ci_lower <- exp(beta - 1.96 * se)
  or_ci_upper <- exp(beta + 1.96 * se)
  p_value <- coefficients[predictor_var, "Pr(>|z|)"]

  # Generating predictions for 0% to 100% K13 resistance
  pred_x <- seq(0, 100, by = 5)

  # Creating prediction data frame
  pred_data <- data.frame(pred_x)
  names(pred_data)[1] <- predictor_var

  pred_logit <- predict(model, newdata = pred_data, se.fit = TRUE)

  predictions <- data.frame(
    k13_resistance = pred_x,
    predicted_prevalence = plogis(pred_logit$fit) * 100,
    se = pred_logit$se.fit
  )
  predictions$ci_lower <- plogis(pred_logit$fit - 1.96 * predictions$se) * 100
  predictions$ci_upper <- plogis(pred_logit$fit + 1.96 * predictions$se) * 100

  # Calculating transmission multipliers
  pred_data_0 <- data.frame(0)
  names(pred_data_0)[1] <- predictor_var
  pred_data_100 <- data.frame(100)
  names(pred_data_100)[1] <- predictor_var

  pred_at_0 <- plogis(predict(model, newdata = pred_data_0)) * 100
  pred_at_100 <- plogis(predict(model, newdata = pred_data_100)) * 100

  # Getting CI bounds for multiplier calculation
  pred_0_ci <- predictions[predictions$k13_resistance == 0, ]
  pred_100_ci <- predictions[predictions$k13_resistance == 100, ]

  multiplier_best <- pred_at_100 / pred_at_0
  multiplier_conservative <- pred_100_ci$ci_lower / pred_0_ci$ci_upper
  multiplier_optimistic <- pred_100_ci$ci_upper / pred_0_ci$ci_lower

  return(list(
    model = model,
    data = complete_data,
    n_studies = nrow(complete_data),
    total_participants = sum(complete_data[[weight_var]]),
    odds_ratio = odds_ratio,
    or_ci_lower = or_ci_lower,
    or_ci_upper = or_ci_upper,
    p_value = p_value,
    predictions = predictions,
    pred_at_0 = pred_at_0,
    pred_at_100 = pred_at_100,
    multiplier_best = multiplier_best,
    multiplier_conservative = multiplier_conservative,
    multiplier_optimistic = multiplier_optimistic,
    aic = AIC(model)
  ))
}

# Weighted correlation function
calculate_weighted_correlation <- function(data, x_var, y_var, weight_var = "study_n") {
  complete_data <- data[!is.na(data[[x_var]]) & !is.na(data[[y_var]]) & !is.na(data[[weight_var]]), ]

  if(nrow(complete_data) < 3) {
    return(list(error = "Insufficient data"))
  }

  # Weighted correlation
  weighted_cor <- wtd.cor(complete_data[[x_var]], complete_data[[y_var]],
                          weight = complete_data[[weight_var]])

  # Unweighted for comparison
  unweighted_cor <- cor.test(complete_data[[x_var]], complete_data[[y_var]])

  return(list(
    n_studies = nrow(complete_data),
    total_participants = sum(complete_data[[weight_var]]),
    weighted_r = weighted_cor[1,1],
    unweighted_r = unweighted_cor$estimate,
    unweighted_p = unweighted_cor$p.value
  ))
}


# Main analyse ------------------------------------------------------------


# 1. Baseline Gametocyte Prevalence Analysis

baseline_analysis <- perform_weighted_binomial(malaria_data, "baseline_gam_prevalence")

if(!"error" %in% names(baseline_analysis)) {
  cat("\nRESULTS:\n")
  cat("Number of studies:", baseline_analysis$n_studies, "\n")
  cat("Total participants:", format(baseline_analysis$total_participants, big.mark = ","), "\n")
  cat("Odds ratio per 1% increase in K13:", round(baseline_analysis$odds_ratio, 4), "\n")
  cat("95% CI:", round(baseline_analysis$or_ci_lower, 4), "to", round(baseline_analysis$or_ci_upper, 4), "\n")
  cat("P-value:", format.pval(baseline_analysis$p_value), "\n")
  cat("AIC:", round(baseline_analysis$aic, 2), "\n")

  cat("\nTRANSMISSION MULTIPLIERS (100% vs 0% K13 resistance):\n")
  cat("Best estimate:", round(baseline_analysis$multiplier_best, 2), "x\n")
  cat("Conservative:", round(baseline_analysis$multiplier_conservative, 2), "x\n")
  cat("Optimistic:", round(baseline_analysis$multiplier_optimistic, 2), "x\n")

  cat("Predicted prevalence at 0% K13:", round(baseline_analysis$pred_at_0, 1), "%\n")
  cat("Predicted prevalence at 100% K13:", round(baseline_analysis$pred_at_100, 1), "%\n")
} else {
  cat("Error:", baseline_analysis$error, "\n")
}

# 2. Day 7 Gametocyte Prevalence Analysis

day7_analysis <- perform_weighted_binomial(malaria_data, "day7_gam_prevalence")

if(!"error" %in% names(day7_analysis)) {
  cat("\nRESULTS:\n")
  cat("Number of studies:", day7_analysis$n_studies, "\n")
  cat("Total participants:", format(day7_analysis$total_participants, big.mark = ","), "\n")
  cat("Odds ratio per 1% increase in K13:", round(day7_analysis$odds_ratio, 4), "\n")
  cat("95% CI:", round(day7_analysis$or_ci_lower, 4), "to", round(day7_analysis$or_ci_upper, 4), "\n")
  cat("P-value:", format.pval(day7_analysis$p_value), "\n")
  cat("AIC:", round(day7_analysis$aic, 2), "\n")

  cat("\nTRANSMISSION MULTIPLIERS (100% vs 0% K13 resistance):\n")
  cat("Best estimate:", round(day7_analysis$multiplier_best, 2), "x\n")
  cat("Conservative:", round(day7_analysis$multiplier_conservative, 2), "x\n")
  cat("Optimistic:", round(day7_analysis$multiplier_optimistic, 2), "x\n")

  cat("Predicted prevalence at 0% K13:", round(day7_analysis$pred_at_0, 1), "%\n")
  cat("Predicted prevalence at 100% K13:", round(day7_analysis$pred_at_100, 1), "%\n")
} else {
  cat("Error:", day7_analysis$error, "\n")
}

# 3. Day 3 Gametocyte Prevalence Analysis (for completeness)

day3_analysis <- perform_weighted_binomial(malaria_data, "day3_gam_prevalence")

if(!"error" %in% names(day3_analysis)) {
  cat("\nRESULTS:\n")
  cat("Number of studies:", day3_analysis$n_studies, "\n")
  cat("Total participants:", format(day3_analysis$total_participants, big.mark = ","), "\n")
  cat("Odds ratio per 1% increase in K13:", round(day3_analysis$odds_ratio, 4), "\n")
  cat("95% CI:", round(day3_analysis$or_ci_lower, 4), "to", round(day3_analysis$or_ci_upper, 4), "\n")
  cat("P-value:", format.pval(day3_analysis$p_value), "\n")
  cat("AIC:", round(day3_analysis$aic, 2), "\n")

  cat("\nEFFECT MULTIPLIER:\n")
  cat("Best estimate:", round(day3_analysis$multiplier_best, 2), "x\n")
  cat("Conservative:", round(day3_analysis$multiplier_conservative, 2), "x\n")
  cat("Optimistic:", round(day3_analysis$multiplier_optimistic, 2), "x\n")
  cat("INTERPRETATION: K13 resistance DECREASES Day 3 gametocytes (faster initial clearance)\n")
} else {
  cat("Error:", day3_analysis$error, "\n")
}


# Weighted correlations ---------------------------------------------------


# Key relationships
correlations <- list(
  baseline = calculate_weighted_correlation(malaria_data, "k13_mutation_pct", "baseline_gam_prevalence"),
  day7 = calculate_weighted_correlation(malaria_data, "k13_mutation_pct", "day7_gam_prevalence"),
  day3 = calculate_weighted_correlation(malaria_data, "k13_mutation_pct", "day3_gam_prevalence")
)

for(name in names(correlations)) {
  corr <- correlations[[name]]
  if(!"error" %in% names(corr)) {
    cat("\nK13 vs", name, "gametocytes:\n")
    cat("  Weighted r =", round(corr$weighted_r, 3), "\n")
    cat("  Unweighted r =", round(corr$unweighted_r, 3), "\n")
    cat("  P-value =", format.pval(corr$unweighted_p), "\n")
    cat("  N studies =", corr$n_studies, ", N participants =", format(corr$total_participants, big.mark = ","), "\n")
  } else {
    cat("\nK13 vs", name, ": Insufficient data\n")
  }
}


# Asexual density analysis ------------------------------------------------


# Filtering data
asexual_data <- malaria_data[!is.na(malaria_data$k13_mutation_pct) &
                               !is.na(malaria_data$baseline_asexual_density) &
                               !is.na(malaria_data$study_n), ]

if(nrow(asexual_data) >= 3) {
  # Log transform
  asexual_data$log_asexual <- log10(asexual_data$baseline_asexual_density)

  cat("Study populations included in asexual density analysis:\n")
  for(i in 1:nrow(asexual_data)) {
    cat(sprintf("  %s: K13=%.1f%%, Density=%s, n=%d\n",
                asexual_data$Study[i],
                asexual_data$k13_mutation_pct[i],
                format(asexual_data$baseline_asexual_density[i], scientific = FALSE, big.mark = ","),
                asexual_data$study_n[i]))
  }

  # Weighted linear regression on log-transformed data
  asexual_model <- lm(log_asexual ~ k13_mutation_pct,
                      data = asexual_data,
                      weights = study_n)

  asexual_summary <- summary(asexual_model)

  # Calculating fold-change from 0% to 100% K13 resistance
  pred_0 <- predict(asexual_model, newdata = data.frame(k13_mutation_pct = 0))
  pred_100 <- predict(asexual_model, newdata = data.frame(k13_mutation_pct = 100))
  fold_change <- 10^pred_100 / 10^pred_0

  cat("\nRESULTS:\n")
  cat("Number of study populations:", nrow(asexual_data), "\n")
  cat("Total participants:", format(sum(asexual_data$study_n), big.mark = ","), "\n")
  cat("Coefficient (log10 scale):", round(asexual_summary$coefficients["k13_mutation_pct", "Estimate"], 6), "\n")
  cat("P-value:", format.pval(asexual_summary$coefficients["k13_mutation_pct", "Pr(>|t|)"]), "\n")
  cat("R-squared:", round(asexual_summary$r.squared, 4), "\n")
  cat("Adjusted R-squared:", round(asexual_summary$adj.r.squared, 4), "\n")

  cat("Fold change in asexual density (0% to 100% K13):", round(fold_change, 2), "x\n")
  cat("Predicted asexual density at 0% K13:", format(round(10^pred_0, 0), big.mark = ","), "parasites/μL\n")
  cat("Predicted asexual density at 100% K13:", format(round(10^pred_100, 0), big.mark = ","), "parasites/μL\n")

  # Prediction data for plotting
  pred_x <- seq(0, 100, by = 5)
  pred_data_asexual <- data.frame(k13_mutation_pct = pred_x)
  pred_log <- predict(asexual_model, newdata = pred_data_asexual, se.fit = TRUE)

  pred_data_asexual$predicted_density <- 10^pred_log$fit
  pred_data_asexual$ci_lower <- 10^(pred_log$fit - 1.96 * pred_log$se.fit)
  pred_data_asexual$ci_upper <- 10^(pred_log$fit + 1.96 * pred_log$se.fit)

  # Storing results for plotting
  asexual_analysis <- list(
    model = asexual_model,
    data = asexual_data,
    predictions = pred_data_asexual,
    n_studies = nrow(asexual_data),
    total_participants = sum(asexual_data$study_n),
    fold_change = fold_change,
    p_value = asexual_summary$coefficients["k13_mutation_pct", "Pr(>|t|)"],
    r_squared = asexual_summary$r.squared
  )
} else {
  cat("Insufficient data for asexual density analysis (n =", nrow(asexual_data), "studies)\n")
  asexual_analysis <- list(error = "Insufficient data")
}


# Improving plots ---------------------------------------------------------


# Cleaner appearance
theme_malaria_clean <- theme_minimal() +
  theme(
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    legend.position = "bottom",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    legend.background = element_rect(fill = "white", color = NA)
  )


create_clean_resistance_plot <- function(data, y_var, analysis_result, x_label, y_label) {

  if("error" %in% names(analysis_result)) {
    return(ggplot() +
             labs(x = x_label, y = y_label) +
             theme_malaria_clean +
             annotate("text", x = 50, y = 50, label = "Insufficient Data",
                      size = 6, color = "grey50"))
  }

  plot_data <- data[!is.na(data[["k13_mutation_pct"]]) &
                      !is.na(data[[y_var]]) &
                      !is.na(data[["study_n"]]), ]

  # Calculating legend breaks
  weight_values <- plot_data$study_n[!is.na(plot_data$study_n)]

  if(length(weight_values) > 0) {
    breaks <- c(
      round(quantile(weight_values, 0.1, na.rm = TRUE), 0),
      round(quantile(weight_values, 0.5, na.rm = TRUE), 0),
      round(quantile(weight_values, 0.9, na.rm = TRUE), 0)
    )
    breaks <- unique(breaks)
  } else {
    breaks <- c(50, 100, 500)
  }

  p <- ggplot(plot_data, aes(x = k13_mutation_pct, y = .data[[y_var]])) +

    # Adding prediction ribbon from weighted binomial model
    geom_ribbon(data = analysis_result$predictions,
                aes(x = k13_resistance, ymin = ci_lower, ymax = ci_upper),
                alpha = 0.12, fill = "#2E86AB", inherit.aes = FALSE) +

    # Adding prediction line from weighted binomial model
    geom_line(data = analysis_result$predictions,
              aes(x = k13_resistance, y = predicted_prevalence),
              color = "#2E86AB", linetype = "solid", linewidth = 1.2, inherit.aes = FALSE) +

    # Adding simple linear trend line
    geom_smooth(method = "lm", se = TRUE, color = "#A23B72", fill = "#A23B72",
                alpha = 0.15, linetype = "dashed", linewidth = 0.8) +

    # Adding data points
    geom_point(aes(size = study_n), alpha = 0.8, color = "#F18F01", stroke = 0.5) +

    scale_size_continuous(
      name = "Sample Size",
      range = c(3, 15),
      breaks = breaks,
      labels = function(x) ifelse(x >= 1000, paste0(round(x/1000, 1), "k"), as.character(x))
    ) +

    scale_x_continuous(breaks = seq(0, 100, 20), limits = c(-2, 102)) +

    labs(x = x_label, y = y_label) +

    theme_malaria_clean

  # Adding statistical annotation in top corner
  multiplier_text <- if(y_var %in% c("baseline_gam_prevalence", "day7_gam_prevalence")) {
    paste0("Transmission multiplier: ", round(analysis_result$multiplier_best, 2), "×")
  } else ""

  or_text <- paste0("OR = ", round(analysis_result$odds_ratio, 3),
                    " (95% CI: ", round(analysis_result$or_ci_lower, 3),
                    "–", round(analysis_result$or_ci_upper, 3), ")")

  p_text <- paste0("p = ", format.pval(analysis_result$p_value))

  # Positioning annotation in top-left or top-right based on data distribution
  x_range <- range(plot_data$k13_mutation_pct, na.rm = TRUE)
  y_range <- range(plot_data[[y_var]], na.rm = TRUE)

  # Choosing position based on where there's more data points
  if(mean(plot_data$k13_mutation_pct, na.rm = TRUE) > 50) {
    x_pos <- x_range[1] + (x_range[2] - x_range[1]) * 0.05
    hjust_val <- 0
  } else {
    x_pos <- x_range[2] - (x_range[2] - x_range[1]) * 0.05
    hjust_val <- 1
  }

  y_pos <- y_range[2] - (y_range[2] - y_range[1]) * 0.05

  # Adding clean annotation box
  if(multiplier_text != "") {
    annotation_text <- paste(or_text, multiplier_text, p_text, sep = "\n")
  } else {
    annotation_text <- paste(or_text, p_text, sep = "\n")
  }

  p <- p +
    annotate("label", x = x_pos, y = y_pos,
             label = annotation_text,
             hjust = hjust_val, vjust = 1, size = 3.8, color = "black",
             fontface = "bold", fill = "white", alpha = 0.9,
             label.padding = unit(0.4, "lines"),
             label.r = unit(0.15, "lines"))

  return(p)
}

# Creating asexual density plot
create_clean_asexual_plot <- function(asexual_analysis) {

  if("error" %in% names(asexual_analysis)) {
    return(ggplot() +
             labs(x = "K13 Resistance (%)", y = "Baseline Asexual Density (parasites/μL)") +
             theme_malaria_clean +
             annotate("text", x = 50, y = 1000, label = "Insufficient Data",
                      size = 6, color = "grey50"))
  }

  # Calculating breaks for legend
  asexual_weights <- asexual_analysis$data$study_n[!is.na(asexual_analysis$data$study_n)]

  if(length(asexual_weights) > 0) {
    asexual_breaks <- c(
      round(quantile(asexual_weights, 0.1, na.rm = TRUE), 0),
      round(quantile(asexual_weights, 0.5, na.rm = TRUE), 0),
      round(quantile(asexual_weights, 0.9, na.rm = TRUE), 0)
    )
    asexual_breaks <- unique(asexual_breaks)
  } else {
    asexual_breaks <- c(50, 100, 500)
  }

  p <- ggplot(asexual_analysis$data, aes(x = k13_mutation_pct, y = baseline_asexual_density)) +

    # Adding prediction ribbon from weighted model
    geom_ribbon(data = asexual_analysis$predictions,
                aes(x = k13_mutation_pct, ymin = ci_lower, ymax = ci_upper),
                alpha = 0.12, fill = "#2E86AB", inherit.aes = FALSE) +

    # Adding prediction line from weighted model
    geom_line(data = asexual_analysis$predictions,
              aes(x = k13_mutation_pct, y = predicted_density),
              color = "#2E86AB", linetype = "solid", linewidth = 1.2, inherit.aes = FALSE) +

    # Adding simple linear trend line
    geom_smooth(method = "lm", se = TRUE, color = "#A23B72", fill = "#A23B72",
                alpha = 0.15, linetype = "dashed", linewidth = 0.8) +

    # Adding data points
    geom_point(aes(size = study_n), alpha = 0.8, color = "#F18F01", stroke = 0.5) +

    scale_y_log10(labels = function(x) ifelse(x >= 1000, paste0(x/1000, "k"), comma_format()(x))) +
    scale_x_continuous(breaks = seq(0, 100, 20), limits = c(-2, 102)) +

    scale_size_continuous(
      name = "Sample Size",
      range = c(3, 15),
      breaks = asexual_breaks,
      labels = function(x) ifelse(x >= 1000, paste0(round(x/1000, 1), "k"), as.character(x))
    ) +

    labs(
      x = "K13 Resistance (%)",
      y = "Baseline Asexual Density (parasites/μL, log scale)"
    ) +

    theme_malaria_clean

  # Adding statistical annotation
  annotation_text <- paste0("Fold change = ", round(asexual_analysis$fold_change, 2), "×\n",
                            "R² = ", round(asexual_analysis$r_squared, 3), "\n",
                            "p = ", format.pval(asexual_analysis$p_value))

  # Positioning in top-left
  x_pos <- 5
  y_pos <- max(asexual_analysis$data$baseline_asexual_density) * 0.8

  p <- p +
    annotate("label", x = x_pos, y = y_pos,
             label = annotation_text,
             hjust = 0, vjust = 1, size = 3.8, color = "black",
             fontface = "bold", fill = "white", alpha = 0.9,
             label.padding = unit(0.4, "lines"),
             label.r = unit(0.15, "lines"))

  return(p)
}

# Improved plots

if(!"error" %in% names(baseline_analysis)) {
  p_baseline_clean <- create_clean_resistance_plot(
    malaria_data, "baseline_gam_prevalence", baseline_analysis,
    "K13 Resistance (%)", "Baseline Gametocyte Prevalence (%)"
  )
  cat("Baseline gametocyte plot created\n")
}

if(!"error" %in% names(day7_analysis)) {
  p_day7_clean <- create_clean_resistance_plot(
    malaria_data, "day7_gam_prevalence", day7_analysis,
    "K13 Resistance (%)", "Day 7 Gametocyte Prevalence (%)"
  )
  cat("Day 7 gametocyte plot created\n")
}

if(!"error" %in% names(day3_analysis)) {
  p_day3_clean <- create_clean_resistance_plot(
    malaria_data, "day3_gam_prevalence", day3_analysis,
    "K13 Resistance (%)", "Day 3 Gametocyte Prevalence (%)"
  )
  cat("Day 3 gametocyte plot created\n")
}

if(!"error" %in% names(asexual_analysis)) {
  p_asexual_clean <- create_clean_asexual_plot(asexual_analysis)
  cat("Asexual density plot created\n")
}


if(exists("p_baseline_clean") && exists("p_day7_clean")) {
  # Main combined plot
  combined_clean <- p_baseline_clean / p_day7_clean
  print(combined_clean)

  # Saving combined plot
  ggsave("k13_gametocyte_analysis_clean.png", combined_clean,
         width = 12, height = 10, dpi = 300, bg = "white")
  cat("Main combined plot saved as 'k13_gametocyte_analysis_clean.png'\n")
}

if(exists("p_day3_clean")) {
  print(p_day3_clean)
  ggsave("k13_day3_analysis_clean.png", p_day3_clean,
         width = 12, height = 6, dpi = 300, bg = "white")
  cat("Day 3 plot saved as 'k13_day3_analysis_clean.png'\n")
}

if(exists("p_asexual_clean")) {
  print(p_asexual_clean)
  ggsave("k13_asexual_analysis_clean.png", p_asexual_clean,
         width = 12, height = 6, dpi = 300, bg = "white")
  cat("Asexual density plot saved as 'k13_asexual_analysis_clean.png'\n")
}

# 4-panel figure
if(exists("p_baseline_clean") && exists("p_day7_clean") &&
   exists("p_day3_clean") && exists("p_asexual_clean")) {

  comprehensive_plot <- (p_baseline_clean | p_day7_clean) /
    (p_day3_clean | p_asexual_clean)

  print(comprehensive_plot)
  ggsave("k13_comprehensive_analysis_clean.png", comprehensive_plot,
         width = 16, height = 12, dpi = 300, bg = "white")
  cat("4-panel plot saved as 'k13_comprehensive_analysis_clean.png'\n")
}


# Final summary -----------------------------------------------------------


if(!"error" %in% names(baseline_analysis) && !"error" %in% names(day7_analysis)) {
  cat("\nRECOMMENDED TRANSMISSION MULTIPLIERS:\n")
  cat("\n1. BASELINE GAMETOCYTE ADVANTAGE:\n")
  cat("   - Best estimate:", round(baseline_analysis$multiplier_best, 2), "x\n")
  cat("   - Range:", round(baseline_analysis$multiplier_conservative, 2), "-", round(baseline_analysis$multiplier_optimistic, 2), "x\n")
  cat("   - Evidence: Strong (", baseline_analysis$n_studies, " studies, ", format(baseline_analysis$total_participants, big.mark = ","), " participants)\n")

  cat("\n2. DAY 7 GAMETOCYTE ADVANTAGE:\n")
  cat("   - Best estimate:", round(day7_analysis$multiplier_best, 2), "x\n")
  cat("   - Range:", round(day7_analysis$multiplier_conservative, 2), "-", round(day7_analysis$multiplier_optimistic, 2), "x\n")
  cat("   - Evidence: Moderate (", day7_analysis$n_studies, " studies, ", format(day7_analysis$total_participants, big.mark = ","), " participants)\n")

  cat("\n3. ASEXUAL DENSITY RELATIONSHIP:\n")
  if(!"error" %in% names(asexual_analysis)) {
    cat("   - Fold change (0% to 100% K13):", round(asexual_analysis$fold_change, 2), "x\n")
    cat("   - R-squared:", round(asexual_analysis$r_squared, 3), "\n")
    cat("   - P-value:", format.pval(asexual_analysis$p_value), "\n")
    cat("   - Evidence: Moderate (", asexual_analysis$n_studies, " studies, ", format(asexual_analysis$total_participants, big.mark = ","), " participants)\n")
  } else {
    cat("   - Insufficient data\n")
  }

  cat("\nSENSITIVITY ANALYSIS RANGES:\n")
  cat("- Conservative scenario: 1.0x (no advantage)\n")
  cat("- Moderate scenario:", round(mean(c(baseline_analysis$multiplier_best, day7_analysis$multiplier_best)), 2), "x\n")
  cat("- Optimistic scenario:", round(max(baseline_analysis$multiplier_optimistic, day7_analysis$multiplier_optimistic), 2), "x\n")

  # Statistical significance check
  cat("\nSTATISTICAL SIGNIFICANCE:\n")
  cat("- Baseline gametocyte analysis p-value:", format.pval(baseline_analysis$p_value), "\n")
  cat("- Day 7 gametocyte analysis p-value:", format.pval(day7_analysis$p_value), "\n")
  if(!"error" %in% names(day3_analysis)) {
    cat("- Day 3 gametocyte analysis p-value:", format.pval(day3_analysis$p_value), "\n")
  }
  if(!"error" %in% names(asexual_analysis)) {
    cat("- Asexual density analysis p-value:", format.pval(asexual_analysis$p_value), "\n")
  }
}

