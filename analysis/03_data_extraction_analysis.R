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
    #geom_smooth(method = "lm", se = TRUE, color = "#A23B72", fill = "#A23B72",
                #alpha = 0.15, linetype = "dashed", linewidth = 0.8) +

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
    #geom_smooth(method = "lm", se = TRUE, color = "#A23B72", fill = "#A23B72",
                #alpha = 0.15, linetype = "dashed", linewidth = 0.8) +

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


# Introducing malaria prevalence as covariate  ----------------------------

# ========================
# K13 RESISTANCE ANALYSIS WITH PREVALENCE COVARIATE - CORRECTED VERSION
# ========================

# Required libraries
library(ggplot2)
library(dplyr)
library(broom)
library(weights)
library(patchwork)
library(readxl)
library(scales)
library(readr)

# Read in your malaria data from Excel file (or use data frame below)
malaria_data <- data.frame(
  Study = c(
    "Ladeia-Andrade 2016", "Olivera 2019", "Teklemariam 2015", "Tesfaye 2024",
    "Ippolito 2020", "Getnet 2013", "Atroosh 2014", "Chang 2016", "Chang 2016",
    "Djimde 2016", "Kyaw Myo Tun 2015", "Kyaw Myo Tun 2015", "Vantaux 2020",
    "Vantaux 2020", "Thriemer 2014", "Rovira-Vallbona 2020", "Rovira-Vallbona 2020",
    "Ashley 2014", "Ashley 2014", "Ashley 2014", "Ashley 2014", "Ashley 2014",
    "Ashley 2014", "Ashley 2014", "Carrara 2009", "Lek 2022"
  ),
  # Add country/location identifier for Ashley 2014 and other multi-location studies
  Location = c(
    "Brazil", "Colombia", "Ethiopia", "Ethiopia", "Zambia", "Ethiopia", "Yemen",
    "Uganda", "Uganda", "Mali", "Myanmar", "Myanmar", "Cambodia", "Cambodia",
    "Indonesia", "Cambodia", "Cambodia",
    # Ashley 2014 - 7 countries
    "Cambodia", "Thailand", "Vietnam", "Myanmar", "Bangladesh", "DRC", "Nigeria",
    "Thailand", "Cambodia"
  ),
  k13_mutation_pct = c(
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 46, 46, 63, 63, 80.7, 87, 87, 44.7, 73.6,
    27.5, 26, 0, 1.7, 0, 21.9, 57.9
  ),
  baseline_gam_prevalence = c(
    32.5, 12.5, 7.6, 6.3, 13, 10, 40.7, 37, 50, 22.8, 22, 21, 49, 37.5, 18.9,
    NA, 76, 11.8, 6.7, 5.9, 11.3, 0, 30.8, 2.9, 6.3, 46.9
  ),
  baseline_asexual_density = c(
    4700.1, 3527, 27798, 10627, 13000, NA, 8199, NA, NA, NA, 3084, 2140, NA, NA,
    8233, 10166, 15832, 55623, 36529, 49738, 66066, 27130, 60037, 52250, 6982, 11873
  ),
  baseline_gam_density = c(
    NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
    NA, NA, NA, NA, NA, NA, NA
  ),
  day3_gam_prevalence = c(
    23.6, 18.6, 1, 2.5, 3, 37.5, NA, NA, NA, 17.4, 8, 13, NA, NA, NA, NA, NA,
    NA, NA, NA, NA, NA, NA, NA, NA, NA
  ),
  day7_gam_prevalence = c(
    NA, 14.3, NA, 2.5, 2, 0, 34.9, 7.84, 0, 23.9, 0, 1.4, NA, NA, 20.5, NA,
    NA, NA, NA, NA, NA, NA, NA, NA, NA, NA
  ),
  sample_size = c(
    162, 84, 92, 73, 94, 80, 89, 61, 9, 92, 78, 76, 53, 56, 95, 33, 27, 460,
    185, 120, 80, 56, 120, 40, 3264, 211
  ),
  stringsAsFactors = FALSE
)

# Clean column names
colnames(malaria_data) <- trimws(colnames(malaria_data))
malaria_data <- malaria_data %>%
  rename(study_n = sample_size)

# Convert "NA" strings to actual NA values and ensure numeric columns
malaria_data[malaria_data == "NA"] <- NA
numeric_cols <- c("k13_mutation_pct", "baseline_gam_prevalence", "baseline_asexual_density",
                  "baseline_gam_density", "day3_gam_prevalence", "day7_gam_prevalence", "study_n")

for(col in numeric_cols) {
  if(col %in% colnames(malaria_data)) {
    malaria_data[[col]] <- as.numeric(malaria_data[[col]])
  }
}

# Read the study locations data
study_locations <- read_csv("analysis/data-derived/study_locations_processed_complete.csv")

# Check the structure of study locations
cat("Study locations structure:\n")
print(study_locations %>%
        group_by(study) %>%
        summarise(n_locations = n(), .groups = 'drop') %>%
        arrange(desc(n_locations)))

# Create the proper mapping based on your extraction description
# Ashley 2014 should have 7 country entries (Cambodia and Thailand averaged from multiple sites)
# Chang 2016 should be excluded (no study dates)

# Create study-location mapping for prevalence
study_location_mapping <- study_locations %>%
  # Create a study-location key that matches your malaria data structure
  mutate(
    study_location_key = case_when(
      # For Ashley 2014, we need to map to countries
      study == "Ashley 2014" & grepl("Cambodia", location) ~ "Ashley 2014_Cambodia",
      study == "Ashley 2014" & grepl("Thailand", location) ~ "Ashley 2014_Thailand",
      study == "Ashley 2014" & grepl("Vietnam", location) ~ "Ashley 2014_Vietnam",
      study == "Ashley 2014" & grepl("Myanmar", location) ~ "Ashley 2014_Myanmar",
      study == "Ashley 2014" & grepl("Bangladesh", location) ~ "Ashley 2014_Bangladesh",
      study == "Ashley 2014" & grepl("DRC", location) ~ "Ashley 2014_DRC",
      study == "Ashley 2014" & grepl("Nigeria", location) ~ "Ashley 2014_Nigeria",
      # For other multi-location studies, take first entry or average
      TRUE ~ paste0(study, "_", "main")
    )
  ) %>%
  # Average prevalence by study-location key (this handles Cambodia/Thailand averaging)
  group_by(study_location_key) %>%
  summarise(
    study = first(study),
    malaria_prevalence = mean(prev, na.rm = TRUE),
    n_sites = n(),
    .groups = 'drop'
  ) %>%
  # Extract location from key for matching
  mutate(
    location_match = case_when(
      grepl("_Cambodia$", study_location_key) ~ "Cambodia",
      grepl("_Thailand$", study_location_key) ~ "Thailand",
      grepl("_Vietnam$", study_location_key) ~ "Vietnam",
      grepl("_Myanmar$", study_location_key) ~ "Myanmar",
      grepl("_Bangladesh$", study_location_key) ~ "Bangladesh",
      grepl("_DRC$", study_location_key) ~ "DRC",
      grepl("_Nigeria$", study_location_key) ~ "Nigeria",
      TRUE ~ "main"
    )
  )

cat("\nProcessed study-location mapping:\n")
print(study_location_mapping)

# Now merge with malaria data
# Create matching key in malaria data
malaria_data_with_key <- malaria_data %>%
  mutate(
    location_match = case_when(
      Study == "Ashley 2014" ~ Location,  # Use the Location column for Ashley
      TRUE ~ "main"  # For other studies, use main
    ),
    study_location_key = paste0(Study, "_", location_match)
  )

# Merge with prevalence data
malaria_data_enhanced <- malaria_data_with_key %>%
  left_join(
    study_location_mapping %>% select(study_location_key, malaria_prevalence, n_sites),
    by = "study_location_key"
  ) %>%
  # Convert prevalence from proportion to percentage
  mutate(malaria_prevalence_pct = malaria_prevalence * 100) %>%
  # Exclude Chang 2016 from prevalence analysis (mark as excluded)
  mutate(
    exclude_from_prevalence = Study == "Chang 2016",
    malaria_prevalence_pct = ifelse(exclude_from_prevalence, NA, malaria_prevalence_pct)
  )

# Check which studies have prevalence data
cat("\nStudies with prevalence data:\n")
prevalence_studies <- malaria_data_enhanced %>%
  filter(!is.na(malaria_prevalence_pct)) %>%
  select(Study, Location, malaria_prevalence_pct, n_sites) %>%
  distinct()
print(prevalence_studies)

cat("\nStudies WITHOUT prevalence data (excluded from prevalence analysis):\n")
excluded_studies <- malaria_data_enhanced %>%
  filter(is.na(malaria_prevalence_pct)) %>%
  select(Study, Location) %>%
  distinct()
print(excluded_studies)

cat("\nSample size check - Ashley 2014 countries:\n")
ashley_check <- malaria_data_enhanced %>%
  filter(Study == "Ashley 2014") %>%
  select(Study, Location, k13_mutation_pct, malaria_prevalence_pct, n_sites)
print(ashley_check)

# Enhanced weighted binomial regression with prevalence covariate
perform_weighted_binomial_with_prevalence <- function(data, outcome_var,
                                                      predictor_var = "k13_mutation_pct",
                                                      prevalence_var = "malaria_prevalence_pct",
                                                      weight_var = "study_n") {

  # Filter complete cases including prevalence (excludes Chang 2016 automatically)
  complete_data <- data[!is.na(data[[outcome_var]]) &
                          !is.na(data[[predictor_var]]) &
                          !is.na(data[[prevalence_var]]) &
                          !is.na(data[[weight_var]]), ]

  if(nrow(complete_data) < 4) {  # Need more data for multiple predictors
    return(list(error = "Insufficient data"))
  }

  cat("\nStudies included in", outcome_var, "analysis with prevalence covariate:\n")
  for(i in 1:nrow(complete_data)) {
    cat(sprintf("  %s (%s): K13=%.1f%%, %s=%.1f%%, Malaria prev=%.1f%%, n=%d\n",
                complete_data$Study[i],
                complete_data$Location[i],
                complete_data[[predictor_var]][i],
                outcome_var,
                complete_data[[outcome_var]][i],
                complete_data[[prevalence_var]][i],
                complete_data[[weight_var]][i]))
  }

  # Create binomial response
  complete_data$successes <- round((complete_data[[outcome_var]]/100) * complete_data[[weight_var]])
  complete_data$failures <- complete_data[[weight_var]] - complete_data$successes

  # Ensure valid counts
  complete_data$successes <- pmax(0, pmin(complete_data$successes, complete_data[[weight_var]]))
  complete_data$failures <- complete_data[[weight_var]] - complete_data$successes

  # Fit models with and without prevalence
  formula_base <- paste("cbind(successes, failures) ~", predictor_var)
  formula_full <- paste("cbind(successes, failures) ~", predictor_var, "+", prevalence_var)

  model_base <- glm(as.formula(formula_base), family = binomial, data = complete_data)
  model_full <- glm(as.formula(formula_full), family = binomial, data = complete_data)

  # Compare models using AIC and likelihood ratio test
  aic_comparison <- data.frame(
    Model = c("K13 only", "K13 + Prevalence"),
    AIC = c(AIC(model_base), AIC(model_full)),
    df = c(model_base$df.residual, model_full$df.residual)
  )

  # Likelihood ratio test
  lrt <- anova(model_base, model_full, test = "Chisq")

  # Extract coefficients from full model
  summary_full <- summary(model_full)
  coefficients <- summary_full$coefficients

  # K13 effects
  k13_beta <- coefficients[predictor_var, "Estimate"]
  k13_se <- coefficients[predictor_var, "Std. Error"]
  k13_or <- exp(k13_beta)
  k13_or_ci_lower <- exp(k13_beta - 1.96 * k13_se)
  k13_or_ci_upper <- exp(k13_beta + 1.96 * k13_se)
  k13_p_value <- coefficients[predictor_var, "Pr(>|z|)"]

  # Prevalence effects
  prev_beta <- coefficients[prevalence_var, "Estimate"]
  prev_se <- coefficients[prevalence_var, "Std. Error"]
  prev_or <- exp(prev_beta)
  prev_or_ci_lower <- exp(prev_beta - 1.96 * prev_se)
  prev_or_ci_upper <- exp(prev_beta + 1.96 * prev_se)
  prev_p_value <- coefficients[prevalence_var, "Pr(>|z|)"]

  return(list(
    model_base = model_base,
    model_full = model_full,
    data = complete_data,
    n_studies = nrow(complete_data),
    total_participants = sum(complete_data[[weight_var]]),

    # K13 results
    k13_odds_ratio = k13_or,
    k13_or_ci_lower = k13_or_ci_lower,
    k13_or_ci_upper = k13_or_ci_upper,
    k13_p_value = k13_p_value,

    # Prevalence results
    prevalence_odds_ratio = prev_or,
    prevalence_or_ci_lower = prev_or_ci_lower,
    prevalence_or_ci_upper = prev_or_ci_upper,
    prevalence_p_value = prev_p_value,

    # Model comparison
    aic_comparison = aic_comparison,
    lrt_p_value = if(length(lrt$`Pr(>Chi)`) >= 2) lrt$`Pr(>Chi)`[2] else NA,

    aic_base = AIC(model_base),
    aic_full = AIC(model_full)
  ))
}

# Enhanced plotting function
create_enhanced_resistance_plot <- function(data, y_var, analysis_result,
                                            x_label, y_label, prevalence_var = "malaria_prevalence_pct") {

  if("error" %in% names(analysis_result)) {
    return(ggplot() +
             labs(x = x_label, y = y_label) +
             theme_minimal() +
             annotate("text", x = 50, y = 50, label = "Insufficient Data",
                      size = 6, color = "grey50"))
  }

  plot_data <- analysis_result$data

  # Get unique studies and create enough shapes
  unique_studies <- unique(plot_data$Study)
  n_studies <- length(unique_studies)

  # Create a comprehensive set of shapes
  all_shapes <- c(16, 17, 18, 15, 3, 4, 8, 11, 12, 13, 14, 19, 20, 21, 22, 23, 24, 25, 0, 1, 2, 5, 6, 7, 9, 10)
  study_shapes <- all_shapes[1:n_studies]
  names(study_shapes) <- unique_studies

  # Generate predictions from the full model for smooth regression line
  k13_range <- seq(0, 100, by = 1)
  mean_prevalence <- mean(plot_data[[prevalence_var]], na.rm = TRUE)

  # Create prediction data
  pred_data <- data.frame(
    k13_mutation_pct = k13_range,
    malaria_prevalence_pct = mean_prevalence
  )
  names(pred_data) <- c("k13_mutation_pct", prevalence_var)

  # Generate predictions from the full model
  pred_logit <- predict(analysis_result$model_full, newdata = pred_data, se.fit = TRUE)

  # Convert to response scale (percentages)
  pred_df <- data.frame(
    k13_mutation_pct = k13_range,
    predicted = plogis(pred_logit$fit) * 100,
    se = pred_logit$se.fit
  )
  pred_df$ci_lower <- plogis(pred_logit$fit - 1.96 * pred_df$se) * 100
  pred_df$ci_upper <- plogis(pred_logit$fit + 1.96 * pred_df$se) * 100

  # Create base plot
  p <- ggplot(plot_data, aes(x = k13_mutation_pct, y = .data[[y_var]])) +

    # Add prediction ribbon and line FIRST (so points appear on top)
    geom_ribbon(data = pred_df,
                aes(x = k13_mutation_pct, ymin = ci_lower, ymax = ci_upper),
                alpha = 0.2, fill = "red", inherit.aes = FALSE) +

    geom_line(data = pred_df,
              aes(x = k13_mutation_pct, y = predicted),
              color = "red", linewidth = 1.2, inherit.aes = FALSE) +

    # Color points by malaria prevalence, shape by study
    geom_point(aes(size = study_n, color = .data[[prevalence_var]], shape = Study),
               alpha = 0.8, stroke = 0.5) +

    scale_color_viridis_c(name = "Malaria\nPrevalence (%)", option = "plasma") +
    scale_size_continuous(name = "Sample Size", range = c(3, 15)) +
    scale_shape_manual(name = "Study", values = study_shapes) +

    scale_x_continuous(breaks = seq(0, 100, 20), limits = c(-2, 102)) +
    scale_y_continuous(breaks = seq(0, 80, 10), limits = c(0, max(plot_data[[y_var]], na.rm = TRUE) * 1.1)) +

    labs(x = x_label, y = y_label,
         title = paste(y_label, "with Malaria Prevalence Covariate"),
         subtitle = paste("Red line: Predicted relationship (prevalence =", round(mean_prevalence, 1), "%)")) +

    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      legend.position = "right",
      legend.box = "vertical",
      panel.grid.minor = element_blank()
    ) +

    # Adjust legend guides to prevent overcrowding
    guides(
      color = guide_colorbar(barwidth = 1, barheight = 6),
      size = guide_legend(override.aes = list(shape = 16)),
      shape = guide_legend(override.aes = list(size = 4), ncol = 1)
    )

  # Add model comparison annotation
  if(!is.null(analysis_result$aic_comparison)) {
    model_text <- paste0(
      "Model Results:\n",
      "K13 only AIC: ", round(analysis_result$aic_base, 1), "\n",
      "K13 + Prevalence AIC: ", round(analysis_result$aic_full, 1), "\n",
      "AIC improvement: ", round(analysis_result$aic_base - analysis_result$aic_full, 1), "\n",
      "K13 OR: ", round(analysis_result$k13_odds_ratio, 3),
      " (p=", format.pval(analysis_result$k13_p_value), ")\n",
      "Prevalence OR: ", round(analysis_result$prevalence_odds_ratio, 3),
      " (p=", format.pval(analysis_result$prevalence_p_value), ")"
    )

    p <- p + annotate("label", x = 5, y = max(plot_data[[y_var]], na.rm = TRUE) * 0.95,
                      label = model_text,
                      hjust = 0, vjust = 1, size = 3.2,
                      fill = "white", alpha = 0.9,
                      label.padding = unit(0.5, "lines"))
  }

  return(p)
}

# Create a simpler plot function for just the regression relationship
create_simple_regression_plot <- function(data, y_var, analysis_result, x_label, y_label) {

  if("error" %in% names(analysis_result)) {
    return(ggplot() +
             labs(x = x_label, y = y_label) +
             theme_minimal() +
             annotate("text", x = 50, y = 50, label = "Insufficient Data", size = 6))
  }

  plot_data <- analysis_result$data

  # Generate predictions for smooth line
  k13_range <- seq(0, 100, by = 1)
  mean_prevalence <- mean(plot_data$malaria_prevalence_pct, na.rm = TRUE)

  pred_data <- data.frame(
    k13_mutation_pct = k13_range,
    malaria_prevalence_pct = mean_prevalence
  )

  pred_logit <- predict(analysis_result$model_full, newdata = pred_data, se.fit = TRUE)

  pred_df <- data.frame(
    k13_mutation_pct = k13_range,
    predicted = plogis(pred_logit$fit) * 100,
    ci_lower = plogis(pred_logit$fit - 1.96 * pred_logit$se.fit) * 100,
    ci_upper = plogis(pred_logit$fit + 1.96 * pred_logit$se.fit) * 100
  )

  p <- ggplot(plot_data, aes(x = k13_mutation_pct, y = .data[[y_var]])) +

    # Prediction ribbon and line
    geom_ribbon(data = pred_df,
                aes(x = k13_mutation_pct, ymin = ci_lower, ymax = ci_upper),
                alpha = 0.15, fill = "#2E86AB", inherit.aes = FALSE) +

    geom_line(data = pred_df,
              aes(x = k13_mutation_pct, y = predicted),
              color = "#2E86AB", linewidth = 1.5, inherit.aes = FALSE) +

    # Data points
    geom_point(aes(size = study_n), alpha = 0.7, color = "#F18F01") +

    scale_size_continuous(name = "Sample Size", range = c(3, 12)) +
    scale_x_continuous(breaks = seq(0, 100, 20)) +

    labs(x = x_label, y = y_label,
         title = paste("K13 Resistance vs", gsub(" \\(%\\)", "", y_label)),
         subtitle = "Controlling for malaria prevalence") +

    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      legend.position = "bottom"
    )

  # Add statistics annotation
  stats_text <- paste0(
    "Adjusted OR = ", round(analysis_result$k13_odds_ratio, 3), "\n",
    "95% CI: ", round(analysis_result$k13_or_ci_lower, 3), "-", round(analysis_result$k13_or_ci_upper, 3), "\n",
    "p = ", format.pval(analysis_result$k13_p_value)
  )

  p <- p + annotate("label", x = 75, y = max(plot_data[[y_var]], na.rm = TRUE) * 0.9,
                    label = stats_text, hjust = 0, vjust = 1, size = 3.5,
                    fill = "white", alpha = 0.9)

  return(p)
}

# Run enhanced analyses
cat("\n=== BASELINE GAMETOCYTE ANALYSIS WITH PREVALENCE ===\n")
baseline_analysis_enhanced <- perform_weighted_binomial_with_prevalence(
  malaria_data_enhanced, "baseline_gam_prevalence"
)

cat("\n=== DAY 3 GAMETOCYTE ANALYSIS WITH PREVALENCE ===\n")
day3_analysis_enhanced <- perform_weighted_binomial_with_prevalence(
  malaria_data_enhanced, "day3_gam_prevalence"
)

cat("\n=== DAY 7 GAMETOCYTE ANALYSIS WITH PREVALENCE ===\n")
day7_analysis_enhanced <- perform_weighted_binomial_with_prevalence(
  malaria_data_enhanced, "day7_gam_prevalence"
)

# Create enhanced plots
if(!"error" %in% names(baseline_analysis_enhanced)) {
  # Enhanced plot with all details
  p_baseline_enhanced <- create_enhanced_resistance_plot(
    malaria_data_enhanced, "baseline_gam_prevalence", baseline_analysis_enhanced,
    "K13 Resistance (%)", "Baseline Gametocyte Prevalence (%)"
  )
  print(p_baseline_enhanced)
  ggsave("baseline_gametocyte_enhanced_corrected.png", p_baseline_enhanced,
         width = 16, height = 9, dpi = 300)

  # Simple regression plot
  p_baseline_simple <- create_simple_regression_plot(
    malaria_data_enhanced, "baseline_gam_prevalence", baseline_analysis_enhanced,
    "K13 Resistance (%)", "Baseline Gametocyte Prevalence (%)"
  )
  print(p_baseline_simple)
  ggsave("baseline_gametocyte_regression.png", p_baseline_simple,
         width = 10, height = 7, dpi = 300)
}

if(!"error" %in% names(day3_analysis_enhanced)) {
  # Enhanced plot with all details
  p_day3_enhanced <- create_enhanced_resistance_plot(
    malaria_data_enhanced, "day3_gam_prevalence", day3_analysis_enhanced,
    "K13 Resistance (%)", "Day 3 Gametocyte Prevalence (%)"
  )
  print(p_day3_enhanced)
  ggsave("day3_gametocyte_enhanced_corrected.png", p_day3_enhanced,
         width = 16, height = 9, dpi = 300)

  # Simple regression plot
  p_day3_simple <- create_simple_regression_plot(
    malaria_data_enhanced, "day3_gam_prevalence", day3_analysis_enhanced,
    "K13 Resistance (%)", "Day 3 Gametocyte Prevalence (%)"
  )
  print(p_day3_simple)
  ggsave("day3_gametocyte_regression.png", p_day3_simple,
         width = 10, height = 7, dpi = 300)
}

if(!"error" %in% names(day7_analysis_enhanced)) {
  # Enhanced plot with all details
  p_day7_enhanced <- create_enhanced_resistance_plot(
    malaria_data_enhanced, "day7_gam_prevalence", day7_analysis_enhanced,
    "K13 Resistance (%)", "Day 7 Gametocyte Prevalence (%)"
  )
  print(p_day7_enhanced)
  ggsave("day7_gametocyte_enhanced_corrected.png", p_day7_enhanced,
         width = 16, height = 9, dpi = 300)

  # Simple regression plot
  p_day7_simple <- create_simple_regression_plot(
    malaria_data_enhanced, "day7_gam_prevalence", day7_analysis_enhanced,
    "K13 Resistance (%)", "Day 7 Gametocyte Prevalence (%)"
  )
  print(p_day7_simple)
  ggsave("day7_gametocyte_regression.png", p_day7_simple,
         width = 10, height = 7, dpi = 300)
}

# Print enhanced results function
print_enhanced_results <- function(analysis_result, outcome_name) {
  if("error" %in% names(analysis_result)) {
    cat("\n", outcome_name, "- Insufficient data\n")
    return()
  }

  cat("\n=== ENHANCED RESULTS FOR", toupper(outcome_name), "===\n")
  cat("Studies included:", analysis_result$n_studies, "\n")
  cat("Total participants:", format(analysis_result$total_participants, big.mark = ","), "\n\n")

  cat("K13 RESISTANCE EFFECTS (adjusted for malaria prevalence):\n")
  cat("Odds ratio:", round(analysis_result$k13_odds_ratio, 4), "\n")
  cat("95% CI:", round(analysis_result$k13_or_ci_lower, 4), "to",
      round(analysis_result$k13_or_ci_upper, 4), "\n")
  cat("P-value:", format.pval(analysis_result$k13_p_value), "\n\n")

  cat("MALARIA PREVALENCE EFFECTS:\n")
  cat("Odds ratio per 1% increase in prevalence:", round(analysis_result$prevalence_odds_ratio, 4), "\n")
  cat("95% CI:", round(analysis_result$prevalence_or_ci_lower, 4), "to",
      round(analysis_result$prevalence_or_ci_upper, 4), "\n")
  cat("P-value:", format.pval(analysis_result$prevalence_p_value), "\n\n")

  cat("MODEL COMPARISON:\n")
  print(analysis_result$aic_comparison)
  if(!is.na(analysis_result$lrt_p_value)) {
    cat("Likelihood ratio test p-value:", format.pval(analysis_result$lrt_p_value), "\n")
  }
  aic_improvement <- analysis_result$aic_base - analysis_result$aic_full
  cat("AIC improvement:", round(aic_improvement, 2),
      if(aic_improvement > 0) "(better fit)" else "(worse fit)", "\n")
}

# Print results
print_enhanced_results(baseline_analysis_enhanced, "Baseline Gametocytes")
print_enhanced_results(day3_analysis_enhanced, "Day 3 Gametocytes")
print_enhanced_results(day7_analysis_enhanced, "Day 7 Gametocytes")

# Summary of findings
cat("\n=== SUMMARY OF FINDINGS ===\n")
cat("This analysis properly handles:\n")
cat("- Ashley 2014: 7 countries with Cambodia/Thailand averaged from multiple sites\n")
cat("- Chang 2016: Excluded from prevalence analysis (no study dates)\n")
cat("- Other studies: Single locations or appropriately handled\n\n")

if(!"error" %in% names(baseline_analysis_enhanced)) {
  aic_change_baseline <- baseline_analysis_enhanced$aic_base - baseline_analysis_enhanced$aic_full
  cat("BASELINE GAMETOCYTES:\n")
  cat("- Including malaria prevalence",
      if(aic_change_baseline > 2) "substantially improves" else if(aic_change_baseline > 0) "improves" else "does not improve",
      "model fit\n")
  cat("- AIC change:", round(aic_change_baseline, 2), "\n")
  cat("- K13 effect (adjusted):", round(baseline_analysis_enhanced$k13_odds_ratio, 3), "\n")
}

if(!"error" %in% names(day3_analysis_enhanced)) {
  aic_change_day3 <- day3_analysis_enhanced$aic_base - day3_analysis_enhanced$aic_full
  cat("\nDAY 3 GAMETOCYTES:\n")
  cat("- Including malaria prevalence",
      if(aic_change_day3 > 2) "substantially improves" else if(aic_change_day3 > 0) "improves" else "does not improve",
      "model fit\n")
  cat("- AIC change:", round(aic_change_day3, 2), "\n")
  cat("- K13 effect (adjusted):", round(day3_analysis_enhanced$k13_odds_ratio, 3), "\n")
  cat("- INTERPRETATION: K13 resistance",
      if(day3_analysis_enhanced$k13_odds_ratio < 1) "DECREASES" else "INCREASES",
      "Day 3 gametocytes\n")
}

if(!"error" %in% names(day7_analysis_enhanced)) {
  aic_change_day7 <- day7_analysis_enhanced$aic_base - day7_analysis_enhanced$aic_full
  cat("\nDAY 7 GAMETOCYTES:\n")
  cat("- Including malaria prevalence",
      if(aic_change_day7 > 2) "substantially improves" else if(aic_change_day7 > 0) "improves" else "does not improve",
      "model fit\n")
  cat("- AIC change:", round(aic_change_day7, 2), "\n")
  cat("- K13 effect (adjusted):", round(day7_analysis_enhanced$k13_odds_ratio, 3), "\n")
}


# Transmission Multipliers ------------------------------------------------


# Function to calculate transmission multipliers from analysis results
calculate_transmission_multipliers <- function(analysis_result, outcome_name, prevalence_var = "malaria_prevalence_pct") {

  if("error" %in% names(analysis_result)) {
    cat("\n", outcome_name, "- Insufficient data for multiplier calculation\n")
    return(NULL)
  }

  cat("\n=== TRANSMISSION MULTIPLIERS FOR", toupper(outcome_name), "===\n")

  # Use mean malaria prevalence from the data for predictions
  mean_prevalence <- mean(analysis_result$data[[prevalence_var]], na.rm = TRUE)
  cat("Using mean malaria prevalence:", round(mean_prevalence, 2), "%\n")

  # Create prediction data for 0% and 100% K13 resistance
  pred_data_0 <- data.frame(
    k13_mutation_pct = 0,
    malaria_prevalence_pct = mean_prevalence
  )
  names(pred_data_0) <- c("k13_mutation_pct", prevalence_var)

  pred_data_100 <- data.frame(
    k13_mutation_pct = 100,
    malaria_prevalence_pct = mean_prevalence
  )
  names(pred_data_100) <- c("k13_mutation_pct", prevalence_var)

  # Get predictions from the FULL model (adjusted for malaria prevalence)
  pred_0_logit <- predict(analysis_result$model_full, newdata = pred_data_0, se.fit = TRUE)
  pred_100_logit <- predict(analysis_result$model_full, newdata = pred_data_100, se.fit = TRUE)

  # Convert to prevalence scale (percentages)
  prev_at_0 <- plogis(pred_0_logit$fit) * 100
  prev_at_100 <- plogis(pred_100_logit$fit) * 100

  # Calculate confidence intervals
  prev_0_ci_lower <- plogis(pred_0_logit$fit - 1.96 * pred_0_logit$se.fit) * 100
  prev_0_ci_upper <- plogis(pred_0_logit$fit + 1.96 * pred_0_logit$se.fit) * 100

  prev_100_ci_lower <- plogis(pred_100_logit$fit - 1.96 * pred_100_logit$se.fit) * 100
  prev_100_ci_upper <- plogis(pred_100_logit$fit + 1.96 * pred_100_logit$se.fit) * 100

  # Calculate transmission multipliers
  multiplier_best <- prev_at_100 / prev_at_0
  multiplier_conservative <- prev_100_ci_lower / prev_0_ci_upper  # Most conservative
  multiplier_optimistic <- prev_100_ci_upper / prev_0_ci_lower   # Most optimistic

  # Print results
  cat("\nPREDICTED PREVALENCES (adjusted for malaria prevalence):\n")
  cat("At 0% K13 resistance:", round(prev_at_0, 2), "% (95% CI:",
      round(prev_0_ci_lower, 2), "-", round(prev_0_ci_upper, 2), "%)\n")
  cat("At 100% K13 resistance:", round(prev_at_100, 2), "% (95% CI:",
      round(prev_100_ci_lower, 2), "-", round(prev_100_ci_upper, 2), "%)\n")

  cat("\nTRANSMISSION MULTIPLIERS:\n")
  cat("Best estimate:", round(multiplier_best, 3), "x\n")
  cat("Conservative (lower CI):", round(multiplier_conservative, 3), "x\n")
  cat("Optimistic (upper CI):", round(multiplier_optimistic, 3), "x\n")

  # Statistical significance interpretation
  k13_p <- analysis_result$k13_p_value
  cat("\nSTATISTICAL INTERPRETATION:\n")
  cat("K13 effect p-value:", format.pval(k13_p), "\n")
  if(k13_p > 0.05) {
    cat("=> K13 resistance has NO SIGNIFICANT effect on", outcome_name, "\n")
    cat("=> Recommended multiplier: 1.0x (no transmission advantage)\n")
  } else {
    cat("=> K13 resistance has SIGNIFICANT effect on", outcome_name, "\n")
    cat("=> Use calculated multipliers above\n")
  }

  # Return results for further use
  return(list(
    outcome = outcome_name,
    prev_at_0 = prev_at_0,
    prev_at_100 = prev_at_100,
    multiplier_best = multiplier_best,
    multiplier_conservative = multiplier_conservative,
    multiplier_optimistic = multiplier_optimistic,
    k13_p_value = k13_p,
    significant = k13_p <= 0.05,
    mean_prevalence_used = mean_prevalence
  ))
}

# Calculate multipliers for all outcomes
baseline_multipliers <- calculate_transmission_multipliers(
  baseline_analysis_enhanced, "Baseline Gametocytes"
)

day3_multipliers <- calculate_transmission_multipliers(
  day3_analysis_enhanced, "Day 3 Gametocytes"
)

day7_multipliers <- calculate_transmission_multipliers(
  day7_analysis_enhanced, "Day 7 Gametocytes"
)

# Create summary table of all multipliers
if(!is.null(baseline_multipliers) || !is.null(day3_multipliers) || !is.null(day7_multipliers)) {

  cat("\n" , rep("=", 60), "\n")
  cat("SUMMARY: REVISED TRANSMISSION MULTIPLIERS\n")
  cat(rep("=", 60), "\n")

  # Collect all valid results
  all_multipliers <- list()
  if(!is.null(baseline_multipliers)) all_multipliers[["Baseline"]] <- baseline_multipliers
  if(!is.null(day3_multipliers)) all_multipliers[["Day 3"]] <- day3_multipliers
  if(!is.null(day7_multipliers)) all_multipliers[["Day 7"]] <- day7_multipliers

  if(length(all_multipliers) > 0) {
    # Create summary data frame
    summary_df <- data.frame(
      Outcome = names(all_multipliers),
      Best_Estimate = sapply(all_multipliers, function(x) round(x$multiplier_best, 3)),
      Conservative = sapply(all_multipliers, function(x) round(x$multiplier_conservative, 3)),
      Optimistic = sapply(all_multipliers, function(x) round(x$multiplier_optimistic, 3)),
      P_Value = sapply(all_multipliers, function(x) format.pval(x$k13_p_value)),
      Significant = sapply(all_multipliers, function(x) x$significant),
      stringsAsFactors = FALSE
    )

    cat("\nSUMMARY TABLE:\n")
    print(summary_df, row.names = FALSE)

    # Overall recommendation
    cat("\nOVERALL RECOMMENDATIONS:\n")

    significant_outcomes <- sum(summary_df$Significant)
    if(significant_outcomes == 0) {
      cat("=> NO significant K13 effects found after controlling for malaria prevalence\n")
      cat("=> RECOMMENDED TRANSMISSION MULTIPLIER: 1.0x (no advantage)\n")
      cat("=> RANGE FOR SENSITIVITY ANALYSIS: 0.95x - 1.05x\n")
    } else {
      # Use the significant outcomes to make recommendations
      sig_multipliers <- summary_df[summary_df$Significant, ]
      if(nrow(sig_multipliers) > 0) {
        mean_best <- mean(sig_multipliers$Best_Estimate)
        mean_conservative <- mean(sig_multipliers$Conservative)
        mean_optimistic <- mean(sig_multipliers$Optimistic)

        cat("=>", significant_outcomes, "of", nrow(summary_df), "outcomes show significant K13 effects\n")
        cat("=> RECOMMENDED TRANSMISSION MULTIPLIER:", round(mean_best, 2), "x\n")
        cat("=> RANGE FOR SENSITIVITY ANALYSIS:", round(mean_conservative, 2), "x -", round(mean_optimistic, 2), "x\n")
      }
    }

    # Comparison with original analysis
    cat("\nCOMPARISON WITH ORIGINAL ANALYSIS:\n")
    cat("Original analysis (without prevalence covariate):\n")
    cat("- Baseline: ~2.9x transmission advantage\n")
    cat("- Day 7: ~0.9x transmission advantage\n")
    cat("\nRevised analysis (with prevalence covariate):\n")
    cat("- Baseline:", round(baseline_multipliers$multiplier_best, 2), "x (p =", format.pval(baseline_multipliers$k13_p_value), ")\n")
    if(!is.null(day7_multipliers)) {
      cat("- Day 7:", round(day7_multipliers$multiplier_best, 2), "x (p =", format.pval(day7_multipliers$k13_p_value), ")\n")
    }
    cat("\n=> Controlling for malaria prevalence ELIMINATES the apparent transmission advantage\n")
    cat("=> Original effects were likely due to CONFOUNDING by transmission intensity\n")
  }
}

# Function to create a visual comparison plot
create_multiplier_comparison_plot <- function() {

  if(exists("baseline_multipliers") && exists("day7_multipliers")) {

    # Data for comparison plot
    comparison_data <- data.frame(
      Analysis = rep(c("Original\n(unadjusted)", "Revised\n(adjusted)"), each = 2),
      Outcome = rep(c("Baseline", "Day 7"), 2),
      Multiplier = c(
        2.9, 0.9,  # Original estimates (approximate from your earlier results)
        baseline_multipliers$multiplier_best,
        ifelse(!is.null(day7_multipliers), day7_multipliers$multiplier_best, 1.0)
      ),
      Significant = c(TRUE, TRUE,
                      baseline_multipliers$significant,
                      ifelse(!is.null(day7_multipliers), day7_multipliers$significant, FALSE))
    )

    library(ggplot2)

    p <- ggplot(comparison_data, aes(x = Analysis, y = Multiplier, fill = Outcome, alpha = Significant)) +
      geom_col(position = "dodge", color = "black", size = 0.5) +
      geom_hline(yintercept = 1, linetype = "dashed", color = "red", size = 1) +
      scale_fill_manual(values = c("Baseline" = "#2E86AB", "Day 7" = "#F18F01")) +
      scale_alpha_manual(values = c("TRUE" = 1.0, "FALSE" = 0.4), guide = "none") +
      labs(
        title = "Transmission Multiplier Comparison",
        subtitle = "Effect of controlling for malaria prevalence",
        x = "Analysis Type",
        y = "Transmission Multiplier",
        fill = "Gametocyte Outcome",
        caption = "Dashed red line = no transmission advantage (1.0x)\nFaded bars = non-significant effects"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(face = "bold", size = 14),
        legend.position = "bottom"
      )

    print(p)
    ggsave("transmission_multiplier_comparison.png", p, width = 10, height = 7, dpi = 300)

    return(p)
  }
}

# Create the comparison plot
create_multiplier_comparison_plot()



# Bet-binomial model ------------------------------------------------------


install.packages("glmmTMB")

library(glmmTMB)
library(dplyr)
library(broom.mixed)

malaria_data <- malaria_data %>%
  mutate(
    baseline_positive = round((baseline_gam_prevalence / 100) * study_n),
    day3_positive = round((day3_gam_prevalence / 100) * study_n),
    day7_positive = round((day7_gam_prevalence / 100) * study_n)
  )

library(tidyr)

malaria_long <- malaria_data %>%
  pivot_longer(
    cols = c(baseline_gam_prevalence, day3_gam_prevalence, day7_gam_prevalence,
             baseline_positive, day3_positive, day7_positive),
    names_to = c("timepoint", ".value"),
    names_pattern = "(baseline|day3|day7)_(.*)"
  ) %>%
  mutate(
    timepoint = recode(timepoint, baseline = 0, day3 = 3, day7 = 7),
    timepoint = as.integer(timepoint)
  )

library(glmmTMB)

model_by_day <- function(data, day) {
  df <- data %>% filter(timepoint == day)
  glmmTMB(
    cbind(positive, study_n - positive) ~ k13_mutation_pct,
    family = betabinomial(),
    data = df
  )
}

model_day0 <- model_by_day(malaria_long, 0)
model_day3 <- model_by_day(malaria_long, 3)
model_day7 <- model_by_day(malaria_long, 7)


library(dplyr)

# Function to generate prediction data + fitted values
get_predictions <- function(model, timepoint_label) {
  new_data <- data.frame(k13_mutation_pct = seq(0, 100, by = 1))
  preds <- predict(model, newdata = new_data, se.fit = TRUE, type = "response")

  tibble(
    k13_mutation_pct = new_data$k13_mutation_pct,
    prevalence = preds$fit * 100,
    lower = (preds$fit - 1.96 * preds$se.fit) * 100,
    upper = (preds$fit + 1.96 * preds$se.fit) * 100,
    timepoint = timepoint_label
  )
}

# Get predicted curves for each model
pred_day0 <- get_predictions(model_day0, "Day 0")
pred_day3 <- get_predictions(model_day3, "Day 3")
pred_day7 <- get_predictions(model_day7, "Day 7")

predictions_all <- bind_rows(pred_day0, pred_day3, pred_day7)


library(ggplot2)

ggplot(predictions_all, aes(x = k13_mutation_pct, y = prevalence, color = timepoint)) +
  geom_line(size = 1.2) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = timepoint), alpha = 0.2, color = NA) +
  labs(
    title = "Predicted Gametocyte Prevalence by K13 Resistance",
    x = "K13 Mutation Prevalence (%)",
    y = "Predicted Gametocyte Prevalence (%)",
    color = "Timepoint",
    fill = "Timepoint"
  ) +
  theme_minimal(base_size = 14) +
  scale_color_manual(values = c("Day 0" = "#2E86AB", "Day 3" = "#A23B72", "Day 7" = "#F18F01")) +
  scale_fill_manual(values = c("Day 0" = "#2E86AB", "Day 3" = "#A23B72", "Day 7" = "#F18F01"))


ggplot() +
  geom_point(data = malaria_long, aes(x = k13_mutation_pct, y = (positive / study_n) * 100, color = factor(timepoint)), alpha = 0.5) +
  geom_line(data = predictions_all, aes(x = k13_mutation_pct, y = prevalence, color = timepoint), size = 1.2) +
  geom_ribbon(data = predictions_all, aes(x = k13_mutation_pct, ymin = lower, ymax = upper, fill = timepoint), alpha = 0.2, color = NA) +
  labs(
    title = "Observed and Predicted Gametocyte Prevalence",
    x = "K13 Mutation Prevalence (%)",
    y = "Prevalence (%)",
    color = "Timepoint",
    fill = "Timepoint"
  ) +
  theme_minimal(base_size = 14)

ggplot(predictions_all, aes(x = k13_mutation_pct, y = prevalence)) +
  geom_line(color = "#2E86AB", size = 1.2) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#2E86AB", alpha = 0.2) +
  facet_wrap(~ timepoint, ncol = 1) +
  labs(
    title = "Predicted Gametocyte Prevalence by K13 Resistance",
    x = "K13 Mutation Prevalence (%)",
    y = "Predicted Prevalence (%)"
  ) +
  theme_minimal(base_size = 14)



# Figures for write-up ----------------------------------------------------

# Create horizontal point scale legend
improved_size_scale_horizontal <- scale_size_continuous(
  name = "Sample Size",
  range = c(2.5, 11),
  breaks = c(50, 75, 100, 150, 200, 300, 500, 1000, 2000, 3000),
  labels = c("50", "75", "100", "150", "200", "300", "500", "1K", "2K", "3K"),
  trans = "sqrt",
  guide = guide_legend(
    override.aes = list(alpha = 0.8, color = "#F18F01"),
    nrow = 1,  # Single horizontal row
    title.position = "left",  # Title on the left
    label.position = "bottom",  # Numbers below the points
    direction = "horizontal",
    title.hjust = 0.5
  )
)

# Apply to your plots with horizontal legend at bottom
p_baseline_horizontal <- p_baseline_simple +
  improved_size_scale_horizontal +
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.margin = margin(t = 10),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

p_day3_horizontal <- p_day3_simple +
  improved_size_scale_horizontal +
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.margin = margin(t = 10),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

p_day7_horizontal <- p_day7_simple +
  improved_size_scale_horizontal +
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.margin = margin(t = 10)
  )

# Alternative: Custom breaks that focus more on your smaller studies
focused_size_scale <- scale_size_continuous(
  name = "Sample Size",
  range = c(3, 10),
  breaks = c(50, 75, 100, 150, 200, 300, 500, 1000, 3000),  # Removed 2000 to save space
  labels = c("50", "75", "100", "150", "200", "300", "500", "1K", "3K"),
  trans = "sqrt",
  guide = guide_legend(
    override.aes = list(alpha = 0.8, color = "#F18F01"),
    nrow = 1,
    title.position = "left",
    label.position = "bottom",
    direction = "horizontal"
  )
)

# Apply focused version
p_baseline_focused <- p_baseline_simple +
  focused_size_scale +
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

p_day3_focused <- p_day3_simple +
  focused_size_scale +
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

p_day7_focused <- p_day7_simple +
  focused_size_scale +
  theme(
    legend.position = "bottom",
    legend.box = "horizontal"
  )

# Create patchwork with individual horizontal legends
combined_horizontal <- p_baseline_horizontal / p_day3_horizontal / p_day7_horizontal +
  plot_layout(ncol = 1, guides = "keep")

# Or with focused legends (fewer points, cleaner)
combined_focused <- p_baseline_focused / p_day3_focused / p_day7_focused +
  plot_layout(ncol = 1, guides = "keep")

# View the results
print("Horizontal legend version:")
print(combined_horizontal)

print("Focused horizontal legend version:")
print(combined_focused)

# Save both versions
ggsave("combined_horizontal_legends.png", combined_horizontal,
       width = 12, height = 18, dpi = 300)
ggsave("combined_focused_legends.png", combined_focused,
       width = 12, height = 16, dpi = 300)

# If you want even more control, create a minimal version with just key sizes:
minimal_size_scale <- scale_size_continuous(
  name = "n = ",
  range = c(3, 9),
  breaks = c(60, 100, 200, 500, 1000, 3000),
  labels = c("60", "100", "200", "500", "1K", "3K"),
  trans = "sqrt",
  guide = guide_legend(
    override.aes = list(alpha = 0.8, color = "#F18F01"),
    nrow = 1,
    title.position = "left",
    label.position = "bottom",
    direction = "horizontal"
  )
)

# Quick test of minimal version
p_baseline_minimal <- p_baseline_simple +
  minimal_size_scale +
  theme(legend.position = "bottom", legend.box = "horizontal")

print("Minimal version:")
print(p_baseline_minimal)


# REMOVE ALL TITLES FROM PLOTS

# Clean versions without any titles or subtitles
p_baseline_no_titles <- p_baseline_focused +
  labs(title = NULL, subtitle = NULL) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

p_day3_no_titles <- p_day3_focused +
  labs(title = NULL, subtitle = NULL) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

p_day7_no_titles <- p_day7_focused +
  labs(title = NULL, subtitle = NULL)

# Combined plot without titles
combined_no_titles <- p_baseline_no_titles / p_day3_no_titles / p_day7_no_titles +
  plot_layout(ncol = 1, guides = "keep")

# View and save
print(combined_no_titles)
ggsave("combined_no_titles.png", combined_no_titles,
       width = 12, height = 16, dpi = 300)
