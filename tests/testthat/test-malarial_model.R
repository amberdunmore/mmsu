# tests/testthat/test-malarial_model.R

test_that("malaria_model runs without errors and produces expected output", {
  params <- list(S0 = 0.9, Ds0 = 0.05, As0 = 0.02, Ts0 = 0.03, DR0 = 0, AR0 = 0,
                 TR0 = 0, Sv0 = 0.9, Ev_s0 = 0.08, Iv_s0 = 0.02, Ev_r0 = 0, Iv_r0 = 0,
                 m = 1.271483, a = 0.3, b = 0.5876259, Phi = 0.7, fT = 0.1, rD = 0.2,
                 rA = 0.005, rTs = 0.2, rTR_true = 0.01, e = 0.132, mu = 0.132,
                 n = 10, c_A = 0.05, c_D = 0.06, c_T = 0.02, ton = 10000,
                 toff = 20000, res_time = 100)

  model <- malaria_model(params)
  times <- seq(0, 5000, by = 1)
  model_results <- model$run(times)

  expect_type(model_results, "list")
  expect_true(all(c("t", "prevalence", "prevalence_res", "population") %in% names(model_results)))
})

plot_malaria_model <- function(model_results, res_time) {
  results_df <- data.frame(
    time = model_results[, "t"],
    prevalence = model_results[, "prevalence"],
    prevalence_res = model_results[, "prevalence_res"],
    N = model_results[, "population"]
  )

  p <- ggplot(results_df, aes(x = time)) +
    geom_line(aes(y = prevalence, color = "Prevalence"), linewidth = 1.5) +
    geom_line(aes(y = prevalence_res, color = "Prevalence (Resistant)"), linewidth = 1.5) +
    geom_line(aes(y = N, color = "Population"), linetype = "dashed", linewidth = 1.5) +
    geom_vline(xintercept = res_time, linetype = "dashed", color = "#000000", linewidth = 1.5) +
    labs(x = "Time", y = "Value", color = "Variable",
         title = "Malaria Model Results",
         subtitle = paste("Resistance introduced at time", res_time)) +
    scale_color_manual(values = c("Prevalence" = "#FE938C", "Prevalence (Resistant)" = "#5CA4A9", "Population" = "#B8BEDD")) +
    theme_minimal()

  print(p)
}

test_that("plot_malaria_model generates a plot", {
  params <- list(S0 = 0.9, Ds0 = 0.05, As0 = 0.02, Ts0 = 0.03, DR0 = 0, AR0 = 0,
                 TR0 = 0, Sv0 = 0.9, Ev_s0 = 0.08, Iv_s0 = 0.02, Ev_r0 = 0, Iv_r0 = 0,
                 m = 1.271483, a = 0.3, b = 0.5876259, Phi = 0.7, fT = 0.1, rD = 0.2,
                 rA = 0.005, rTs = 0.2, rTR_true = 0.01, e = 0.132, mu = 0.132,
                 n = 10, c_A = 0.05, c_D = 0.06, c_T = 0.02, ton = 10000,
                 toff = 20000, res_time = 100)

  model <- malaria_model(params)
  times <- seq(0, 5000, by = 1)
  model_results <- model$run(times)

  expect_silent({
    plot_malaria_model(model_results, res_time = 100)
  })
})
