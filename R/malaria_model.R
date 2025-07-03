#' Create Malaria Model
#'
#' @param params A list of parameters to pass to the model. If NULL, parameters will be generated based on EIR and ft.
#' @param EIR Entomological Inoculation Rate (used if params is NULL)
#' @param ft Treatment rate (used if params is NULL)
#' @param ton Time at which resistance is turned on
#' @param toff Time at which resistance is turned off
#' @param init_res Initial resistance level at res_time
#' @param day0_res Resistant at Day 0. Default = 0.01
#' @param res_time Time at which resistance is introduced
#' @param rT_r_true True treatment rate for resistant parasites
#' @param resistance_trans_mult transmission multiplier for resistant parasites (default = 1)
#' @param resistance_dur_mult duration multiplier for resistant infections (default = 1)
#' @param resistance_baseline_ratio baseline infectiousness ratio for untreated resistant parasites (default = 1)
#' @param resistance_treated_ratio post-treatment infectiousness ratio for treated resistant parasites (default = 1)
#' @param verbose Logical. If TRUE, prints detailed logs. Default is FALSE.
#' @return An object of class `odin_model`.
#' @export
malaria_model <- function(params = NULL, EIR = NULL, ft = NULL,
                          ton = 365, toff = 4015, day0_res = 0.01,
                          init_res = 0.0, res_time = 0, rT_r_true = 0.0,
                          resistance_trans_mult = 1, resistance_dur_mult = 1,
                          resistance_baseline_ratio = 1, resistance_treated_ratio = 1,
                          verbose = FALSE) {

    if (is.null(params)) {
      if (is.null(EIR) || is.null(ft)) {
        stop("Either params or both EIR and ft must be provided")
      }
      params <- phi_eir_rel(EIR, ft, ton, toff, init_res, res_time, rT_r_true,
                            resistance_trans_mult = resistance_trans_mult,
                            resistance_dur_mult = resistance_dur_mult,
                            resistance_baseline_ratio = resistance_baseline_ratio,
                            resistance_treated_ratio = resistance_treated_ratio)
    }

    # Update parameters
    params$ton <- ton
    params$toff <- toff
    params$init_res <- init_res
    params$res_time <- res_time
    params$rT_r_true <- rT_r_true
    params$resistance_trans_mult <- resistance_trans_mult  # add transmission multiplier (default 1 if not provided)
    params$resistance_dur_mult <- resistance_dur_mult      # add duration multiplier (default 1 if not provided)
    params$resistance_baseline_ratio <- resistance_baseline_ratio
    params$resistance_treated_ratio <- resistance_treated_ratio

    # Generate initial parameters if not already present
    if (!"S0" %in% names(params)) {
      params$S0 <- params$S
      params$D_s0 <- params$D * (1 - day0_res)
      params$D_r0 <- params$D * day0_res
      params$A_s0 <- params$A * (1 - day0_res)
      params$A_r0 <- params$A * day0_res
      params$T_s0 <- params$T * (1 - day0_res)
      params$T_r0 <- params$T * day0_res
      params$Sv0 <- params$Sv
      params$Ev_s0 <- params$Ev * (1 - day0_res)
      params$Iv_s0 <- params$Iv * (1 - day0_res)
      params$Ev_r0 <- params$Ev * day0_res
      params$Iv_r0 <- params$Iv * day0_res
    }

    # Ensure all parameters are scalar
    params <- lapply(params, function(x) if(length(x) > 1) x[1] else x)

    # Check if all required parameters are present and numeric
    required_params <- c("S0", "D_s0", "A_s0", "T_s0", "D_r0", "A_r0", "T_r0",
                         "Sv0", "Ev_s0", "Iv_s0", "Ev_r0", "Iv_r0",
                         "m", "a", "b", "phi", "ft", "rD", "rA", "rT_s", "rT_r",
                         "mu", "n", "cA", "cD", "cT",
                         "ton", "toff", "res_time", "init_res",
                         "resistance_trans_mult", "resistance_dur_mult",
                         "resistance_baseline_ratio", "resistance_treated_ratio")
    missing_params <- setdiff(required_params, names(params))
    if (length(missing_params) > 0) {
      stop("Missing required parameters: ", paste(missing_params, collapse = ", "))
    }

    # Calculate EIR if not provided
    if (!"EIR" %in% names(params)) {
      params$EIR <- params$m * params$a * params$Iv_s0 * 365  # Annual EIR
    }

    non_numeric_params <- names(params)[!sapply(params, is.numeric)]
    if (length(non_numeric_params) > 0) {
      stop("Non-numeric parameters: ", paste(non_numeric_params, collapse = ", "))
    }

    if (verbose) {
      cat("Creating odin model...\n")
    }

    if (verbose) {
      cat("Initializing model with parameters...\n")
    }

    # create our model
    model_instance <- model$new(user = params, unused_user_action = "ignore")

    if (verbose) {
      cat("Model initialized successfully.\n")
    }

    return(model_instance)

}

# Helper function to provide a default value
`%||%` <- function(x, y) if (is.null(x)) y else x
