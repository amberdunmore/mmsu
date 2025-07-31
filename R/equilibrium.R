#' @importFrom magrittr %>%
#' @importFrom dplyr bind_rows
#' @importFrom tidyr expand_grid
#' @importFrom stats weighted.mean
NULL

#' @import ICDMM
NULL

#' Create equilibrium initial conditions
#'
#' @param par A list of parameters
#' @return A data frame of equilibrium conditions
#' @export
equilibrium_init_create <- function(par) {
  ## EIR
  EIRY_eq <- par$EIR  # initial annual EIR
  EIRd_eq <- EIR_eq <- EIRY_eq/365

  # FOI
  FOI_eq <- EIR_eq * par$b

  # FOI to T_cleared and T_failed (split based on failure rate)
  aT_cleared <- FOI_eq * par$phi * par$ft * (1 - par$treatment_failure_rate) / par$rT_r_cleared
  aT_failed <- FOI_eq * par$phi * par$ft * par$treatment_failure_rate / par$rT_r_failed
  aD <- FOI_eq * par$phi * (1 - par$ft) / par$rD

  # Updated equilibrium calculations
  Z_eq <- rep(NA, 4)  # Now 4 states: S, D, T_cleared, T_failed
  Z_eq[1] <- 1/(1 + aT_cleared + aT_failed + aD)  # S
  Z_eq[2] <- aD * Z_eq[1]  # D
  Z_eq[3] <- aT_cleared * Z_eq[1]  # T_cleared
  Z_eq[4] <- aT_failed * Z_eq[1]  # T_failed

  Y_eq <- Z_eq[1]
  D_eq <- Z_eq[2]
  T_cleared_eq <- Z_eq[3]
  T_failed_eq <- Z_eq[4]

  betaS <- FOI_eq
  betaA <- FOI_eq * par$phi + par$rA

  # A includes flows from failed treatments (not cleared)
  # T_r_cleared flows back to S, so doesn't contribute to A equilibrium
  A_eq <- (FOI_eq * (1 - par$phi) * Y_eq + par$rD * D_eq +
             par$rT_r_failed * T_failed_eq) / (betaA + FOI_eq * (1 - par$phi))

  # S equilibrium includes people recovering from T_r_cleared
  S_eq <- Y_eq - A_eq + par$rT_r_cleared * T_cleared_eq / (FOI_eq)

  # Updated FOIv calculation including both treatment compartments
  FOIv_eq <- par$a * (par$cT * (T_cleared_eq + T_failed_eq) + par$cD * D_eq + par$cA * A_eq)

  # mosquito states
  Iv_eq <- FOIv_eq * exp(-par$mu * par$n)/(FOIv_eq + par$mu)
  Sv_eq <- par$mu * Iv_eq/(FOIv_eq * exp(-par$mu * par$n))
  Ev_eq <- 1 - Sv_eq - Iv_eq

  # mosquito density needed to give this EIR
  mv0 <- EIRd_eq/(Iv_eq * par$a)

  ## collate init
  list(
    EIR = par$EIR, ft = par$ft,
    S = S_eq, D = D_eq, A = A_eq,
    T_cleared = T_cleared_eq, T_failed = T_failed_eq,
    phi = par$phi, b = par$b,
    m = mv0, Sv = Sv_eq, Ev = Ev_eq, Iv = Iv_eq, a = par$a,
    cA = par$cA, cD = par$cD, cT = par$cT,
    n = par$n,
    mu = par$mu,
    rD = par$rD,
    rA = 1/(1/par$rA +1/par$rU),
    rT_s = par$rT_s,
    rT_r_cleared = par$rT_r_cleared,
    rT_r_failed = par$rT_r_failed,
    treatment_failure_rate = par$treatment_failure_rate,
    resistance_trans_mult = ifelse(is.null(par$resistance_trans_mult), 1, par$resistance_trans_mult),
    resistance_dur_mult = ifelse(is.null(par$resistance_dur_mult), 1, par$resistance_dur_mult),
    resistance_baseline_ratio = ifelse(is.null(par$resistance_baseline_ratio), 1, par$resistance_baseline_ratio),
    resistance_cleared_ratio = ifelse(is.null(par$resistance_cleared_ratio), 1, par$resistance_cleared_ratio),
    resistance_failed_ratio = ifelse(is.null(par$resistance_failed_ratio), 1, par$resistance_failed_ratio)
  ) %>%
    as.data.frame()
}

#' Generate equilibrium parameters based on EIR and ft
#'
#' @param EIR Numeric. The Entomological Inoculation Rate.
#' @param ft Numeric. The treatment rate.
#' @param ton Time at which treatment is turned on
#' @param toff Time at which treatment is turned off
#' @param day0_res Resistant at Day 0. Default = 0
#' @param init_res Initial resistance level at res_time
#' @param res_time Time at which resistance is introduced
#' @param treatment_failure_rate Proportion of treatments that fail completely (default = 0.3)
#' @param rT_r_cleared Recovery rate for delayed clearance treatment (default = 0.1)
#' @param rT_r_failed Rate of transition from failed treatment to asymptomatic (default = 0.1)
#' @param resistance_trans_mult transmission multiplier for resistant parasites (default = 1)
#' @param resistance_dur_mult duration multiplier for resistant infections (default = 1)
#' @param resistance_baseline_ratio baseline infectiousness ratio for untreated resistant parasites (default = 1)
#' @param resistance_cleared_ratio infectiousness ratio for partial treatment response (default = 1)
#' @param resistance_failed_ratio infectiousness ratio for treatment failure (default = 1)
#' @return A list of generated parameters.
#' @export
phi_eir_rel <- function(EIR, ft, ton = 5000, toff = 50000, init_res = 0.01, res_time = 3000,
                        treatment_failure_rate = 0.3, rT_r_cleared = 0.1, rT_r_failed = 0.1,
                        day0_res = 0.01, resistance_trans_mult = 1, resistance_dur_mult = 1,
                        resistance_baseline_ratio = 1, resistance_cleared_ratio = 1, resistance_failed_ratio = 1) {
  mpl <- ICDMM::model_param_list_create(rho=0, rA = 1/(250), rU = Inf, rP = Inf, sigma2 = 0)
  eq <- ICDMM::equilibrium_init_create(
    age_vector=c(0,0.25,0.5,0.75,1,1.25,1.5,1.75,2,3.5,5,7.5,10,15,20,30,40,50,60,70,80),
    EIR=EIR, ft=ft,
    model_param_list = mpl, het_brackets=2,
    country = NULL,
    admin_unit = NULL)

  # Safe function to calculate weighted mean
  phi <- weighted.mean(
    apply(eq$phi_eq, 1, weighted.mean, eq$het_wt),
    eq$den
  )

  c_A <- weighted.mean(
    apply(eq$cA_eq, 1, weighted.mean, eq$het_wt),
    rowMeans(eq$init_A)
  )

  b <- weighted.mean(rowMeans(eq$b0 * ((1 - eq$b1)/(1 + (eq$init_IB/eq$IB0)^eq$kB) + eq$b1)), eq$den)

  S <- sum(eq$init_S) + sum(eq$init_P)
  D <- sum(eq$init_D)
  A <- sum(eq$init_A + eq$init_U)
  T_total <- sum(eq$init_T)

  # Initial populations (split treatments by failure rate)
  T_cleared <- T_total * (1 - treatment_failure_rate)
  T_failed <- T_total * treatment_failure_rate

  lambda_v_scale <- ((eq$av0 * (c_A*A + eq$cD*D + eq$cT*(T_cleared + T_failed)))/eq$FOIv_eq)

  par <- list(
    EIR = EIR, ft = ft,
    S = S, D = D, A = A, T_cleared = T_cleared, T_failed = T_failed, phi = phi, b = b,
    m = eq$mv0, Sv = eq$init_Sv, Ev = eq$init_Ev, Iv = eq$init_Iv, a = eq$av0,
    cA = c_A, cD = mean(eq$cD, na.rm = TRUE), cT = mean(eq$cT, na.rm = TRUE),
    n = eq$delayMos,
    mu = eq$mu0,
    rD = eq$rD,
    rU = eq$rU,
    rA = eq$rA,
    rT_s = eq$rT,
    rT_r_cleared = rT_r_cleared,
    rT_r_failed = rT_r_failed,
    lambda_v_scale = lambda_v_scale,
    ton = ton,
    toff = toff,
    res_time = res_time,
    day0_res = day0_res,
    init_res = init_res,
    treatment_failure_rate = treatment_failure_rate,
    resistance_trans_mult = resistance_trans_mult,   # include user-specified resistance multipliers
    resistance_dur_mult = resistance_dur_mult,       # (defaults are 1 if not provided)
    resistance_baseline_ratio = resistance_baseline_ratio,
    resistance_cleared_ratio = resistance_cleared_ratio,
    resistance_failed_ratio = resistance_failed_ratio
  )

  equilibrium_init_create(par)
}

#' Generate starting parameters for a range of EIR and ft values
#'
#' @param EIRs Numeric vector. The range of Entomological Inoculation Rates.
#' @param fts Numeric vector. The range of treatment rates.
#' @return A data frame of starting parameters for each combination of EIR and ft.
#' @export
starting_params <- function(EIRs = c(0.25, 0.5, 0.75, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10, 15, 20, 30, 50, 100, 200),
                            fts = seq(0.1, 0.9, 0.1)) {
  pars <- expand_grid("EIR" = EIRs, "ft" = fts)
  pars$n <- seq_along(pars$EIR)

  starting_params <- lapply(split(pars, pars$n),
                            function(x) {
                              phi_eir_rel(x$EIR, x$ft)
                            }) %>% bind_rows()

  return(starting_params)
}
