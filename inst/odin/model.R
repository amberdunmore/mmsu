# Model definitions and differential equations

# Human Equations
deriv(S) <- (-S * lambda_s * (phi * ft + phi * (1 - ft) + (1 - phi)) -
               S * lambda_r * (phi * ft + phi * (1 - ft) + (1 - phi)) +
               T_s * rT_s + A_s * rA + A_r * rA + T_r_cleared * rT_r_cleared)

deriv(D_s) <- (S * lambda_s * phi * (1 - ft) +
                lambda_s * A_r * phi * (1 - ft) +
                lambda_s * A_s * phi * (1 - ft) -
                D_s * rD) - invading_D_r

deriv(A_s) <- (S * lambda_s * (1 - phi) +
                lambda_s * A_r * (1 - phi) +
                D_s * rD -
                lambda_s * A_s * phi * (1 - ft) -
                lambda_s * A_s * phi * ft -
                lambda_r * A_s * phi * (1 - ft) -
                lambda_r * A_s * phi * ft -
                lambda_r * A_s * (1 - phi) -
                A_s * rA - invading_A_r)

deriv(T_s) <- (S * lambda_s * phi * ft +
                lambda_s * A_r * phi * ft +
                lambda_s * A_s * phi * ft -
                T_s * rT_s) - invading_T_r_cleared - invading_T_r_failed

deriv(D_r) <- invading_D_r + (S * lambda_r * phi * (1 - ft) +
                lambda_r * A_s * phi * (1 - ft) +
                lambda_r * A_r * phi * (1 - ft) -
                D_r * rD_r)

deriv(A_r) <- invading_A_r + (S * lambda_r * (1 - phi) +
                              lambda_r * A_s * (1 - phi) +
                              D_r * rD_r +
                              T_r_failed * rT_r_failed -
                              lambda_r * A_r * phi * (1 - ft) -
                              lambda_r * A_r * phi * ft -
                              lambda_s * A_r * phi * (1 - ft) -
                              lambda_s * A_r * phi * ft -
                              lambda_s * A_r * (1 - phi) -
                              A_r * rA_r)

# NEW: Split resistant treatment compartments

deriv(T_r_cleared) <- invading_T_r_cleared + (S * lambda_r * phi * ft * (1 - treatment_failure_rate) +
                lambda_r * A_r * phi * ft * (1 - treatment_failure_rate) +
                lambda_r * A_s * phi * ft * (1 - treatment_failure_rate) -
                T_r_cleared * rT_r_cleared)

deriv(T_r_failed) <- invading_T_r_failed + (S * lambda_r * phi * ft * treatment_failure_rate +
                lambda_r * A_r * phi * ft * treatment_failure_rate +
                lambda_r * A_s * phi * ft * treatment_failure_rate -
                T_r_failed * rT_r_failed)


# Mosquito Equations
deriv(Sv) <- mu - (lambda_v_s + lambda_v_r) * Sv - mu * Sv

delayed_lambda_v_s_Sv <- delay(lambda_v_s * Sv * exp(-mu * n), n)
deriv(Ev_s) <- lambda_v_s * Sv - delayed_lambda_v_s_Sv - mu * Ev_s - invading_Ev_r
deriv(Iv_s) <- delayed_lambda_v_s_Sv - mu * Iv_s - invading_Iv_r

delayed_lambda_v_r_Sv <- delay(lambda_v_r * Sv * exp(-mu * n), n)
deriv(Ev_r) <- invading_Ev_r + (lambda_v_r * Sv - delayed_lambda_v_r_Sv - mu * Ev_r)
deriv(Iv_r) <- invading_Iv_r + (delayed_lambda_v_r_Sv - mu * Iv_r)


# Outputs
output(prevalence) <- A_s + D_s + T_s + A_r + D_r + T_r_cleared + T_r_failed
output(prevalence_res) <- (A_r + D_r + T_r_cleared + T_r_failed) / (A_s + D_s + T_s + A_r + D_r + T_r_cleared + T_r_failed)
output(population) <- S + D_s + A_s + T_s + D_r + A_r + T_r_cleared + T_r_failed
output(population_v) <- Sv + Ev_s + Iv_s + Ev_r + Iv_r
output(prevalence_sensitive) <- A_s + D_s + T_s
output(prevalence_resistant_treated) <- T_r_cleared + T_r_failed
output(prevalence_resistant_cleared) <- T_r_cleared
output(prevalence_resistant_failed) <- T_r_failed

# EIR calculations
EIR_s <- m * a * Iv_s * 365
EIR_r <- m * a * Iv_r * 365
output(EIR_s) <- EIR_s
output(EIR_r) <- EIR_r


#EIR_global <- (1 - resistant_ratio) * EIR_s + resistant_ratio * EIR_r
EIR_global <- EIR_s + EIR_r
output(EIR_global) <- EIR_global

# Resistance introduction
invading_A_r <- if(t < res_time || t > (res_time+1)) 0 else A_s*log(1/(1-init_res))
invading_T_r_cleared <- if(t < res_time || t > (res_time+1)) 0 else T_s*log(1/(1-init_res))*(1-treatment_failure_rate)
invading_T_r_failed <- if(t < res_time || t > (res_time+1)) 0 else T_s*log(1/(1-init_res))*treatment_failure_rate
invading_D_r <- if(t < res_time || t > (res_time+1)) 0 else D_s*log(1/(1-init_res))
invading_Ev_r <- if(t < res_time || t > (res_time+1)) 0 else Ev_s*log(1/(1-init_res))
invading_Iv_r <- if(t < res_time || t > (res_time+1)) 0 else Iv_s*log(1/(1-init_res))
output(invading_A_r_out) <- invading_A_r

# Initial conditions
initial(S) <- S0
initial(D_s) <- D_s0
initial(A_s) <- A_s0
initial(T_s) <- T_s0
initial(D_r) <- D_r0
initial(A_r) <- A_r0
initial(T_r_cleared) <- T_r_cleared0
initial(T_r_failed) <- T_r_failed0
initial(Sv) <- Sv0
initial(Ev_s) <- Ev_s0
initial(Iv_s) <- Iv_s0
initial(Ev_r) <- Ev_r0
initial(Iv_r) <- Iv_r0

# User-defined parameters
S0 <- user()
D_s0 <- user()
A_s0 <- user()
T_s0 <- user()
D_r0 <- user()
A_r0 <- user()
T_r_cleared0 <- user()
T_r_failed0 <- user()
m <- user()
a <- user()
b <- user()
lambda_s <- m * a * b * Iv_s
lambda_r <- m * a * b * Iv_r * resistance_trans_mult
lambda_v_s <- a * (cA_s * A_s + cD_s * D_s + cT_s * T_s) # Updated to incorporate separate infectiousness
lambda_v_r <- a * (cA_r * A_r + cD_r * D_r + cT_r_cleared * T_r_cleared + cT_r_failed * T_r_failed)
phi <- user()
ft <- user()
rD <- user()
rA <- user()
rT_s <- user()
rT_r_cleared <- user()
rT_r_failed <- user()
Sv0 <- user()
Ev_s0 <- user()
Iv_s0 <- user()
Ev_r0 <- user()
Iv_r0 <- user()
mu <- user()
n <- user()
cA <- user()
cD <- user()
cT <- user()
ton <- user()
toff <- user()
res_time <- user()
init_res <- user()


# Treatment failure and resistance parameters
treatment_failure_rate <- user()
resistance_trans_mult <- user()
resistance_dur_mult <- user()
resistance_baseline_ratio <- user()
resistance_cleared_ratio <- user()
resistance_failed_ratio <- user()

# Defining separate infectiousness parameters for sensitive and resistant strains
cA_s <- cA  # baseline asymptomatic infectiousness (sensitive)
cD_s <- cD  # baseline diseased infectiousness (sensitive)
cT_s <- cT  # baseline treated infectiousness (sensitive)

# Resistant infectiousness parameters (time-dependent activation)
cA_r <- cA * if (t > ton && t < toff) resistance_baseline_ratio else 1
cD_r <- cD * if (t > ton && t < toff) resistance_baseline_ratio else 1
cT_r_cleared <- cT * if (t > ton && t < toff) resistance_cleared_ratio else 1
cT_r_failed <- cT * if (t > ton && t < toff) resistance_failed_ratio else 1


# FIXED: Actually use resistance_dur_mult in the resistant recovery rates
rD_r <- rD / if (t > ton && t < toff) resistance_dur_mult else 1   # symptomatic (untreated) infection duration for resistant
rA_r <- rA / if (t > ton && t < toff) resistance_dur_mult else 1   # asymptomatic infection duration for resistant
