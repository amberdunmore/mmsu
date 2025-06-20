# Model definitions and differential equations

# Human Equations
deriv(S) <- (-S * lambda_s * (phi * ft + phi * (1 - ft) + (1 - phi)) -
               S * lambda_r * (phi * ft + phi * (1 - ft) + (1 - phi)) +
               T_s * rT_s + A_s * rA + A_r * rA + T_r * rT_r)

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
                T_s * rT_s) - invading_T_r

deriv(D_r) <- invading_D_r + (S * lambda_r * phi * (1 - ft) +
                lambda_r * A_s * phi * (1 - ft) +
                lambda_r * A_r * phi * (1 - ft) -
                D_r * rD)

deriv(A_r) <- invading_A_r + (S * lambda_r * (1 - phi) +
                              lambda_r * A_s * (1 - phi) +
                              D_r * rD -
                              lambda_r * A_r * phi * (1 - ft) -
                              lambda_r * A_r * phi * ft -
                              lambda_s * A_r * phi * (1 - ft) -
                              lambda_s * A_r * phi * ft -
                              lambda_s * A_r * (1 - phi) -
                              A_r * rA)

deriv(T_r) <- invading_T_r + (S * lambda_r * phi * ft +
                lambda_r * A_r * phi * ft +
                lambda_r * A_s * phi * ft -
                T_r * rT_r)

# Mosquito Equations
deriv(Sv) <- mu - (lambda_v_s + lambda_v_r) * Sv - mu * Sv

delayed_lambda_v_s_Sv <- delay(lambda_v_s * Sv * exp(-mu * n), n)
deriv(Ev_s) <- lambda_v_s * Sv - delayed_lambda_v_s_Sv - mu * Ev_s - invading_Ev_r
deriv(Iv_s) <- delayed_lambda_v_s_Sv - mu * Iv_s - invading_Iv_r

delayed_lambda_v_r_Sv <- delay(lambda_v_r * Sv * exp(-mu * n), n)
deriv(Ev_r) <- invading_Ev_r + (lambda_v_r * Sv - delayed_lambda_v_r_Sv - mu * Ev_r)
deriv(Iv_r) <- invading_Iv_r + (delayed_lambda_v_r_Sv - mu * Iv_r)


# Outputs
output(prevalence) <- A_s + D_s + T_s + A_r + D_r + T_r
output(prevalence_res) <- (A_r + D_r + T_r) / (A_s + D_s + T_s + A_r + D_r + T_r)
output(population) <- S + D_s + A_s + T_s + D_r + A_r + T_r
output(population_v) <- Sv + Ev_s + Iv_s + Ev_r + Iv_r
output(prevalence_sensitive) <- A_s + D_s + T_s

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
invading_T_r <- if(t < res_time || t > (res_time+1)) 0 else T_s*log(1/(1-init_res))
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
initial(T_r) <- T_r0
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
T_r0 <- user()
m <- user()
a <- user()
b <- user()
lambda_s <- m * a * b * Iv_s
lambda_r <- m * a * b * Iv_r * resistance_trans_multdev
lambda_v_s <- a * (cA * A_s + cD * D_s + cT * T_s)
lambda_v_r <- a * (cA * A_r + cD * D_r + cT * T_r)
phi <- user()
ft <- user()
rD <- user()
rA <- user()
rT_s <- user()
rT_r_true <- user()
rT_r <- if (t > ton && t < toff) rT_r_true else (rT_s / resistance_dur_mult)
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


# Sexual commitment/transmission advantage parameters
resistance_trans_mult <- user()  # Transmission advantage multiplier for resistant strain
resistance_dur_mult <- user()    # Duration advantage multiplier for resistant strain
