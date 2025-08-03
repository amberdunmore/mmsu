
<!-- README.md is generated from README.Rmd. Please edit that file -->

# mmsu

**mmsu** is an R package designed to create and run deterministic models
of the spread of antimalarial drug resistance. This package provides a
range of functionalities to simulate malaria transmission dynamics and
visualize the results.

## Installation

To install the package, you can use the `devtools` package to install
directly from GitHub:

``` r
# devtools::install_github("OJWatson/mmsu")
library(mmsu)
```

## Dependencies

The package imports several other R packages, including: **odin**,
**ICDMM**, **tidyr**, **stats**, **ggplot2**, **data.table**, **utils**,
**dplyr**, **Rcpp**, **magrittr**. These packages will be installed
automatically when you install the **mmsu** package.

for ICDMM though please install this branch:

`devtools::install_github("mrc-ide/deterministic-malaria-model", ref = "fix/compilers")`

## Usage

### Creating a Malaria Model

The **malaria_model** function is used to create a model of malaria
transmission dynamics. You can either provide a list of parameters
directly or generate parameters based on the Entomological Inoculation
Rate (EIR) and treatment rate (ft).

#### Example 1: Create a model using EIR and treatment rate

``` r
# Create the model using EIR and ft
model <- malaria_model(EIR = 10, ft = 0.5)

# Run the model for 1000 days
results <- model$run(0:1000)

# Plot the results
plot <- plot_model(results, c("S", "prevalence"))
```

<img src="man/figures/README-unnamed-chunk-3-1.png" width="100%" />

The model is set up to start at steady state equilibrium for the
parameters provided, which is why outputs are flat over time.

#### Example 2: Create a model by specfiying all parameters and initial conditions

If we provide parameters that are not the steady state equilibrium the
model will then be seen to move towards this. In general, we will always
intialise a model using just EIR and ft so that we do not have to wait
for the model to reach equilibrium.

``` r

# Define parameters
params <- list(
  S = 1000, D = 10, A = 5, T = 3, Sv = 200, Ev = 10, Iv = 5,
  m = 0.1, a = 0.3, b = 0.55, phi = 0.3, ft = 0.5, rD = 0.1, rA = 0.1,
  rT_s = 0.1, rT_r = 0.1, e = 0.2, mu = 0.1, n = 10, cA = 0.1, cD = 0.2, cT = 0.3,
  ton = 5000, toff = 50000, res_time = 3000, init_res = 0.01
)

# Create the model
model <- malaria_model(params = params)

# Run the model for 1000 days
results <- model$run(0:1000)

# Plot the results
plot <- plot_model(results, c("S", "prevalence", "prevalence_res"))
```

<img src="man/figures/README-unnamed-chunk-4-1.png" width="100%" />

### Running the model to look at impact of resistance

The dafults set the model up with the effect of resistant parasites not
being modelled. Currently, resistant parasites have an extra parameter
which governs the speed at which a treated indiviudal recovers. If this
is slower than wild type parasites we would expect a selective benefit
to occur and the model system to move to a new equilibrium.

#### Example: Resistance 100 times slower treatment rates

We set up the model by passing in what the initial resistance prevalence
is at t = 0, what the treatment rate is for resistance when we turn it
on, and the time we turn on the effect of resitsance (ton). We use this
so that we can check that our model is at steady state before this and
we can more easily see the impact of resistance when it is turned on

``` r

# Create the model using EIR and ft
model <- malaria_model(EIR = 10, ft = 0.5, day0_res = 0.1, ton = 500)

# Run the model for 5000 days
results <- model$run(0:5000)

# Plot the results
plot <- plot_model(results, c("S", "prevalence", "prevalence_res"))
```

<img src="man/figures/README-unnamed-chunk-5-1.png" width="100%" />

### Creating Equilibrium Initial Conditions

We have discussed the concept of the equilibrium state and solution. The
malaria model has been created so that it has a steady state
equilibrium. To see a bit more about how this is created internally by
some of the funcitons above.

#### Example: Generate parameters based on EIR and ft using `phi_eir_rel`

``` r
# Generate parameters
params <- phi_eir_rel(EIR = 10, ft = 0.5)

# Print generated parameters
print(params)
#>   EIR  ft         S           D         A           T      phi        b
#> 1  10 0.5 0.3142287 0.004491399 0.6767885 0.004491399 0.157997 0.418798
#>          m        Sv         Ev        Iv         a         cA        cD
#> 1 5.508219 0.9392847 0.04449609 0.0162192 0.3066667 0.04051686 0.0676909
#>           cT  n    mu  rD    rA rT_s treatment_failure_rate
#> 1 0.02179647 10 0.132 0.2 0.004  0.2                    0.3
#>   resistance_trans_mult resistance_dur_mult resistance_baseline_ratio
#> 1                     1                   1                         1
#>   resistance_cleared_ratio resistance_failed_ratio
#> 1                        1                       1
```

## Contributing

Contributions to the package are welcome. You can contribute by
reporting issues, suggesting new features, or submitting pull requests.

## License

This package is licensed under the GPL-3 license.
