# 01_explore_mmsu_model.R
# Exploring the mmsu malaria model for systematic review

# Load the package
library(mmsu)

# =============================================================================
# 1. EXPLORE AVAILABLE FUNCTIONS
# =============================================================================

# See all exported functions
ls("package:mmsu")

# Key functions available:
# - malaria_model(): Creates the main model
# - phi_eir_rel(): Generates equilibrium parameters
# - plot_model(): Visualizes results
# - equilibrium_init_create(): Creates initial conditions
# - starting_params(): Generate parameter sets

# =============================================================================
# 2. BASIC MODEL EXPLORATION
# =============================================================================

# Create a basic model
model <- malaria_model(EIR = 10, ft = 0.5)
print("Model created successfully!")

# Run for 1 year
results <- model$run(0:365)
print("Model run completed!")

# Basic visualization
plot_model(results, c("S", "prevalence"))

# =============================================================================
# 3. UNDERSTAND MODEL STRUCTURE
# =============================================================================

# Look at what the model outputs
print("Model output variables:")
print(colnames(results))

# Key variables for your systematic review:
# - S: Susceptible individuals
# - prevalence: Total prevalence
# - prevalence_res: Resistant prevalence
# - EIR_s, EIR_r: Transmission rates (sensitive/resistant)

# =============================================================================
# 4. PARAMETER EXPLORATION
# =============================================================================

# Generate equilibrium parameters
params <- phi_eir_rel(EIR = 10, ft = 0.5)
print("Parameter structure:")
str(params)

# Continue with your exploration...
