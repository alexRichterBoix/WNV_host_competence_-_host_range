#This R script analyzes host competence for West Nile Virus (WNV) by modeling viremia dynamics, vector transmission efficiency, and host survival, with the goal of estimating species-specific transmission potential. It fits mathematical models to experimental data, calculates summary metrics, and generates plots and Excel outputs for further analysis.

# STEP 1: Fit Logistic Curve for Vector Transmission Efficiency
#     Fit two logistic models using experimental Culex pipiens data.
#     Compare models using RSS and AIC.

# STEP 2: Get the Host Viremic Curves
#     Fit an exponential growth-and-decay model per species.
#     Calculate key metrics:
#     Peak titer & day
#     Duration of viremia (titer > 4)
#     Viral load (AUC)
#     Tmax, Vmax
#     Transmission potential (via logistic model)
#     Save results to Excel.
# STEP 2.1: Plot One Species Estimated Curve
# STEP 2.2: Plot Curves at the Family or Order Level

# STEP 3: Calculate the Daily Survival Ratio
#     Compute daily survival (Sample_Size_Survival / Sample_Size_Original).
#     Calculate weighted mean survival ratio per species/day.
#     Aggregate and summarize by species.
# STEP 3.1: Plot Daily Survival Rate
#     With or without confidence intervals
#     Optionally grouped by host family
#     Optionally filtered for selected species/families

#STEP 4: Calculate Daily Transmission Efficiency (Corrected for Mortality)
#   Predict daily transmission efficiency from titers.
#   Combine with survival probability and apply infectious threshold.
#   Compute corrected daily transmission efficiency per species.


# load necessary packages
library(writexl)
library(readxl)
library(dplyr)
library(ggplot2)
library(nlstools)
library(caper)
library(lme4)
library(minpack.lm)
library(gridExtra)
library(boot)

#****************************************************************************************
#****************************************************************************************
#****************************************************************************************
#**** STEP 1
#****# FIT THE LOGISTIC CURVE FOR VECTOR TRANSMISSION EFFICIENCY = number of mosquitoes with the virus in the head or saliva glands after feeding at different doses.

# Import Data for Culex pipiens experiments using a titer gradient
data_culex <- data.frame(
  dose_trans = c(8.3, 5.95, 5.6, 6.65, 7.7, 4.9, 7.2, 6.1, 5.25, 4.2, 5.3, 4.55, 5.7, 3.85, 4.7, 1.75, 3.15, 3.5),
  transmission_efficiency = c(0.67, 0.58, 0.54, 0.47, 0.44, 0.43, 0.4, 0.34, 0.32, 0.18, 0.12, 0.08, 0.06, 0.04, 0.04, 0, 0, 0)
)

# Define the first logistic function: a traditional logistic curve with an S-shaped sigmoidal growth, symmetrical around the inflection point, and ensure y = 0 when x < or = to 0
logistic1 <- function(x, L, x0, k) {
  L / (1 + exp(-k * (x - x0))) * (x > 0)
}

# Define the modified second logistic function: adds a linear dependence on X in the numerator, making the curve flatter near X = 0; the initial growth is less steep but can ramp up faster for larger X. It is suitable when the response variable grows slowly at small X but increases rapidly.
logistic2 <- function(x, L, x0, k) {
  L * (x / (1 + exp(-k * (x - x0))))
}

# Fit the first logistic curve
fit1 <- nls(
  transmission_efficiency ~ logistic1(dose_trans, L, x0, k),
  data = data_culex,
  start = list(L = 1, x0 = 5, k = 1)
)

# Fit the modified second logistic curve
fit2 <- nls(
  transmission_efficiency ~ logistic2(dose_trans, L, x0, k),
  data = data_culex,
  start = list(L = 1, x0 = 5, k = 1)
)

# Extract parameters
fit1_params <- coef(fit1)
fit2_params <- coef(fit2)

# Generate smoother predictions for plotting
x_values <- seq(0, max(data_culex$dose_trans), length.out = 100)
curve1 <- logistic1(x_values, fit1_params["L"], fit1_params["x0"], fit1_params["k"])
curve2 <- logistic2(x_values, fit2_params["L"], fit2_params["x0"], fit2_params["k"])

# Combine data for plotting
plot_data <- data.frame(
  x = rep(x_values, 2),
  y = c(curve1, curve2),
  curve = rep(c("Logistic 1", "Modified Logistic 2"), each = length(x_values))
)

# Plot the data and fitted curves
ggplot(data_culex, aes(x = dose_trans, y = transmission_efficiency)) +
  geom_point(color = "black", size = 2) +
  geom_line(data = plot_data, aes(x = x, y = y, color = curve), linewidth = 1) +
  labs(
    title = "Comparison of Logistic Fits for Culex pipiens",
    x = "Dose (dose_trans)",
    y = "Transmission Efficiency"
  ) +
  scale_color_manual(values = c("red", "blue"), labels = c("Fit 1", "Modified Fit 2")) +
  theme_classic()

# Print the formulas for each logistic function
cat("Logistic 1 Formula: Transmission Efficiency = ",
    fit1_params["L"], "/(1 + exp(-", fit1_params["k"], " * (x - ", fit1_params["x0"], ")))\n")
    
# Logistic 1 Formula: Transmission Efficiency =  0.5370646 /(1 + exp(- 1.351117  * (x -  5.275871 )))

cat("Modified Logistic 2 Formula: Transmission Efficiency = ",
    fit2_params["L"], " * (x / (1 + exp(-", fit2_params["k"], " * (x - ", fit2_params["x0"], "))))\n")
    
# Modified Logistic 2 Formula: Transmission Efficiency =  0.06764112  * (x / (1 + exp(- 1.964664  * (x -  4.591495 ))))

# Compare models using RSS
rss1 <- sum(residuals(fit1)^2)
rss2 <- sum(residuals(fit2)^2)

cat("Residual Sum of Squares (RSS):\n")
cat("Fit 1 (Logistic 1):", rss1, "\n")
cat("Fit 2 (Modified Logistic 2):", rss2, "\n")

# Residual Sum of Squares (RSS):
# Fit 1 (Logistic 1): 0.3144834
# Fit 2 (Modified Logistic 2): 0.3032254

# Compare models using AIC
aic1 <- AIC(fit1)
aic2 <- AIC(fit2)

cat("Akaike Information Criterion (AIC):\n")
cat("Fit 1 (Logistic 1):", aic1, "\n")
cat("Fit 2 (Modified Logistic 2):", aic2, "\n")

# Akaike Information Criterion (AIC):
#Fit 1 (Logistic 1): -13.76773
#Fit 2 (Modified Logistic 2): -14.42392




#****************************************************************************************
#****************************************************************************************
#****************************************************************************************
#******* STEP 2
#******* # Get the host viremic curves.

# Load the dataset
file_path <- "WNV_Host_Competence.xlsx"  # Replace with your file path
data <- read_excel(file_path)

# Ensure column names are syntactically valid
colnames(data) <- make.names(colnames(data))

# Filter data for curves with at least 3 data points and remove problematic values
filtered_data <- data %>%
  filter(!is.na(Day) & !is.na(Titer.PFU.) & is.finite(Day) & is.finite(Titer.PFU.)) %>%
  group_by(curve) %>%
  filter(n() >= 3) %>%
  ungroup()

# Define the exponential growth and decay model
exponential_growth_decay_model <- function(t, a1, b1, c1) {
  return(a1 * (t^b1) * exp(-c1 * t))
}

# Define the transmission efficiency formula
transmission_efficiency <- function(x) {
  0.06764112 * (x / (1 + exp(-1.964664 * (x - 4.591495))))
}

# Ensure the `species` column exists
if (!"species" %in% colnames(filtered_data)) {
  stop("The dataset must contain a 'species' column.")
}

# Get unique species
species_list <- unique(filtered_data$species)

# Prepare an empty list to store results
species_results_list <- list()

# Loop through each species
for (species_name in species_list) {
  
  # Filter data for the species
  species_data <- filtered_data %>%
    filter(species == species_name)
  
  if (nrow(species_data) == 0) next  # Skip if no data for the species
  
  # Aggregate the data by Day and compute the median of Titer.PFU
  aggregated_data <- species_data %>%
    group_by(Day) %>%
    summarize(Titer.PFU. = median(Titer.PFU., na.rm = TRUE), .groups = "drop")
  
  # Fit the exponential growth and decay model
  fit <- tryCatch({
    nlsLM(
      Titer.PFU. ~ exponential_growth_decay_model(Day, a1, b1, c1),
      data = aggregated_data,
      start = list(a1 = 1, b1 = 0.5, c1 = 0.1),
      lower = c(0, 0, 0),  # Constrain parameters to be positive
      upper = c(Inf, 20, 10),  # Logical upper bounds for parameters
      control = nls.lm.control(maxiter = 500)
    )
  }, error = function(e) {
    message("Error fitting model for species:", species_name, " Error:", e$message)
    return(NULL)
  })
  
  if (is.null(fit)) next
  
  # Extract model parameters
  params <- coef(fit)
  a1 <- params["a1"]
  b1 <- params["b1"]
  c1 <- params["c1"]
  
  # Generate predicted values for days
  days <- seq(0, max(aggregated_data$Day) + 5, by = 0.1)
  predicted_titers <- exponential_growth_decay_model(days, a1, b1, c1)
  predicted_titers[predicted_titers < 0] <- 0  # Clamp to 0
  
  # Calculate additional metrics
  peak_titer <- max(predicted_titers, na.rm = TRUE)
  peak_day <- days[which.max(predicted_titers)]
  duration_viremia <- sum(predicted_titers > 4) * 0.1  # Duration in days
  start_day <- min(days[predicted_titers > 4], na.rm = TRUE)
  end_day <- max(days[predicted_titers > 4], na.rm = TRUE)
  auc_viremia <- sum(predicted_titers[predicted_titers > 4] * 0.1, na.rm = TRUE)
  viral_load <- sum(predicted_titers[predicted_titers > 0] * 0.1, na.rm = TRUE)
  
  # Calculate Tmax
  Tmax <- tryCatch({
    Tmax_value <- (b1 + sqrt(b1)) / c1
    if (Tmax_value < 0 || Tmax_value > max(aggregated_data$Day) + 5) {
      NA  # Flag as NA if Tmax is outside the logical range
    } else {
      Tmax_value
    }
  }, error = function(e) {
    NA  # If calculation fails, return NA
  })
  
  # Calculate Vmax
  Vmax <- tryCatch({
    abs(
      a1 * sqrt(b1) * ((b1 + sqrt(b1)) / c1)^(b1 - 1) *
        exp(-((b1 + sqrt(b1)) / c1))
    )
  }, error = function(e) {
    NA  # If calculation fails, return NA
  })
  
  # Calculate transmission efficiency
  transmission_efficiency_values <- transmission_efficiency(predicted_titers)
  sum_transmission_efficiency <- sum(transmission_efficiency_values[predicted_titers > 4] * 0.1, na.rm = TRUE)
  
  # Annotate with number of curves used for species (assuming data has `curve` column)
  num_curves <- n_distinct(species_data$curve)
  
  # Save results with species information, host data, and number of curves used
  species_results_list[[species_name]] <- data.frame(
    species = species_name,
    a1 = a1,
    b1 = b1,
    c1 = c1,
    Tmax = Tmax,
    Vmax = Vmax,
    peak_titer = peak_titer,
    peak_day = peak_day,
    duration_viremia = duration_viremia,
    start_day = start_day,
    end_day = end_day,
    viral_load = viral_load,
    auc_viremia = auc_viremia,
    sum_transmission_efficiency = sum_transmission_efficiency,
    Host_Family = species_data$Host_Family[1],  # Assuming Host_Family is a column
    Host_Order = species_data$Host_Order[1],    # Assuming Host_Order is a column
    Class = species_data$Class[1],              # Assuming Class is a column
    num_viremia_curves = num_curves
  )
}

# Combine all results into a single data frame
all_species_results <- bind_rows(species_results_list)

# Save results to an Excel file
write_xlsx(all_species_results, "species_decay_results_with_metrics_and_host_info_viral_load.xlsx")

# Print the results
glimpse(all_species_results)



###-----------------------------------------------------------------------------
#******* STEP 2.1
##---Plot one species estimated curve

# Load the dataset
file_path <- "WNV_Host_Competence.xlsx"
data <- read_excel(file_path)

# Ensure column names are syntactically valid
colnames(data) <- make.names(colnames(data))

# Filter data for curves with at least 3 data points and remove problematic values
filtered_data <- data %>%
  filter(!is.na(Day) & !is.na(Titer.PFU.) & is.finite(Day) & is.finite(Titer.PFU.)) %>%
  group_by(curve) %>%
  filter(n() >= 3) %>%
  ungroup()

# Define the exponential growth and decay model
exponential_growth_decay_model <- function(t, a1, b1, c1) {
  return(a1 * (t^b1) * exp(-c1 * t))  # Exponential growth and decay model
}

# Ensure the `species` column exists and get the first species
if (!"species" %in% colnames(filtered_data)) {
  stop("The dataset must contain a 'species' column.")
}

# Select the species by name, for example "Pica pica"
species_name <- "Turdus migratorius"

# Filter data for the selected species
species_data <- filtered_data %>%
  filter(species == species_name)

# Aggregate the data by Day and compute the median of Titer.PFU
aggregated_data <- species_data %>%
  group_by(Day) %>%
  summarize(Titer.PFU. = median(Titer.PFU., na.rm = TRUE), .groups = "drop")

# Fit the exponential growth and decay model
fit_exponential_growth_decay <- tryCatch({
  nlsLM(
    Titer.PFU. ~ exponential_growth_decay_model(Day, a1, b1, c1),
    data = aggregated_data,
    start = list(a1 = 1, b1 = 0.5, c1 = 0.1),  # Starting parameters
    control = nls.lm.control(maxiter = 500)
  )
}, error = function(e) {
  stop("Error fitting model for species:", first_species_name, " Error:", e$message)
})

# Generate predicted values for plotting
days <- seq(0, 25, by = 0.1)  # Range of days from 0 to 25
predicted_titers <- predict(fit_exponential_growth_decay, newdata = data.frame(Day = days))

# Calculate the peak_day from the fitted model
params <- coef(fit_exponential_growth_decay)
peak_day <- days[which.max(predicted_titers)] 

# Create a data frame for plotting
plot_data <- data.frame(
  Day = days,
  Predicted_Titer = predicted_titers
)

# Plot the observed and predicted values with a vertical dashed line at peak_day
ggplot() +
  geom_point(data = aggregated_data, aes(x = Day, y = Titer.PFU.), color = "red4", size = 4) +  # Plot the original data points
  geom_line(data = plot_data, aes(x = Day, y = Predicted_Titer), color = "red4", size = 1.5) +  # Plot the fitted curve
  geom_vline(xintercept = peak_day, linetype = "dashed", color = "gray", linewidth = 0.8) +  # Add a vertical dashed line at peak_day
  geom_hline(yintercept = 4, linetype = "dashed", color = "gray", linewidth = 0.8) +
  annotate("text", x = peak_day + 0.5, y = max(predicted_titers) * 0.9, 
           label = paste("Peak Day:", round(peak_day, 2)), color = "white") +  # Annotate the peak day
  labs(
    title = paste("Estimated Survival Curve for", species_name),
    x = "Day",
    y = "Titer (PFU)"
  ) +
  theme_classic() +
  scale_x_continuous(limits = c(0, 10)) +  # Fix X-axis range from 0 to 25
  scale_y_continuous(limits = c(0, 10)) +  # Adjust Y-axis based on predicted values
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(size = 30),
    axis.text = element_text(size = 15)
  )


###-----------------------------------------------------------------------------
#******* STEP 2.2
#----- Plot curves at the family or order level

file_path <- "WNV_Host_Competence.xlsx"
data <- read_excel(file_path)

# Ensure column names are syntactically valid
colnames(data) <- make.names(colnames(data))

# Filter data for curves with at least 3 data points and remove problematic values
filtered_data <- data %>%
  filter(!is.na(Day) & !is.na(Titer.PFU.) & is.finite(Day) & is.finite(Titer.PFU.)) %>%
  group_by(curve) %>%
  filter(n() >= 3) %>%  # Keep curves with more than one observation
  ungroup()

# Define the exponential growth and decay model
exponential_growth_decay_model <- function(t, a1, b1, c1) {
  return(a1 * (t^b1) * exp(-c1 * t))  # Exponential growth and decay model
}

# Calculate the predicted viremia curves for each species
all_species_results_list <- list()

# Group data by Host_Family and loop through each family
for (family_name in unique(filtered_data$Host_Family)) {
  
  # Filter data for the current Host_Family
  family_data <- filtered_data %>%
    filter(Host_Family == family_name)
  
  # List to store results for each species within the family
  species_results_list <- list()
  
  # Loop through species within the family
  for (species_name in unique(family_data$species)) {
    
    # Filter data for the current species
    species_data <- family_data %>%
      filter(species == species_name)
    
    # Aggregate the data by Day and compute the median of Titer.PFU
    aggregated_data <- species_data %>%
      group_by(Day) %>%
      summarize(Titer.PFU. = median(Titer.PFU., na.rm = TRUE), .groups = "drop")
    
    # Fit the exponential growth and decay model
    fit_exponential_growth_decay <- tryCatch({
      nlsLM(
        Titer.PFU. ~ exponential_growth_decay_model(Day, a1, b1, c1),
        data = aggregated_data,
        start = list(a1 = 1, b1 = 0.5, c1 = 0.1),  # Starting parameters
        control = nls.lm.control(maxiter = 500)
      )
    }, error = function(e) {
      message("Error fitting model for species:", species_name, " Error:", e$message)
      return(NULL)
    })
    
    if (is.null(fit_exponential_growth_decay)) next  # Skip if fitting fails
    
    # Extract model parameters
    params <- coef(fit_exponential_growth_decay)
    a1 <- params["a1"]
    b1 <- params["b1"]
    c1 <- params["c1"]
    
    # Generate predicted values for plotting
    days <- seq(0, 25, by = 0.1)  # Range of days from 0 to 25
    predicted_titers <- predict(fit_exponential_growth_decay, newdata = data.frame(Day = days))
    
    # Store the results in the species_results_list
    species_results_list[[species_name]] <- data.frame(
      species = species_name,
      Host_Family = family_name,
      days = days,
      predicted_titers = predicted_titers
    )
  }
  
  # Combine all species results for the current family into one data frame
  family_results <- bind_rows(species_results_list)
  
  # Add the results to the all_species_results_list
  all_species_results_list[[family_name]] <- family_results
}

# Combine all species results across all families
all_species_results <- bind_rows(all_species_results_list)

# Calculate the mean viremia curve for each family
mean_curve_by_family <- all_species_results %>%
  group_by(Host_Family, days) %>%
  summarise(
    mean_titer = mean(predicted_titers, na.rm = TRUE),
    ci_lower = mean_titer - 1.96 * sd(predicted_titers) / sqrt(n()),  # 95% confidence interval
    ci_upper = mean_titer + 1.96 * sd(predicted_titers) / sqrt(n()),  # 95% confidence interval
    .groups = "drop"
  )

# Specify the families to plot (for example, you can choose any family you want)
selected_families <- c("Corvidae")  # Example families to plot

# Filter the data to include only the selected families
filtered_family_data <- all_species_results %>%
  filter(Host_Family %in% selected_families)

# Filter the mean curve data by the selected families
filtered_mean_curve_by_family <- mean_curve_by_family %>%
  filter(Host_Family %in% selected_families)

# Create the plot: species lines in light grey and mean curve by family with shaded CI
p1 <- ggplot() +
  # Plot individual species viremia curves in light grey
  geom_line(data = filtered_family_data, aes(x = days, y = predicted_titers, color = species), size = 0.5, alpha = 0.3) +  
  # Plot the mean viremia curve for each family
  geom_line(data = filtered_mean_curve_by_family, aes(x = days, y = mean_titer, color = Host_Family), size = 1) +
  # Add shaded confidence intervals for the mean viremia curve
  geom_ribbon(data = filtered_mean_curve_by_family, aes(x = days, ymin = ci_lower, ymax = ci_upper, fill = Host_Family), 
              alpha = 0.2) +  
  # Add a dashed line at Y = 4 (light grey)
  geom_hline(yintercept = 4, linetype = "dashed", color = "lightgrey", size = 0.8) +
  labs(
    title = "Viremia Curves Corvidae",
    x = "Days",
    y = "Predicted Viral Titer (PFU)"
  ) +
  theme_classic() +
  theme(
    legend.position = "none"  # Remove the legend as it's redundant in a faceted plot
  ) +
  coord_cartesian(ylim = c(0, 15), xlim = c(0, 20))  # Restrict x-axis to 12 days

# Display the plot
print(p1)






#****************************************************************************************
#****************************************************************************************
#****************************************************************************************
#******* STEP 3
#******* # Calculate the daily survival ratio for each day

# Load the dataset
file_path <- "WNV_Host_Competence.xlsx"
data <- read_excel(file_path)

# Strip leading/trailing spaces from column names if necessary
colnames(data) <- trimws(colnames(data))

# Clean column names and make them valid for selection
colnames(data) <- make.names(colnames(data), unique = TRUE)

# Filter out rows with NA values in 'Sample_Size_Survival' and 'Sample_Size_Original' columns
filtered_data <- data %>%
  filter(!is.na(Sample_Size_Survival) & !is.na(Sample_Size_Original)) %>%
  group_by(curve) %>%
  filter(n() >= 2) %>%  # Keep curves with more than one observation
  ungroup()

# Calculate the daily survival ratio for each day, curve, and species
daily_survival_data <- filtered_data %>%
  mutate(daily_survival_ratio = Sample_Size_Survival / Sample_Size_Original) %>%
  dplyr::select(species, curve, Day, daily_survival_ratio)

# Join the unique rows with survival data
data_unique <- data %>%
  dplyr::select(curve, Sample_Size_Original, Host_Family, Host_Order, Class) %>%
  distinct()

# Join the survival data with host family, order, and class information
daily_survival_data <- daily_survival_data %>%
  left_join(data_unique, by = "curve")

# View the data after the join
head(daily_survival_data)

# Now calculate the weighted daily survival ratio by grouping by species, day, and including host info
weighted_daily_survival_data <- daily_survival_data %>%
  group_by(species, Day, Host_Family, Host_Order, Class) %>%
  summarise(weighted_daily_survival_ratio = sum(daily_survival_ratio * Sample_Size_Original, na.rm = TRUE) / 
              sum(Sample_Size_Original, na.rm = TRUE)) %>%  # Calculate the weighted mean
  ungroup()

# View the weighted survival data with host information
head(weighted_daily_survival_data)

# Save the results to a new Excel file with host information
write_xlsx(weighted_daily_survival_data, "daily_survival_data.xlsx.xlsx")

# Calculate mean survival and standard error at the species level
species_survival_summary <- weighted_daily_survival_data %>%
  group_by(species, Host_Family, Host_Order, Class) %>%
  summarise(
    mean_daily_survival = mean(weighted_daily_survival_ratio, na.rm = TRUE),
    se_daily_survival = sd(weighted_daily_survival_ratio, na.rm = TRUE) / sqrt(n())
  ) %>%
  ungroup()

# View the summary
head(species_survival_summary)

# Final survival per curve: take the last observation for each curve
final_curve_survival <- filtered_data %>%
  group_by(curve) %>%
  arrange(Day) %>%
  slice_tail(n = 1) %>%  # Get the last day
  mutate(final_survival_ratio = Sample_Size_Survival / Sample_Size_Original) %>%
  dplyr::select(species, curve, final_survival_ratio, Sample_Size_Original, Host_Family, Host_Order, Class) %>%
  ungroup()

# Calculate mean survival ratio and standard error by species
species_final_survival_summary <- final_curve_survival %>%
  group_by(species, Host_Family, Host_Order, Class) %>%
  summarise(
    mean_survival_ratio = mean(final_survival_ratio, na.rm = TRUE),
    se_survival_ratio = sd(final_survival_ratio, na.rm = TRUE) / sqrt(n()),
    n_curves = n()
  ) %>%
  ungroup()

# View the summary
head(species_final_survival_summary)

# Sum of individuals tested per species (across all curves)
individuals_tested_per_species <- final_curve_survival %>%
  group_by(species) %>%
  summarise(
    total_individuals_tested = sum(Sample_Size_Original, na.rm = TRUE),
    .groups = "drop"
  )

# Merge all species-level summaries
species_combined_summary <- species_final_survival_summary %>%
  left_join(species_survival_summary, by = c("species", "Host_Family", "Host_Order", "Class")) %>%
  left_join(individuals_tested_per_species, by = "species")

# Calculate total individuals tested across all species
total_row <- species_combined_summary %>%
  summarise(
    species = "TOTAL",
    Host_Family = NA,
    Host_Order = NA,
    Class = NA,
    mean_survival_ratio = NA,
    se_survival_ratio = NA,
    n_curves = sum(n_curves, na.rm = TRUE),
    mean_daily_survival = NA,
    se_daily_survival = NA,
    total_individuals_tested = sum(total_individuals_tested, na.rm = TRUE)
  )

# Combine with the species-level summary
final_output <- bind_rows(species_combined_summary, total_row)

# Save the combined result to an Excel file
write_xlsx(final_output, "species_survival_summary_combined.xlsx")


###-----------------------------------------------------------------------------
#******* STEP 3.1
#----- Plot Daily Survival Rate

# Filter data for the first 10 species (you can change the number as needed)
first_10_species <- unique(weighted_daily_survival_data$species)[1:16]  # Get the first 10 species
filtered_data <- weighted_daily_survival_data %>%
  filter(species %in% first_10_species)

# Create a multiplot with facet_wrap
ggplot(filtered_data, aes(x = Day, y = weighted_daily_survival_ratio)) +
  geom_line(aes(color = species)) +  # Plot survival ratio over time
  geom_point(aes(color = species)) +  # Add points to show individual observations
  facet_wrap(~species, scales = "free_y") +  # Create separate plots for each species
  labs(
    title = "Weighted Daily Survival Ratio Over Time",
    x = "Day",
    y = "Weighted Daily Survival Ratio"
  ) +
  coord_cartesian(ylim = c(0, 1)) +  # Ensure y-axis goes from 0 to 1 for all plots
  theme_classic() +
  theme(
    legend.position = "none"  # Remove the legend as it's redundant in a faceted plot
  )


# Calculate the weighted daily survival ratio and confidence intervals
weighted_daily_survival_data <- daily_survival_data %>%
  group_by(species, Day) %>%
  summarise(
    weighted_daily_survival_ratio = sum(daily_survival_ratio * Sample_Size_Original, na.rm = TRUE) / 
      sum(Sample_Size_Original, na.rm = TRUE),  # Calculate the weighted mean
    # Calculate the standard error (SE) of the weighted mean
    se = sqrt(sum((Sample_Size_Original^2 * (daily_survival_ratio - weighted_daily_survival_ratio)^2) / 
                    sum(Sample_Size_Original), na.rm = TRUE) / sum(Sample_Size_Original, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(
    # Calculate the confidence intervals
    ci_lower = weighted_daily_survival_ratio - 1.96 * se,  # 95% confidence interval lower bound
    ci_upper = weighted_daily_survival_ratio + 1.96 * se   # 95% confidence interval upper bound
  )

# View the first few rows of the data with confidence intervals
head(weighted_daily_survival_data)

# Filter data for the first 10 species (you can change the number as needed)
first_10_species <- unique(weighted_daily_survival_data$species)[1:16]  # Get the first 10 species
filtered_data <- weighted_daily_survival_data %>%
  filter(species %in% first_10_species)

# Create a plot with a survival line and shaded area for confidence interval
ggplot(filtered_data, aes(x = Day, y = weighted_daily_survival_ratio)) +
  geom_line(aes(color = species), size = 1) +  # Plot survival ratio over time
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper, fill = species), alpha = 0.2) +  # Add shadow for CI
  geom_point(aes(color = species), size = 2) +  # Add points to show individual observations
  facet_wrap(~species, scales = "free_y") +  # Create separate plots for each species
  labs(
    title = "Weighted Daily Survival Ratio with Confidence Intervals",
    x = "Day",
    y = "Weighted Daily Survival Ratio"
  ) +
  coord_cartesian(ylim = c(0, 1)) +  # Ensure y-axis goes from 0 to 1 for all plots
  theme_classic() +
  theme(
    legend.position = "none"  # Remove the legend as it's redundant in a faceted plot
  )


##----------------same but grouping species by family

# Calculate the weighted daily survival ratio and confidence intervals by family and species
weighted_daily_survival_data <- daily_survival_data %>%
  group_by(Host_Family, species, Day) %>%
  summarise(
    weighted_daily_survival_ratio = sum(daily_survival_ratio * Sample_Size_Original, na.rm = TRUE) / 
      sum(Sample_Size_Original, na.rm = TRUE),  # Calculate the weighted mean
    # Calculate the standard error (SE) of the weighted mean
    se = sqrt(sum((Sample_Size_Original^2 * (daily_survival_ratio - weighted_daily_survival_ratio)^2) / 
                    sum(Sample_Size_Original), na.rm = TRUE) / sum(Sample_Size_Original, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(
    # Calculate the confidence intervals
    ci_lower = weighted_daily_survival_ratio - 1.96 * se,  # 95% confidence interval lower bound
    ci_upper = weighted_daily_survival_ratio + 1.96 * se   # 95% confidence interval upper bound
  )

# View the first few rows of the data with confidence intervals
head(weighted_daily_survival_data)

# Filter data for the first 6 families (or any other set of families)
first_6_families <- unique(weighted_daily_survival_data$Host_Family)[1:40]  # Get the first 6 families
filtered_data <- weighted_daily_survival_data %>%
  filter(Host_Family %in% first_6_families)

# Create a plot with survival lines and shaded areas for confidence intervals
ggplot(filtered_data, aes(x = Day, y = weighted_daily_survival_ratio)) +
  geom_line(aes(color = species), size = 1) +  # Plot survival ratio over time
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper, fill = species), alpha = 0.2) +  # Add shadow for CI
  geom_point(aes(color = species), size = 2) +  # Add points to show individual observations
  facet_wrap(~Host_Family, scales = "free_y") +  # Create separate plots for each family
  labs(
    title = "Weighted Daily Survival Ratio with Confidence Intervals by Family",
    x = "Day",
    y = "Weighted Daily Survival Ratio"
  ) +
  coord_cartesian(ylim = c(0, 1)) +  # Ensure y-axis goes from 0 to 1 for all plots
  theme_classic() +
  theme(
    legend.position = "none"  # Remove the legend as it's redundant in a faceted plot
  )


###--------------- same for selected families

# Calculate the weighted daily survival ratio and confidence intervals by family and species
weighted_daily_survival_data <- daily_survival_data %>%
  group_by(Host_Family, species, Day) %>%
  summarise(
    weighted_daily_survival_ratio = sum(daily_survival_ratio * Sample_Size_Original, na.rm = TRUE) / 
      sum(Sample_Size_Original, na.rm = TRUE),  # Calculate the weighted mean
    # Calculate the standard error (SE) of the weighted mean
    se = sqrt(sum((Sample_Size_Original^2 * (daily_survival_ratio - weighted_daily_survival_ratio)^2) / 
                    sum(Sample_Size_Original), na.rm = TRUE) / sum(Sample_Size_Original, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(
    # Calculate the confidence intervals
    ci_lower = weighted_daily_survival_ratio - 1.96 * se,  # 95% confidence interval lower bound
    ci_upper = weighted_daily_survival_ratio + 1.96 * se   # 95% confidence interval upper bound
  )

# View the first few rows of the data with confidence intervals
head(weighted_daily_survival_data)

# Specify the two selected families (you can change these as needed)
selected_families <- c("Corvidae")  # Example families

# Filter data for the selected families
filtered_data <- weighted_daily_survival_data %>%
  filter(Host_Family %in% selected_families)

# Create a plot with survival lines and shaded areas for confidence intervals for selected families
ggplot(filtered_data, aes(x = Day, y = weighted_daily_survival_ratio)) +
  geom_line(aes(color = species), size = 1) +  # Plot survival ratio over time
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper, fill = species), alpha = 0.1) +  # Add shadow for CI
  geom_point(aes(color = species), size = 1.5) +  # Add points to show individual observations
  facet_wrap(~Host_Family, scales = "free_y") +  # Create separate plots for each selected family
  labs(
    title = "Weighted Daily Survival Ratio with Confidence Intervals for Selected Families",
    x = "Day",
    y = "Weighted Daily Survival Ratio"
  ) +
  coord_cartesian(ylim = c(0, 1), xlim = c(0, 20)) +  # Restrict x-axis to 12 days
  theme_classic() +
  theme(
    legend.position = "none"  # Remove the legend as it's redundant in a faceted plot
  )


#---------------same for selected species
# Get a list of unique species from the weighted_daily_survival_data
unique_species <- unique(weighted_daily_survival_data$species)

# Print the unique species
print(unique_species)

# Subset data for selected species "Falco rusticolus"
species_data <- weighted_daily_survival_data %>%
  filter(species == "Turdus migratorius")

# Create a plot with a survival line and shaded area for confidence interval
ggplot(species_data, aes(x = Day, y = weighted_daily_survival_ratio)) +
  geom_line(aes(color = species), size = 1) +  # Plot survival ratio over time
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper, fill = species), alpha = 0.2) +  # Add shadow for CI
  geom_point(aes(color = species), size = 4) +  # Add points to show individual observations
  facet_wrap(~species, scales = "free_y") +  # Create separate plots for each species
  labs(
    title = "Weighted Daily Survival Ratio with Confidence Intervals",
    x = "Day",
    y = "Weighted Daily Survival Ratio"
  ) +
  coord_cartesian(ylim = c(0, 1), xlim = c(0, 10)) + # Ensure y-axis goes from 0 to 1 for all plots
  theme_classic() +
  theme(
    legend.position = "none"  # Remove the legend as it's redundant in a faceted plot
  )





#****************************************************************************************
#****************************************************************************************
#****************************************************************************************
#******* STEP 4
#******* # Calculate the daily transmission efficiency once daily & considering mortality

# Load the dataset
file_path <- "WNV_Host_Competence.xlsx"
data <- read_excel(file_path)

# Ensure column names are syntactically valid
colnames(data) <- make.names(colnames(data))

# Filter data for curves with at least 3 data points and remove problematic values
filtered_data <- data %>%
  filter(!is.na(Day) & !is.na(Titer.PFU.) & is.finite(Day) & is.finite(Titer.PFU.)) %>%
  group_by(curve) %>%
  filter(n() >= 3) %>%
  ungroup()

# Define the exponential growth and decay model
exponential_growth_decay_model <- function(t, a1, b1, c1) {
  return(a1 * (t^b1) * exp(-c1 * t))
}

# Define the transmission efficiency formula
transmission_efficiency <- function(x) {
  0.06764112 * (x / (1 + exp(-1.964664 * (x - 4.591495))))
}

# Ensure the `species` column exists
if (!"species" %in% colnames(filtered_data)) {
  stop("The dataset must contain a 'species' column.")
}

# Get unique species
species_list <- unique(filtered_data$species)

# Prepare an empty list to store results
species_results_list <- list()

# Loop through each species
for (species_name in species_list) {
  
  # Filter data for the species
  species_data <- filtered_data %>%
    filter(species == species_name)
  
  if (nrow(species_data) == 0) next  # Skip if no data for the species
  
  # Aggregate the data by Day and compute the median of Titer.PFU
  aggregated_data <- species_data %>%
    group_by(Day) %>%
    summarize(Titer.PFU. = median(Titer.PFU., na.rm = TRUE), .groups = "drop")
  
  # Fit the exponential growth and decay model
  fit <- tryCatch({
    nlsLM(
      Titer.PFU. ~ exponential_growth_decay_model(Day, a1, b1, c1),
      data = aggregated_data,
      start = list(a1 = 1, b1 = 0.5, c1 = 0.1),
      lower = c(0, 0, 0),  # Constrain parameters to be positive
      upper = c(Inf, 20, 10),  # Logical upper bounds for parameters
      control = nls.lm.control(maxiter = 500)
    )
  }, error = function(e) {
    message("Error fitting model for species:", species_name, " Error:", e$message)
    return(NULL)
  })
  
  if (is.null(fit)) next
  
  # Extract model parameters
  params <- coef(fit)
  a1 <- params["a1"]
  b1 <- params["b1"]
  c1 <- params["c1"]
  
  # Generate predicted values for integer days only (no fractional days)
  days <- seq(0, max(aggregated_data$Day), by = 1)  # Integer days only
  predicted_titers <- exponential_growth_decay_model(days, a1, b1, c1)
  predicted_titers[predicted_titers < 0] <- 0  # Clamp to 0
  
  # Calculate daily transmission efficiency values (only for each full day)
  transmission_efficiency_values <- transmission_efficiency(predicted_titers)
  
  # Create a binary column for predicted_titers > 4
  predicted_titers_flag <- ifelse(predicted_titers > 4, 1, 0)
  
  # Combine the results into a data frame for the species
  species_results_list[[species_name]] <- data.frame(
    species = species_name,
    Day = days,
    transmission_efficiency_values = transmission_efficiency_values,
    predicted_titers_flag = predicted_titers_flag
  )
}

# Combine all results into a single data frame
all_species_results <- bind_rows(species_results_list)

# Save the results to an Excel file
write_xlsx(all_species_results, "species_transmission_efficiency_results_daily.xlsx")

# Print the results
glimpse(all_species_results)


# Combine the daily transmission rate with the daily survival rate

# Load the transmission dataset
file_path <- "species_transmission_efficiency_results_daily.xlsx"  # Replace with your file path
transmission <- read_excel(file_path)

# Load the survival dataset
file_path1 <- "daily_survival_data.xlsx"  # Replace with your file path
survival <- read_excel(file_path1)

# Join the two datasets by species and Day
combined_data <- transmission %>%
  left_join(survival, by = c("species", "Day"))

# View the first few rows of the combined dataset
head(combined_data)

# Save the combined dataset to a new Excel file
write_xlsx(combined_data, "combined_transmission_survival_data.xlsx")

# Create the corrected transmission efficiency by multiplying the relevant columns
combined_data <- combined_data %>%
  mutate(corrected_transmission_efficiency = transmission_efficiency_values * 
                                                weighted_daily_survival_ratio * 
                                                predicted_titers_flag)

# View the first few rows of the updated dataset
head(combined_data)

# Save the updated dataset to a new Excel file
write_xlsx(combined_data, "corrected_transmission_efficiency_data.xlsx")

# Calculate the total transmission efficiency by species, 
# and also the sum of transmission_efficiency_values where titers > 4 (predicted_titers_flag = 1)
total_transmission_efficiency <- combined_data %>%
  group_by(species) %>%
  summarise(
    total_transmission_efficiency = sum(corrected_transmission_efficiency, na.rm = TRUE),  # Sum of corrected transmission efficiency
    total_transmission_efficiency_values_above_4 = sum(transmission_efficiency_values * (predicted_titers_flag == 1), na.rm = TRUE)  # Sum only for titers > 4
  ) %>%
  ungroup()

# View the result
head(total_transmission_efficiency)

# Save the result to a new Excel file
write_xlsx(total_transmission_efficiency, "corrected_transmission_efficiency_data.xlsx")



