#This R script analyzes host competence for West Nile Virus (WNV) by modeling viremia dynamics, vector transmission efficiency, and host survival, with the goal of estimating species-specific transmission potential. It fits mathematical models to experimental data, calculates summary metrics, and generates plots and Excel outputs for further analysis.

# Goal: Estimate species-level transmission potential by modeling:
#   (a) viremic dynamics per host species,
#   (b) vector transmission efficiency (Culex pipiens),
#   (c) host survival (daily → cumulative),
#   (d) mortality-corrected infectiousness,
# and generate summary tables + plotting

# STEP 1: Fit Logistic Curve for Vector Transmission Efficiency
#     Fit two logistic models using experimental Culex pipiens data.
#     Compare models using RSS and AIC.

# STEP 2: Data Load & Cleaning
#   - Read WNV_Host_Competence.xlsx
#   - Standardize column names; check required columns
#   - Keep finite Day and Titer.PFU. values
#   - Build:
#       agg_species: mean Titer.PFU. by species–day
#       taxa: species → Host_Family, Host_Order, Class

# Fit Viremic Curves (per species)
#   - Model: vgd(t) = a * t^b * exp(-c * t)  (growth–decay)
#   - Safe nlsLM fitting with adaptive starts + diagnostics
#   - Extract metrics on a fine grid:
#       Tmax = (b + sqrt(b))/c
#       Vmax = | -a * sqrt(b) * ((b + sqrt(b))/c)^(b-1) * exp(-(b + sqrt(b))) |
#       peak_titer, peak_day
#       start_day, end_day, duration_viremia (titer > threshold)
#       viral_load = total AUC; auc_viremia = AUC above threshold
#       sum_transmission_efficiency (uncorrected; above threshold)
#   - Save failures to model_fit_failures.xlsx
#   - species_metrics: metrics + taxonomy

# Survival (daily → weighted → cumulative)
#   - Build surv_daily from per-curve Sample_Size_Original/Survival
#   - Weighted daily survival per species–day (by Sample_Size_Original)
#   - cum_survival = cumprod(weighted_daily_survival_ratio)
#   - Summaries per species:
#       mean_daily_survival
#       mean_survival_ratio & se_survival_ratio (end-of-curve summary)

# Daily Transmission (Mortality-Corrected)
#   - daily_curves: merge predicted viremia (at observed days) + survival
#   - TE_raw_daily = TE(pred_viremia)
#   - Apply infectious mask (titer > threshold)
#   - TE_corrected_daily = TE_raw_daily * mask * cum_survival
#   - transmission_mortality_corrected = sum over days
#   - Clamp so corrected ≤ uncorrected (sanity guard)

# Counts & Sample Sizes
#   - per_species_counts:
#       n_curves = unique curves per species
#       individuals_tested = sum of initial Sample_Size_Original (first record per curve)

# STEP 3: Plotting Utilities
#   - Build viremia predictions from saved fits (fine grid)
#   - Viremia plots:
#       • by species: curve + observed points, shade area above threshold
#       • by family / by order: faint species curves + mean curve with 95% CI
#       • Fixed axes for species facets (x: 0–15, y: 0–10)
#   - Survival plots:
#       • daily weighted survival + 95% CI (optional)
#       • cumulative survival S(t) with 95% CI via delta method
#   - Combined species panels:
#       • Top: viremia (with shaded infectious area)
#       • Bottom: cumulative survival with CI


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

# ===========================
# Host Competence – Unified Pipeline (with diagnostics & extended metrics)
# ===========================
library(readxl)
library(writexl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(minpack.lm)
library(purrr)

# -----------------------------------------------------
# CONFIG
# -----------------------------------------------------
INPUT_XLSX <- "WNV_Host_Competence.xlsx"
OUT_DAILY  <- "daily_timeseries_with_survival_and_transmission.xlsx"
OUT_SPECIES <- "species_summary_parameters_and_indices.xlsx"
DIAG_SPECIES <- "Passer domesticus"
VIREMIA_THRESHOLD <- 4
FINE_DT <- 0.1

# -----------------------------------------------------
# MODEL DEFINITIONS
# -----------------------------------------------------
vgd <- function(t, a, b, c) a * (t^b) * exp(-c * t)

tmax_fun <- function(a, b, c){
  if (b > 0 & c > 0) (b + sqrt(b)) / c else NA_real_
}

vmax_fun <- function(a, b, c){
  if (b <= 0 | c <= 0) return(NA_real_)
  tstar <- (b + sqrt(b)) / c
  abs(-a * sqrt(b) * (tstar^(b - 1)) * exp(-(b + sqrt(b))))
}

# Culex pipiens transmission dose–response (your fitted curve)
transmission_efficiency <- function(x){
  0.06764112 * (x / (1 + exp(-1.964664 * (x - 4.591495))))
}

auc_sum <- function(y, dt) sum(y, na.rm = TRUE) * dt

# -----------------------------------------------------
# LOAD DATA
# -----------------------------------------------------
dat <- read_excel(INPUT_XLSX)
names(dat) <- make.names(names(dat))

required_cols <- c("species","curve","Day","Titer.PFU.","Sample_Size_Original",
                   "Sample_Size_Survival","Host_Family","Host_Order","Class")
missing <- setdiff(required_cols, names(dat))
if (length(missing) > 0) stop("Missing columns: ", paste(missing, collapse=", "))

dat <- dat %>% filter(is.finite(Day), is.finite(Titer.PFU.))

# -----------------------------------------------------
# AGGREGATE VIREMIA (mean by species-day)
# -----------------------------------------------------
agg_species <- dat %>%
  group_by(species, Day) %>%
  summarise(Titer.PFU. = mean(Titer.PFU., na.rm = TRUE), .groups="drop")

taxa <- dat %>%
  group_by(species) %>%
  summarise(across(c(Host_Family,Host_Order,Class), ~first(na.omit(.))), .groups="drop")

# per-species counts & sample sizes
per_species_counts <- dat %>%
  group_by(species) %>%
  summarise(
    n_curves = n_distinct(curve),
    individuals_tested = sum(Sample_Size_Original, na.rm = TRUE),
    .groups = "drop"
  )

# -----------------------------------------------------
# SURVIVAL: DAILY → WEIGHTED → CUMULATIVE
# -----------------------------------------------------
surv_daily <- dat %>%
  filter(!is.na(Sample_Size_Original), Sample_Size_Original>0) %>%
  mutate(daily_survival_ratio = Sample_Size_Survival / Sample_Size_Original) %>%
  dplyr::select(species, curve, Day, daily_survival_ratio, Sample_Size_Original)

surv_weighted <- surv_daily %>%
  group_by(species, Day) %>%
  summarise(
    weighted_daily_survival_ratio =
      sum(daily_survival_ratio * Sample_Size_Original, na.rm = TRUE) /
      sum(Sample_Size_Original, na.rm = TRUE),
    .groups="drop"
  ) %>%
  arrange(species, Day) %>%
  group_by(species) %>%
  mutate(cum_survival = cumprod(pmax(pmin(weighted_daily_survival_ratio,1),0))) %>%
  ungroup()

# species-level survival summaries
survival_species_summary <- surv_weighted %>%
  group_by(species) %>%
  summarise(
    mean_daily_survival = mean(weighted_daily_survival_ratio, na.rm = TRUE),
    mean_survival_ratio = mean(weighted_daily_survival_ratio, na.rm = TRUE),  # same metric
    se_survival_ratio   = sd(weighted_daily_survival_ratio, na.rm = TRUE) / sqrt(sum(!is.na(weighted_daily_survival_ratio))),
    .groups = "drop"
  )

# -----------------------------------------------------
# FIT VIREMIA MODEL PER SPECIES (with diagnostics)
# -----------------------------------------------------
fit_species <- function(df){
  if(nrow(df) < 3){
    return(list(model=NULL, reason="Too few data points (<3 days)"))
  }
  if(sd(df$Titer.PFU., na.rm=TRUE) < 0.1){
    return(list(model=NULL, reason="Viremia shows no variation"))
  }
  
  max_titer <- max(df$Titer.PFU., na.rm=TRUE)
  tmax_day <- max(df$Day, na.rm=TRUE)
  start_list <- list(
    a = max_titer,
    b = 1,
    c = ifelse(tmax_day > 0, 1/tmax_day, 0.1)
  )
  
  fit_attempt <- tryCatch(
    nlsLM(Titer.PFU. ~ vgd(Day,a,b,c),
          data=df,
          start=start_list,
          lower=c(0,0,0),
          control=nls.lm.control(maxiter=500)),
    error=function(e) e
  )
  if(inherits(fit_attempt, "error")){
    return(list(model=NULL, reason=paste("nlsLM failed:", fit_attempt$message)))
  }
  list(model=fit_attempt, reason="OK")
}

fits <- agg_species %>%
  group_by(species) %>%
  group_map(~{
    fit_info <- fit_species(.x)
    tibble(
      species   = .y$species,
      model     = list(fit_info$model),
      data      = list(.x),
      fit_status= fit_info$reason
    )
  }) %>% bind_rows()

fits_success <- fits %>% filter(!sapply(model, is.null))
fits_fail    <- fits %>% filter(sapply(model, is.null))

if(nrow(fits_fail)>0){
  write_xlsx(fits_fail, "model_fit_failures.xlsx")
  message("⚠️ Some species failed to fit. See: model_fit_failures.xlsx")
} else {
  message("✅ All species fitted successfully.")
}

# -----------------------------------------------------
# EXTRACT METRICS (fine grid, uncorrected transmission)
# -----------------------------------------------------
get_metrics <- function(model, df){
  pars <- coef(model)
  a <- pars["a"]; b <- pars["b"]; c <- pars["c"]
  
  # Fine grid for continuous metrics
  max_day <- max(df$Day, na.rm = TRUE)
  days_fine <- seq(0, max_day + 5, by = FINE_DT)
  pred_fine <- pmax(vgd(days_fine, a, b, c), 0)
  
  # Tmax & Vmax
  tmax <- tmax_fun(a, b, c)
  vmax <- vmax_fun(a, b, c)
  
  # Peak
  peak_val <- max(pred_fine, na.rm = TRUE)
  peak_day <- days_fine[which.max(pred_fine)]
  
  # Threshold region
  above_idx <- which(pred_fine > VIREMIA_THRESHOLD)
  start_day <- if (length(above_idx)) days_fine[min(above_idx)] else NA_real_
  end_day   <- if (length(above_idx)) days_fine[max(above_idx)] else NA_real_
  duration  <- if (length(above_idx)) length(above_idx) * FINE_DT else 0
  
  # AUCs
  viral_load  <- auc_sum(pred_fine, FINE_DT)   # TOTAL AUC
  auc_viremia <- auc_sum(pred_fine[above_idx], FINE_DT)  # AUC above threshold
  
  # Transmission potential (uncorrected)
  trans_vals <- transmission_efficiency(pred_fine)
  sum_transmission_efficiency <- auc_sum(trans_vals[above_idx], FINE_DT)
  
  tibble(
    a=a, b=b, c=c,
    Tmax=tmax, Vmax=vmax,
    peak_titer=peak_val, peak_day=peak_day,
    duration_viremia=duration,
    start_day=start_day, end_day=end_day,
    viral_load=viral_load,
    auc_viremia=auc_viremia,
    sum_transmission_efficiency=sum_transmission_efficiency
  )
}

species_metrics <- fits_success %>%
  mutate(params = map2(model, data, get_metrics)) %>%
  unnest(params) %>%
  left_join(taxa, by="species")


# -----------------------------------------------------
# DAILY SERIES + mortality-corrected transmission
# -----------------------------------------------------
daily_curves <- fits_success %>%
  mutate(curve_data = map2(model, data, ~{
    days <- sort(unique(.y$Day))
    pr <- pmax(vgd(days, coef(.x)[1], coef(.x)[2], coef(.x)[3]), 0)
    tibble(Day = days, pred_viremia = pr)
  })) %>%
  dplyr::select(species, curve_data) %>%
  unnest(curve_data) %>%
  left_join(surv_weighted, by=c("species","Day")) %>%
  mutate(
    TE_raw_daily = transmission_efficiency(pred_viremia),
    TE_mask = ifelse(pred_viremia > VIREMIA_THRESHOLD, 1, 0),
    TE_uncorrected_daily = TE_raw_daily * TE_mask,
    TE_corrected_daily   = TE_uncorrected_daily * pmin(pmax(cum_survival, 0), 1)
  )

mortality_corrected_TE <- daily_curves %>%
  group_by(species) %>%
  summarise(
    transmission_mortality_corrected = sum(TE_corrected_daily, na.rm = TRUE),
    .groups = "drop"
  )

# Merge corrected TE back into species metrics
species_metrics <- species_metrics %>%
  left_join(mortality_corrected_TE, by="species") %>%
  mutate(
    transmission_mortality_corrected =
      pmin(transmission_mortality_corrected, sum_transmission_efficiency + 1e-12)
  )


# -----------------------------------------------------
# PER-SPECIES COUNTS: n_curves & individuals_tested (initial N per curve)
# -----------------------------------------------------
# Safety: find the Sample_Size_Original column even if make.names changed it
orig_col <- names(dat)[grepl("^Sample[_\\.]?Size[_\\.]?Original$", names(dat), ignore.case = TRUE)]
if (length(orig_col) == 0) stop("Column 'Sample_Size_Original' not found after make.names().")

# One row per curve with INITIAL sample size
curve_initial_n <- dat %>%
  arrange(curve, Day) %>%
  group_by(species, curve) %>%
  slice_head(n = 1) %>%                    # first record per curve = initial N
  ungroup() %>%
  dplyr::select(species, curve, Sample_Size_Original = all_of(orig_col))

per_species_counts <- curve_initial_n %>%
  group_by(species) %>%
  summarise(
    n_curves = n_distinct(curve),
    individuals_tested = sum(Sample_Size_Original, na.rm = TRUE),
    .groups = "drop"
  )

# -----------------------------------------------------
# SURVIVAL SUMMARIES:
#   - mean_daily_survival (from weighted daily survival over time)
#   - mean_survival_ratio & se_survival_ratio (final survival per curve)
# -----------------------------------------------------

# 1) mean_daily_survival from daily weighted data
survival_daily_summary <- surv_weighted %>%
  group_by(species) %>%
  summarise(
    mean_daily_survival = mean(weighted_daily_survival_ratio, na.rm = TRUE),
    .groups = "drop"
  )

# 2) Final survival per curve -> mean & SE by species
final_curve_survival <- dat %>%
  filter(!is.na(Sample_Size_Survival), !is.na(!!sym(orig_col))) %>%
  arrange(curve, Day) %>%
  group_by(species, curve) %>%
  slice_tail(n = 1) %>%  # last day observed for that curve
  ungroup() %>%
  mutate(
    final_survival_ratio = Sample_Size_Survival / !!sym(orig_col)
  ) %>%
  dplyr::select(species, curve, final_survival_ratio)

survival_curve_summary <- final_curve_survival %>%
  group_by(species) %>%
  summarise(
    mean_survival_ratio = mean(final_survival_ratio, na.rm = TRUE),
    se_survival_ratio   = sd(final_survival_ratio,  na.rm = TRUE) /
      sqrt(sum(!is.na(final_survival_ratio))),
    .groups = "drop"
  )

# -----------------------------------------------------
# MERGE EVERYTHING INTO species_metrics 
# -----------------------------------------------------
species_metrics_final <- species_metrics %>%
  # add counts
  left_join(per_species_counts,     by = "species") %>%
  # add survival summaries
  left_join(survival_daily_summary, by = "species") %>%
  left_join(survival_curve_summary, by = "species") %>%
  # clamp (column already present)
  mutate(
    transmission_mortality_corrected =
      pmin(transmission_mortality_corrected, sum_transmission_efficiency + 1e-12)
  ) %>%
  # tidy column order
  relocate(Host_Family, Host_Order, Class, .after = species) %>%
  relocate(sum_transmission_efficiency, transmission_mortality_corrected,
           .after = auc_viremia) %>%
  relocate(n_curves, individuals_tested, .after = transmission_mortality_corrected) %>%
  relocate(mean_daily_survival, mean_survival_ratio, se_survival_ratio,
           .after = individuals_tested)


# Save outputs
write_xlsx(daily_curves, OUT_DAILY)
write_xlsx(species_metrics_final, OUT_SPECIES)


# -----------------------------------------------------
# DIAGNOSTIC PLOT (uses fitted model from fits_success)
# -----------------------------------------------------
if (DIAG_SPECIES %in% unique(fits_success$species)) {
  diag_row <- fits_success %>% filter(species == DIAG_SPECIES) %>% slice(1)
  sp_v <- agg_species %>% filter(species == DIAG_SPECIES) %>% arrange(Day)
  m <- diag_row$model[[1]]
  pars <- coef(m); a <- unname(pars["a"]); b <- unname(pars["b"]); c <- unname(pars["c"])
  days_fine <- seq(0, max(ceiling(max(sp_v$Day)), 25), by = FINE_DT)
  pred_fine <- vgd(days_fine, a, b, c); pred_fine[pred_fine < 0] <- 0
  tmax <- tmax_fun(a,b,c)
  
  p <- ggplot(sp_v, aes(Day, Titer.PFU.)) +
    geom_point(size = 3, color = "black") +
    geom_line(data = tibble(Day = days_fine, Pred = pred_fine),
              aes(Day, Pred), linewidth = 1.2, color = "#1b9e77") +
    geom_hline(yintercept = VIREMIA_THRESHOLD, linetype = "dashed", color = "grey50") +
    { if (!is.na(tmax)) geom_vline(xintercept = tmax, linetype = "dotted", color = "grey50") } +
    labs(title = paste0("Viremia fit: ", DIAG_SPECIES),
         x = "Day", y = "Titer (log10 PFU)") +
    theme_classic(base_size = 14)
  
  print(p)
} else {
  message("Diagnostic plot: species '", DIAG_SPECIES, "' not found among successful fits.")
}

# -----------------------------------------------------
# SAVE OUTPUTS
# -----------------------------------------------------
write_xlsx(daily_curves, OUT_DAILY)
write_xlsx(species_metrics, OUT_SPECIES)




#****************************************************************************************
#****************************************************************************************
#****************************************************************************************
#**** STEP 3
## -------- PLOTS --------------------------------------------------------------------------
##---------------------------------------------------------------------------------------###

# ==========================================
# VIREMIA PREDICTIONS PLOTTING
# ==========================================

# Build fine-grid predictions per species from saved fits
build_viremia_predictions <- function(fits_tbl, taxa, dt = FINE_DT, extend_days = 5) {
  library(dplyr)
  library(purrr)
  library(tidyr)
  vgd <- function(t, a, b, c) a * (t^b) * exp(-c * t)
  
  preds <- fits_tbl %>%
    mutate(curve_data = map2(model, data, ~{
      pars <- coef(.x); a <- unname(pars["a"]); b <- unname(pars["b"]); c <- unname(pars["c"])
      max_day <- max(.y$Day, na.rm = TRUE)
      days <- seq(0, max_day + extend_days, by = dt)
      tibble(Day = days, Pred = pmax(vgd(days, a, b, c), 0))
    })) %>%
    dplyr::select(species, curve_data) %>%
    unnest(curve_data) %>%
    left_join(taxa, by = "species")
  
  preds
}

viremia_preds <- build_viremia_predictions(fits_success, taxa)

# -------------------------
# Plot by SPECIES (faceted)
# -------------------------
plot_viremia_by_species <- function(species_names,
                                    preds = viremia_preds,
                                    raw = agg_species,
                                    threshold = VIREMIA_THRESHOLD) {
  
  library(dplyr)
  library(ggplot2)
  
  df_pred <- preds %>% filter(species %in% species_names)
  df_raw  <- raw  %>% filter(species %in% species_names)
  
  # Create shading region only where Pred > threshold
  df_shade <- df_pred %>%
    mutate(above = Pred > threshold) %>%
    group_by(species) %>%
    mutate(group_id = cumsum(!above)) %>%   # separate ribbons
    filter(above)
  
  ggplot() +
    # Shaded infectious region
    geom_ribbon(data = df_shade,
                aes(x = Day, ymin = threshold, ymax = Pred, fill = species),
                alpha = 0.25, color = NA) +
    # Fitted viremia curve
    geom_line(data = df_pred,
              aes(Day, Pred, color = species),
              linewidth = 1.2) +
    # Observed data points
    geom_point(data = df_raw,
               aes(Day, Titer.PFU., color = species),
               size = 2, alpha = 0.9) +
    # Horizontal line at threshold
    geom_hline(yintercept = threshold, linetype = "dashed", color = "grey50") +
    facet_wrap(~ species, scales = "fixed") +
    labs(title = "Observed and Predicted Viremia Curves (with infectious area shaded)",
         x = "Day",
         y = "Viremia (log10 PFU/mL)") +
    coord_cartesian(xlim = c(0, 15), ylim = c(0, 10)) +   # << FIXED AXES
    theme_classic(base_size = 14) +
    theme(legend.position = "none")
}


# ------------------------
# Plot by FAMILY (faceted)
#   - light species lines
#   - family mean curve + 95% CI
# ------------------------
plot_viremia_by_family <- function(families, data = viremia_preds, xlim = NULL, ylim = c(0, NA)) {
  library(ggplot2)
  df <- data %>% dplyr::filter(Host_Family %in% families)
  
  fam_mean <- df %>%
    dplyr::group_by(Host_Family, Day) %>%
    dplyr::summarise(
      mean_titer = mean(Pred, na.rm = TRUE),
      se = sd(Pred, na.rm = TRUE) / sqrt(sum(!is.na(Pred))),
      ci_lower = mean_titer - 1.96 * se,
      ci_upper = mean_titer + 1.96 * se,
      .groups = "drop"
    )
  
  ggplot() +
    # species-level curves (faint)
    geom_line(data = df, aes(Day, Pred, group = species, color = species),
              linewidth = 0.6, alpha = 0.35) +
    # family mean + CI
    geom_ribbon(data = fam_mean, aes(Day, ymin = ci_lower, ymax = ci_upper, fill = Host_Family), alpha = 0.20) +
    geom_line(data = fam_mean, aes(Day, mean_titer, color = Host_Family), linewidth = 1.1) +
    geom_hline(yintercept = VIREMIA_THRESHOLD, linetype = "dashed", color = "grey60") +
    facet_wrap(~ Host_Family, scales = "free_y") +
    labs(title = "Predicted viremia curves (by family)",
         x = "Day", y = "Predicted titer (log10 PFU)") +
    coord_cartesian(xlim = xlim, ylim = ylim) +
    theme_classic() +
    theme(legend.position = "none")
}

# -----------------------
# Plot by ORDER (faceted)
# -----------------------
plot_viremia_by_order <- function(orders, data = viremia_preds, xlim = NULL, ylim = c(0, NA)) {
  library(ggplot2)
  df <- data %>% dplyr::filter(Host_Order %in% orders)
  
  ord_mean <- df %>%
    dplyr::group_by(Host_Order, Day) %>%
    dplyr::summarise(
      mean_titer = mean(Pred, na.rm = TRUE),
      se = sd(Pred, na.rm = TRUE) / sqrt(sum(!is.na(Pred))),
      ci_lower = mean_titer - 1.96 * se,
      ci_upper = mean_titer + 1.96 * se,
      .groups = "drop"
    )
  
  ggplot() +
    geom_line(data = df, aes(Day, Pred, group = species, color = species),
              linewidth = 0.6, alpha = 0.35) +
    geom_ribbon(data = ord_mean, aes(Day, ymin = ci_lower, ymax = ci_upper, fill = Host_Order), alpha = 0.20) +
    geom_line(data = ord_mean, aes(Day, mean_titer, color = Host_Order), linewidth = 1.1) +
    geom_hline(yintercept = VIREMIA_THRESHOLD, linetype = "dashed", color = "grey60") +
    facet_wrap(~ Host_Order, scales = "free_y") +
    labs(title = "Predicted viremia curves (by order)",
         x = "Day", y = "Predicted titer (log10 PFU)") +
    coord_cartesian(xlim = xlim, ylim = ylim) +
    theme_classic() +
    theme(legend.position = "none")
}

# -----------------
# examples
# -----------------
plot_viremia_by_species(c("Passer domesticus", "Turdus migratorius"))
plot_viremia_by_family(c("Corvidae"))
plot_viremia_by_order(c("Passeriformes"))



# ============================
# CUMULATIVE SURVIVAL PLOTTING
# ============================

# Build cumulative survival + CI (delta method on log scale)
# Inputs:
#   surv_daily: columns species, curve, Day, daily_survival_ratio, Sample_Size_Original
#   taxa:       species → Host_Family, Host_Order, Class
# Output:
#   one row per species–day with:
#     mu_daily   (weighted daily survival ratio)
#     cum_survival, cum_lower, cum_upper  (monotone)
build_cum_survival_ci <- function(surv_daily, taxa) {
  eps <- 1e-8  # numerical guard
  
  surv_with_taxa <- surv_daily %>%
    left_join(taxa, by = "species")
  
  daily <- surv_with_taxa %>%
    group_by(species, Host_Family, Host_Order, Class, Day) %>%
    summarise(
      w = sum(Sample_Size_Original, na.rm = TRUE),
      mu_daily = sum(daily_survival_ratio * Sample_Size_Original, na.rm = TRUE) /
        pmax(sum(Sample_Size_Original, na.rm = TRUE), eps),
      .groups = "drop"
    ) %>%
    mutate(
      # clamp to (0,1] to avoid log/CI issues
      mu_daily = pmin(pmax(mu_daily, eps), 1)
    ) %>%
    arrange(species, Day)
  
  # Binomial-ish SE for a weighted mean ratio; conservative
  daily <- daily %>%
    mutate(
      se_daily  = sqrt( (mu_daily * (1 - mu_daily)) / pmax(w, 1) ),
      var_logmu = (se_daily^2) / pmax(mu_daily^2, eps)  # delta method: Var(log r)
    )
  
  # Cumulative product and delta-method CI on log scale
  cum <- daily %>%
    group_by(species, Host_Family, Host_Order, Class) %>%
    arrange(Day, .by_group = TRUE) %>%
    mutate(
      logS      = cumsum(log(mu_daily)),
      var_logS  = cumsum(var_logmu),
      cum_survival = exp(logS),
      cum_lower    = exp(logS - 1.96 * sqrt(var_logS)),
      cum_upper    = exp(logS + 1.96 * sqrt(var_logS))
    ) %>%
    ungroup() %>%
    mutate(
      # cumulative survival must be in (0,1]; enforce numerically
      cum_survival = pmin(pmax(cum_survival, eps), 1),
      cum_lower    = pmin(pmax(cum_lower,    eps), 1),
      cum_upper    = pmin(pmax(cum_upper,    eps), 1)
    )
  
  # Add Day 0 = 1 with zero variance for each species (nice for plots)
  day0 <- cum %>%
    group_by(species, Host_Family, Host_Order, Class) %>%
    summarise(Day = 0,
              mu_daily = NA_real_, se_daily = NA_real_,
              cum_survival = 1, cum_lower = 1, cum_upper = 1,
              .groups = "drop")
  
  bind_rows(day0, cum) %>%
    arrange(species, Day)
}

# Build once; reuse for all plots
cum_survival_ci <- build_cum_survival_ci(surv_daily, taxa)

# ----------------------------
# Plotters (CUMULATIVE curves)
# ----------------------------

plot_cum_survival_by_species <- function(species_names, data = cum_survival_ci, xlim = NULL) {
  df <- data %>% filter(species %in% species_names)
  ggplot(df, aes(Day, cum_survival)) +
    geom_ribbon(aes(ymin = cum_lower, ymax = cum_upper, fill = species), alpha = 0.18) +
    geom_line(aes(color = species), linewidth = 1.05) +
    geom_point(aes(color = species), size = 1.6) +
    facet_wrap(~ species, scales = "free_y") +
    labs(
      title = "Cumulative survival with 95% CI (by species)",
      x = "Day", y = "Cumulative survival"
    ) +
    coord_cartesian(ylim = c(0, 1), xlim = xlim) +
    theme_classic() +
    theme(legend.position = "none")
}

plot_cum_survival_by_family <- function(families, data = cum_survival_ci, xlim = NULL) {
  df <- data %>% filter(Host_Family %in% families)
  ggplot(df, aes(Day, cum_survival)) +
    geom_ribbon(aes(ymin = cum_lower, ymax = cum_upper, fill = species), alpha = 0.12) +
    geom_line(aes(color = species), linewidth = 0.9) +
    facet_wrap(~ Host_Family, scales = "free_y") +
    labs(
      title = "Cumulative survival with 95% CI (by family)",
      x = "Day", y = "Cumulative survival"
    ) +
    coord_cartesian(ylim = c(0, 1), xlim = xlim) +
    theme_classic() +
    theme(legend.position = "none")
}

plot_cum_survival_by_order <- function(orders, data = cum_survival_ci, xlim = NULL) {
  df <- data %>% filter(Host_Order %in% orders)
  ggplot(df, aes(Day, cum_survival)) +
    geom_ribbon(aes(ymin = cum_lower, ymax = cum_upper, fill = species), alpha = 0.12) +
    geom_line(aes(color = species), linewidth = 0.9) +
    facet_wrap(~ Host_Order, scales = "free_y") +
    labs(
      title = "Cumulative survival with 95% CI (by order)",
      x = "Day", y = "Cumulative survival"
    ) +
    coord_cartesian(ylim = c(0, 1), xlim = xlim) +
    theme_classic() +
    theme(legend.position = "none")
}

# Quick helper
plot_first_n_species_cum <- function(n = 12, data = cum_survival_ci, xlim = NULL) {
  first_n <- data %>% distinct(species) %>% arrange(species) %>% slice(1:n) %>% pull(species)
  plot_cum_survival_by_species(first_n, data, xlim)
}

# -----------------
# examples
# -----------------
plot_cum_survival_by_species(c("Passer domesticus", "Turdus migratorius"))
plot_first_n_species_cum(16)
plot_cum_survival_by_family(c("Corvidae"))
plot_cum_survival_by_order(c("Passeriformes"))



# ============================
# Viremia + Cumulative Survival for selected species
# ============================

# Reuse objects: surv_daily (species, curve, Day, daily_survival_ratio, Sample_Size_Original)
# and taxa (species → Host_Family, Host_Order, Class)

build_survival_ci <- function(surv_daily, taxa) {
  surv_daily %>%
    dplyr::left_join(taxa, by = "species") %>%
    dplyr::group_by(species, Host_Family, Host_Order, Class, Day) %>%
    dplyr::summarise(
      w = sum(Sample_Size_Original, na.rm = TRUE),
      mu = sum(daily_survival_ratio * Sample_Size_Original, na.rm = TRUE) / pmax(w, 1e-12),
      se = {
        num <- sum((Sample_Size_Original^2) * (daily_survival_ratio - mu)^2, na.rm = TRUE)
        den <- (pmax(w, 1e-12))^2
        sqrt(num / pmax(den, 1e-12))
      },
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      weighted_daily_survival_ratio = pmin(pmax(mu, 0), 1),
      ci_lower = pmax(0, mu - 1.96 * se),
      ci_upper = pmin(1, mu + 1.96 * se)
    ) %>%
    dplyr::select(-mu, -w)
}

survival_ci <- build_survival_ci(surv_daily, taxa)

build_cum_survival_ci <- function(survival_ci_tbl) {
  eps <- 1e-8
  survival_ci_tbl %>%
    dplyr::arrange(species, Day) %>%
    dplyr::group_by(species) %>%
    dplyr::mutate(
      r  = pmin(pmax(weighted_daily_survival_ratio, eps), 1),  # clamp
      lr = log(r),
      var_log_r = (ifelse(se > 0, se, 0) / pmax(r, eps))^2,
      logS = cumsum(lr),
      var_logS = cumsum(var_log_r),
      cum_survival = exp(logS),
      cum_low  = exp(logS - 1.96 * sqrt(var_logS)),
      cum_high = exp(logS + 1.96 * sqrt(var_logS))
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      cum_survival = pmin(pmax(cum_survival, 0), 1),
      cum_low  = pmin(pmax(cum_low, 0), 1),
      cum_high = pmin(pmax(cum_high, 0), 1)
    )
}

cum_survival_ci <- build_cum_survival_ci(survival_ci)

library(patchwork)

plot_viremia_and_cum_survival <- function(species_names,
                                          preds = viremia_preds,   # from previous step
                                          raw   = agg_species,     # observed viremia points
                                          survival_ci_tbl = cum_survival_ci,   # new CI table
                                          threshold = VIREMIA_THRESHOLD) {
  # filter
  df_pred <- preds %>% dplyr::filter(species %in% species_names)
  df_raw  <- raw   %>% dplyr::filter(species %in% species_names)
  df_surv <- survival_ci_tbl %>% dplyr::filter(species %in% species_names)
  
  # palette
  sp_levels <- unique(df_pred$species)
  pal <- setNames(scales::hue_pal()(length(sp_levels)), sp_levels)
  
  # shading for infectious area
  df_shade <- df_pred %>%
    dplyr::mutate(above = Pred > threshold) %>%
    dplyr::group_by(species) %>%
    dplyr::mutate(group_id = cumsum(!above)) %>%
    dplyr::filter(above)
  
  # top: viremia
  p_v <- ggplot() +
    geom_ribbon(data = df_shade,
                aes(x = Day, ymin = threshold, ymax = Pred,
                    fill = species, group = interaction(species, group_id)),
                alpha = 0.25, color = NA) +
    geom_line(data = df_pred, aes(Day, Pred, color = species), linewidth = 1.2) +
    geom_point(data = df_raw, aes(Day, Titer.PFU., color = species), size = 2) +
    geom_hline(yintercept = threshold, linetype = "dashed", color = "grey50") +
    facet_wrap(~ species, scales = "fixed") +
    coord_cartesian(xlim = c(0, 15), ylim = c(0, 10)) +
    scale_color_manual(values = pal) +
    scale_fill_manual(values  = pal) +
    labs(title = "Observed & predicted viremia (infectious area shaded)",
         x = "Day", y = "Viremia (log10 PFU/mL)") +
    theme_classic(base_size = 13) +
    theme(legend.position = "none", strip.text = element_text(face = "bold"))
  
  # bottom: cumulative survival with CI
  p_s <- ggplot(df_surv, aes(Day, cum_survival, color = species, fill = species)) +
    geom_ribbon(aes(ymin = cum_low, ymax = cum_high), alpha = 0.18, color = NA) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 1.6) +
    facet_wrap(~ species, scales = "fixed") +
    coord_cartesian(xlim = c(0, 15), ylim = c(0, 1)) +
    scale_color_manual(values = pal) +
    scale_fill_manual(values  = pal) +
    labs(title = "Cumulative survival with 95% CI",
         x = "Day", y = "S(t)") +
    theme_classic(base_size = 13) +
    theme(legend.position = "none", strip.text = element_text(face = "bold"))
  
  p_v / p_s + plot_layout(heights = c(2, 1))
}

# Examples
plot_viremia_and_cum_survival(c("Passer domesticus"))
plot_viremia_and_cum_survival(c("Passer domesticus","Turdus migratorius"))



