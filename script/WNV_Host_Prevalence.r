# ======================================================================================
# WNV HOST PREVALENCE ANALYSIS
# ======================================================================================
#
# Description:
# This script analyzes the WNV host prevalence database using the revised long-format
# structure, where a single biological observation may correspond to multiple assay rows.
# To avoid double-counting the same sampled record, most descriptive summaries are based
# on an observation-level dataset collapsed by observation_code.
#
# Main data objects:
#   - wnv_long: row-level diagnostic data (one row per assay result)
#   - wnv_obs:  one row per observation_code (used for most descriptive summaries)
#
# General rule:
#   - Use wnv_obs for most counts, summaries, site-level analyses, country analyses,
#     temporal trends, and species tested/positive summaries.
#   - Use wnv_long for method-specific analyses (diagnostic method, assay stage,
#     screening vs confirmatory, etc.).
# ======================================================================================

# ---------------------------------
# LOAD PACKAGES
# ---------------------------------
library(readxl)
library(writexl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(patchwork)
library(RColorBrewer)
library(janitor)
library(scales)
library(stringr)

# ---------------------------------
# CONFIGURATION
# ---------------------------------
INPUT_XLSX <- "WNV_Host_Prevalence.xlsx"
MAP_REALMS <- "maps/newRealms.shp"

# Optional: references without usable prevalence information
exclude_ref_ids <- c(304, 222, 323)

# ---------------------------------
# HELPER FUNCTIONS
# ---------------------------------
first_non_na <- function(x) {
  x <- x[!is.na(x) & x != ""]
  if (length(x) == 0) NA else x[1]
}

safe_max_num <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (length(x) == 0) NA_real_ else max(x, na.rm = TRUE)
}

safe_any_positive <- function(pos, ratio) {
  pos <- suppressWarnings(as.numeric(pos))
  ratio <- suppressWarnings(as.numeric(ratio))
  any((!is.na(pos) & pos > 0) | (!is.na(ratio) & ratio > 0), na.rm = TRUE)
}

build_species_status <- function(data, group_vars) {
  data %>%
    group_by(across(all_of(group_vars)), accepted_scientific_name) %>%
    summarise(
      positive = any(positive_observation, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(status = if_else(positive, "Positive", "Negative"))
}

plot_site_map <- function(sf_data, map_data, metric, color_hex, title_text,
                          subtitle_text, size_range = c(1.5, 8),
                          legend_name = NULL, legend_breaks = waiver()) {
  
  ggplot() +
    geom_sf(data = map_data, fill = "grey90", color = "black") +
    geom_sf(
      data = sf_data,
      aes(size = .data[[metric]]),
      color = color_hex,
      fill = NA,
      linewidth = 1.2,
      alpha = 0.35
    ) +
    scale_size_continuous(
      name = legend_name %||% metric,
      range = size_range,
      breaks = legend_breaks,
      labels = comma
    ) +
    labs(
      title = title_text,
      subtitle = subtitle_text
    ) +
    theme_classic()
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

# ======================================================================================
# STEP 1. LOAD AND PREPARE THE PREVALENCE DATABASE
# ======================================================================================

# Read workbook
wnv_long <- read_excel(INPUT_XLSX) %>%
  janitor::clean_names()

# Check required columns
required_cols <- c(
  "observation_code", "numb", "ref_id",
  "accepted_scientific_name", "family", "order", "class",
  "country", "sampling_year", "zoorealms",
  "lat", "long",
  "total_tested", "positive_individuals", "ratio_positive",
  "diagnostic_method", "assay_stage", "diagnostic_group",
  "biological_interpretation"
)

missing_cols <- setdiff(required_cols, names(wnv_long))
if (length(missing_cols) > 0) {
  stop("Missing columns: ", paste(missing_cols, collapse = ", "))
}

# Clean row-level data
wnv_long <- wnv_long %>%
  filter(!ref_id %in% exclude_ref_ids) %>%
  filter(!is.na(accepted_scientific_name), accepted_scientific_name != "NA") %>%
  mutate(
    ref_id = as.numeric(ref_id),
    numb = as.numeric(numb),
    total_tested = as.numeric(total_tested),
    positive_individuals = as.numeric(positive_individuals),
    ratio_positive = as.numeric(ratio_positive),
    sampling_year = as.numeric(sampling_year),
    lat = as.numeric(lat),
    long = as.numeric(long)
  )

# Collapse to ONE ROW PER observation_code
wnv_obs <- wnv_long %>%
  group_by(observation_code) %>%
  summarise(
    numb = first_non_na(as.character(numb)),
    ref_id = first_non_na(as.character(ref_id)),
    accepted_scientific_name = first_non_na(accepted_scientific_name),
    species_name_corrected = first_non_na(species_name_corrected),
    family = first_non_na(family),
    order = first_non_na(order),
    class = first_non_na(class),
    country = first_non_na(country),
    province_region = first_non_na(province_region),
    type_loc = first_non_na(type_loc),
    sampling_year = safe_max_num(sampling_year),
    sampling_period = first_non_na(sampling_period),
    biogeographic_realm = first_non_na(biogeographic_realm),
    zoorealms = first_non_na(zoorealms),
    bird_upgma = first_non_na(as.character(bird_upgma)),
    mam_upgma = first_non_na(as.character(mam_upgma)),
    lat = safe_max_num(lat),
    long = safe_max_num(long),
    total_tested = safe_max_num(total_tested),
    positive_individuals = safe_max_num(positive_individuals),
    positive_observation = safe_any_positive(positive_individuals, ratio_positive),
    ratio_positive_obs = ifelse(
      !is.na(total_tested) & total_tested > 0,
      positive_individuals / total_tested,
      NA_real_
    ),
    diagnostic_methods = paste(sort(unique(na.omit(diagnostic_method))), collapse = "; "),
    assay_roles = paste(sort(unique(na.omit(assay_stage))), collapse = "; "),
    diagnostic_groups = paste(sort(unique(na.omit(diagnostic_group))), collapse = "; "),
    biological_interpretations = paste(sort(unique(na.omit(biological_interpretation))), collapse = "; "),
    source_doi = first_non_na(source_doi),
    .groups = "drop"
  ) %>%
  mutate(
    ref_id = as.numeric(ref_id),
    numb = as.numeric(numb)
  )

# Load world/realm maps once
Sys.setenv(SHAPE_RESTORE_SHX = "YES")
zooregions <- st_read(MAP_REALMS, quiet = TRUE)
if (is.na(st_crs(zooregions))) st_crs(zooregions) <- 4326

world <- ne_countries(scale = "medium", returnclass = "sf")

# ======================================================================================
# STEP 2. GENERAL DESCRIPTIVE STATISTICS
# ======================================================================================

cat("Total unique studies:", n_distinct(wnv_obs$ref_id), "\n")
cat("Total unique observation_code records:", n_distinct(wnv_obs$observation_code), "\n")

studies_by_class <- wnv_obs %>%
  group_by(class) %>%
  summarise(
    unique_studies = n_distinct(ref_id),
    unique_species = n_distinct(accepted_scientific_name),
    total_tested = sum(total_tested, na.rm = TRUE),
    .groups = "drop"
  )

print(studies_by_class)

species_status <- build_species_status(wnv_obs, "class")

unique_species_status <- species_status %>%
  group_by(class, status) %>%
  summarise(n_species = n_distinct(accepted_scientific_name), .groups = "drop")

print(unique_species_status)

# ======================================================================================
# STEP 3. BIRD ORDER SUMMARY TABLE
# ======================================================================================

birds_obs <- wnv_obs %>%
  filter(class == "Aves", !is.na(order), !is.na(accepted_scientific_name))

bird_species_status_order <- birds_obs %>%
  group_by(order, accepted_scientific_name) %>%
  summarise(
    species_positive = any(positive_observation, na.rm = TRUE),
    .groups = "drop"
  )

species_counts_order <- bird_species_status_order %>%
  group_by(order) %>%
  summarise(
    species_tested = n_distinct(accepted_scientific_name),
    species_pos = n_distinct(accepted_scientific_name[species_positive]),
    species_neg = n_distinct(accepted_scientific_name[!species_positive]),
    .groups = "drop"
  )

individuals_counts_order <- birds_obs %>%
  group_by(order) %>%
  summarise(
    individuals_tested = sum(total_tested, na.rm = TRUE),
    .groups = "drop"
  )

bird_order_summary_table <- species_counts_order %>%
  left_join(individuals_counts_order, by = "order") %>%
  arrange(order) %>%
  rename(
    Order = order,
    `Species tested` = species_tested,
    `Species (pos.)` = species_pos,
    `Species (neg.)` = species_neg,
    `Individuals tested` = individuals_tested
  )

print(bird_order_summary_table, n = Inf)
write_xlsx(bird_order_summary_table, "bird_order_summary_table_prevalence.xlsx")

# Wider order summary by class
species_status_order <- build_species_status(wnv_obs, c("order", "class"))

unique_species_status_order <- species_status_order %>%
  group_by(order, class, status) %>%
  summarise(n_species = n_distinct(accepted_scientific_name), .groups = "drop")

total_tested_order <- wnv_obs %>%
  group_by(order, class) %>%
  summarise(total_tested_individuals = sum(total_tested, na.rm = TRUE), .groups = "drop")

order_summary_wide <- unique_species_status_order %>%
  pivot_wider(names_from = status, values_from = n_species, values_fill = 0) %>%
  left_join(total_tested_order, by = c("order", "class"))

print(order_summary_wide, n = 60)
write_xlsx(order_summary_wide, "unique_species_status_by_order_with_total_tested.xlsx")

# ======================================================================================
# STEP 4. SPECIES STATUS BY FAMILY AND CLASS / ZOOREALM
# ======================================================================================

species_status_family_zoorealm <- build_species_status(
  wnv_obs,
  c("family", "class", "zoorealms")
)

family_zoorealm_summary <- species_status_family_zoorealm %>%
  group_by(family, class, zoorealms, status) %>%
  summarise(n_species = n_distinct(accepted_scientific_name), .groups = "drop") %>%
  arrange(class, family, zoorealms, status)

family_zoorealm_summary_wide <- family_zoorealm_summary %>%
  pivot_wider(names_from = status, values_from = n_species, values_fill = 0)

write_xlsx(family_zoorealm_summary_wide, "unique_species_status_family_zoorealm_wide.xlsx")

species_status_class_zoorealm <- build_species_status(
  wnv_obs,
  c("class", "zoorealms")
)

class_zoorealm_summary_wide <- species_status_class_zoorealm %>%
  group_by(class, zoorealms, status) %>%
  summarise(n_species = n_distinct(accepted_scientific_name), .groups = "drop") %>%
  pivot_wider(names_from = status, values_from = n_species, values_fill = 0) %>%
  mutate(
    percentage_positive = Positive / (Positive + Negative) * 100,
    total_species_tested = Positive + Negative
  )

print(class_zoorealm_summary_wide)

# ======================================================================================
# STEP 5. MAP OF PERCENTAGE POSITIVE SPECIES BY ZOOREALM
# ======================================================================================

zooregions_data <- merge(
  zooregions,
  class_zoorealm_summary_wide,
  by.x = "Realm",
  by.y = "zoorealms",
  all.x = TRUE
)

aves_data <- subset(zooregions_data, class == "Aves")
mammalia_data <- subset(zooregions_data, class == "Mammalia")

plot_aves_tested <- ggplot(aves_data) +
  geom_sf(aes(fill = total_species_tested)) +
  scale_fill_gradientn(colors = brewer.pal(9, "Reds"), name = "Species tested", na.value = "lightgrey") +
  theme_classic() +
  labs(title = "Number of Bird Species Tested for WNV by Zoorealm")

plot_aves_positive <- ggplot(aves_data) +
  geom_sf(aes(fill = percentage_positive)) +
  scale_fill_gradientn(colors = brewer.pal(9, "YlOrRd"), name = "% Positive species",
                       limits = c(0, 100), na.value = "lightgrey") +
  theme_classic() +
  labs(title = "Percentage of WNV-Positive Bird Species by Zoorealm")

plot_mammal_tested <- ggplot(mammalia_data) +
  geom_sf(aes(fill = total_species_tested)) +
  scale_fill_gradientn(colors = brewer.pal(9, "Blues"), name = "Species tested", na.value = "lightgrey") +
  theme_classic() +
  labs(title = "Number of Mammal Species Tested for WNV by Zoorealm")

plot_mammal_positive <- ggplot(mammalia_data) +
  geom_sf(aes(fill = percentage_positive)) +
  scale_fill_gradientn(colors = brewer.pal(9, "YlOrRd"), name = "% Positive species",
                       limits = c(0, 100), na.value = "lightgrey") +
  theme_classic() +
  labs(title = "Percentage of WNV-Positive Mammal Species by Zoorealm")

print(plot_aves_tested / plot_aves_positive)
print(plot_mammal_tested / plot_mammal_positive)

# ======================================================================================
# STEP 6. UNIQUE STUDIES BY COUNTRY
# ======================================================================================

recode_country <- c(
  "United States" = "United States of America",
  "Czech Republic" = "Czechia",
  "Arab Republic of Egypt" = "Egypt",
  "Russian Federation" = "Russia",
  "Islamic Republic of Iran" = "Iran",
  "Democratic Republic of the Congo" = "Dem. Rep. Congo",
  "Cote d’Ivoire" = "Côte d'Ivoire",
  "Dominican Republic" = "Dominican Rep.",
  "Central African Republic" = "Central African Rep.",
  "Republic of Paraguay" = "Paraguay",
  "Republic of Trinidad and Tobago" = "Trinidad and Tobago",
  "Reunion" = "Réunion",
  "Western Indian Ocean" = NA_character_
)

wnv_obs_country <- wnv_obs %>%
  mutate(country_new = recode(country, !!!recode_country, .default = country))

studies_by_country <- wnv_obs_country %>%
  filter(!is.na(country_new)) %>%
  group_by(country_new) %>%
  summarise(unique_studies = n_distinct(ref_id), .groups = "drop") %>%
  arrange(desc(unique_studies))

print(studies_by_country, n = 100)

world_data <- left_join(world, studies_by_country, by = c("name" = "country_new")) %>%
  mutate(
    studies_cat = cut(
      unique_studies,
      breaks = c(-Inf, 5, 10, 20, 50, Inf),
      labels = c("0–5", "6–10", "11–20", "21–50", ">50")
    )
  )

map_plot <- ggplot(world_data) +
  geom_sf(aes(fill = studies_cat)) +
  scale_fill_brewer(palette = "Reds", na.value = "grey95") +
  labs(fill = "Unique Studies", title = "Unique WNV Studies by Country") +
  theme_classic()

print(map_plot)

# ======================================================================================
# STEP 7. SITE-LEVEL SUMMARIES
# ======================================================================================

wnv_sites <- wnv_obs %>%
  filter(class %in% c("Aves", "Mammalia")) %>%
  filter(!is.na(lat), !is.na(long)) %>%
  mutate(site = paste(round(lat, 4), round(long, 4), sep = "_"))

site_summary <- wnv_sites %>%
  distinct(observation_code, site, class, lat, long, total_tested, ref_id) %>%
  group_by(site, class, lat, long) %>%
  summarise(
    unique_studies = n_distinct(ref_id),
    individuals_tested_site = sum(total_tested, na.rm = TRUE),
    n_observations = n_distinct(observation_code),
    .groups = "drop"
  )

sites_aves <- site_summary %>% filter(class == "Aves")
sites_mam  <- site_summary %>% filter(class == "Mammalia")

sites_sf_aves <- st_as_sf(sites_aves, coords = c("long", "lat"), crs = st_crs(zooregions))
sites_sf_mam  <- st_as_sf(sites_mam,  coords = c("long", "lat"), crs = st_crs(zooregions))

# Maps: unique studies by site
p_site_studies_aves <- plot_site_map(
  sf_data = sites_sf_aves,
  map_data = zooregions,
  metric = "unique_studies",
  color_hex = "#ab4901",
  title_text = "Bird WNV Studies by Site",
  subtitle_text = "Circle size indicates the number of unique studies per site",
  size_range = c(1, 10),
  legend_name = "Unique studies"
)

p_site_studies_mam <- plot_site_map(
  sf_data = sites_sf_mam,
  map_data = zooregions,
  metric = "unique_studies",
  color_hex = "#065471",
  title_text = "Mammal WNV Studies by Site",
  subtitle_text = "Circle size indicates the number of unique studies per site",
  size_range = c(1, 4),
  legend_name = "Unique studies"
)

print(p_site_studies_aves / p_site_studies_mam)

# Maps: individuals tested by site
bird_breaks   <- c(200, 2000, 20000, 40000)
mammal_breaks <- c(200, 2000, 6000)

p_site_tested_aves <- plot_site_map(
  sf_data = sites_sf_aves,
  map_data = zooregions,
  metric = "individuals_tested_site",
  color_hex = "#ab4901",
  title_text = "Bird WNV Sampling by Site",
  subtitle_text = "Circle size indicates the total number of individuals tested per site",
  size_range = c(1.5, 12),
  legend_name = "Individuals tested",
  legend_breaks = bird_breaks
)

p_site_tested_mam <- plot_site_map(
  sf_data = sites_sf_mam,
  map_data = zooregions,
  metric = "individuals_tested_site",
  color_hex = "#065471",
  title_text = "Mammal WNV Sampling by Site",
  subtitle_text = "Circle size indicates the total number of individuals tested per site",
  size_range = c(1.5, 8),
  legend_name = "Individuals tested",
  legend_breaks = mammal_breaks
)

print(p_site_tested_aves / p_site_tested_mam)

# ======================================================================================
# STEP 8. TOTAL NUMBER OF ANIMALS TESTED BY CLASS AND YEAR
# ======================================================================================

class_colors <- c(
  "Aves" = "#ab4901",
  "Mammalia" = "#065471",
  "Amphibia" = "#cb4c78",
  "Sauropsida" = "#39b54a"
)

wnv_summary <- wnv_obs %>%
  filter(!is.na(sampling_year), !is.na(class)) %>%
  group_by(sampling_year, class) %>%
  summarise(total_tested = sum(total_tested, na.rm = TRUE), .groups = "drop")

year_breaks <- seq(1950, max(wnv_summary$sampling_year, na.rm = TRUE), by = 10)

ggplot(wnv_summary, aes(x = sampling_year, y = total_tested, fill = class)) +
  geom_col(position = "stack", width = 0.8, alpha = 0.7) +
  scale_fill_manual(values = class_colors) +
  scale_x_continuous(breaks = year_breaks) +
  labs(
    title = "Total Number of Animals Tested for WNV by Class and Year",
    x = "Year",
    y = "Number of Animals Tested",
    fill = "Class"
  ) +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# ======================================================================================
# STEP 9. CUMULATIVE FIRST-TIME TESTED AND POSITIVE SPECIES OVER TIME
# ======================================================================================

wnv_filtered <- wnv_obs %>%
  filter(class %in% c("Aves", "Mammalia")) %>%
  filter(!is.na(sampling_year))

wnv_first_tested <- wnv_filtered %>%
  group_by(class, accepted_scientific_name) %>%
  summarise(first_tested_year = min(sampling_year, na.rm = TRUE), .groups = "drop")

wnv_first_positive <- wnv_filtered %>%
  filter(positive_observation) %>%
  group_by(class, accepted_scientific_name) %>%
  summarise(first_positive_year = min(sampling_year, na.rm = TRUE), .groups = "drop")

cumulative_tested <- wnv_first_tested %>%
  count(class, first_tested_year) %>%
  arrange(class, first_tested_year) %>%
  group_by(class) %>%
  mutate(cumulative_species = cumsum(n)) %>%
  rename(sampling_year = first_tested_year) %>%
  mutate(status = "Tested")

cumulative_positive <- wnv_first_positive %>%
  count(class, first_positive_year) %>%
  arrange(class, first_positive_year) %>%
  group_by(class) %>%
  mutate(cumulative_species = cumsum(n)) %>%
  rename(sampling_year = first_positive_year) %>%
  mutate(status = "WNV-Positive")

cumulative_combined <- bind_rows(cumulative_tested, cumulative_positive)

ggplot(cumulative_combined,
       aes(x = sampling_year, y = cumulative_species, color = class, linetype = status)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  scale_color_manual(values = c("Aves" = "#ab4901", "Mammalia" = "#065471")) +
  scale_linetype_manual(values = c("Tested" = "dotted", "WNV-Positive" = "solid")) +
  scale_x_continuous(breaks = year_breaks) +
  labs(
    title = "Cumulative First-Time Tested and WNV-Positive Species Over Time",
    subtitle = "Observation-based counts (one count per observation_code)",
    x = "Sampling Year",
    y = "Cumulative Number of Species",
    color = "Class",
    linetype = "Species status"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 14),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(linewidth = 0.5),
    legend.position = "top"
  )

plot_cumulative_by_class <- function(target_class, fill_color) {
  
  dat_class <- wnv_obs %>%
    filter(class == target_class, !is.na(sampling_year))
  
  first_tested <- dat_class %>%
    group_by(accepted_scientific_name) %>%
    summarise(first_tested_year = min(sampling_year, na.rm = TRUE), .groups = "drop")
  
  first_positive <- dat_class %>%
    filter(positive_observation) %>%
    group_by(accepted_scientific_name) %>%
    summarise(first_positive_year = min(sampling_year, na.rm = TRUE), .groups = "drop")
  
  cumulative_tested <- first_tested %>%
    count(first_tested_year) %>%
    arrange(first_tested_year) %>%
    mutate(cumulative_species = cumsum(n), status = "Tested") %>%
    rename(sampling_year = first_tested_year)
  
  cumulative_positive <- first_positive %>%
    count(first_positive_year) %>%
    arrange(first_positive_year) %>%
    mutate(cumulative_species = cumsum(n), status = "WNV-Positive") %>%
    rename(sampling_year = first_positive_year)
  
  cumulative_combined <- bind_rows(cumulative_tested, cumulative_positive)
  
  ggplot(cumulative_combined, aes(x = sampling_year, y = cumulative_species, color = status)) +
    geom_col(
      data = subset(cumulative_combined, status == "Tested"),
      fill = fill_color,
      alpha = 0.5,
      color = NA
    ) +
    geom_line(
      data = subset(cumulative_combined, status == "WNV-Positive"),
      aes(group = 1),
      linewidth = 1,
      color = fill_color
    ) +
    geom_point(
      data = subset(cumulative_combined, status == "WNV-Positive"),
      size = 3,
      color = fill_color
    ) +
    labs(
      title = paste("Cumulative First-Time Tested and WNV-Positive Species:", target_class),
      x = "Sampling Year",
      y = "Cumulative Number of Species"
    ) +
    theme_minimal() +
    theme(
      text = element_text(size = 12),
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(linewidth = 0.5),
      legend.position = "none"
    ) +
    scale_x_continuous(limits = c(1950, 2025), breaks = seq(1950, 2025, by = 10))
}

p_aves <- plot_cumulative_by_class("Aves", "#ab4901")
p_mammals <- plot_cumulative_by_class("Mammalia", "#065471")
print(p_aves / p_mammals)

# ======================================================================================
# STEP 10. BIRD METHOD COMBINATIONS (USE wnv_long)
# ======================================================================================

birds_long <- wnv_long %>%
  filter(class == "Aves")

all_methods <- birds_long %>%
  count(diagnostic_method, sort = TRUE)

print(all_methods, n = Inf)

birds_reclassified <- birds_long %>%
  mutate(
    method_family = case_when(
      diagnostic_method %in% c("serological.ELISA", "serological", "HI", "IgC", "IgM", "MIA") ~ "serological.ELISA",
      diagnostic_method %in% c("serological.VNT") ~ "serological.VNT",
      diagnostic_method %in% c("molecular_detection", "IHC") ~ "molecular_detection",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(method_family))

species_method <- birds_reclassified %>%
  group_by(accepted_scientific_name, method_family) %>%
  summarise(
    tested_by_method = TRUE,
    positive_by_method = any(
      (!is.na(positive_individuals) & positive_individuals > 0) |
        (!is.na(ratio_positive) & ratio_positive > 0),
      na.rm = TRUE
    ),
    .groups = "drop"
  )

species_methods_wide <- species_method %>%
  pivot_wider(
    names_from = method_family,
    values_from = c(tested_by_method, positive_by_method),
    values_fill = FALSE
  )

needed_cols <- c(
  "tested_by_method_serological.ELISA",
  "tested_by_method_serological.VNT",
  "tested_by_method_molecular_detection",
  "positive_by_method_serological.ELISA",
  "positive_by_method_serological.VNT",
  "positive_by_method_molecular_detection"
)

for (cc in needed_cols) {
  if (!cc %in% names(species_methods_wide)) species_methods_wide[[cc]] <- FALSE
}

species_methods_wide <- species_methods_wide %>%
  mutate(
    tested_ELISA = tested_by_method_serological.ELISA,
    tested_VNT   = tested_by_method_serological.VNT,
    tested_MOL   = tested_by_method_molecular_detection,
    positive_ELISA = positive_by_method_serological.ELISA,
    positive_VNT   = positive_by_method_serological.VNT,
    positive_MOL   = positive_by_method_molecular_detection,
    combo = case_when(
      tested_ELISA & tested_MOL & tested_VNT ~ "ELISA + Molecular + VNT",
      !tested_ELISA & tested_MOL & tested_VNT ~ "Molecular + VNT",
      tested_ELISA & !tested_MOL & tested_VNT ~ "ELISA + VNT",
      tested_ELISA & tested_MOL & !tested_VNT ~ "ELISA + Molecular",
      tested_ELISA & !tested_MOL & !tested_VNT ~ "ELISA only",
      !tested_ELISA & tested_MOL & !tested_VNT ~ "Molecular only",
      !tested_ELISA & !tested_MOL & tested_VNT ~ "VNT only",
      TRUE ~ NA_character_
    ),
    positive_all_tested_methods = case_when(
      combo == "ELISA + Molecular + VNT" ~ (positive_ELISA & positive_MOL & positive_VNT),
      combo == "Molecular + VNT"         ~ (positive_MOL & positive_VNT),
      combo == "ELISA + VNT"             ~ (positive_ELISA & positive_VNT),
      combo == "ELISA + Molecular"       ~ (positive_ELISA & positive_MOL),
      combo == "ELISA only"              ~ positive_ELISA,
      combo == "Molecular only"          ~ positive_MOL,
      combo == "VNT only"                ~ positive_VNT,
      TRUE ~ FALSE
    )
  )

combo_summary <- species_methods_wide %>%
  filter(!is.na(combo)) %>%
  group_by(combo) %>%
  summarise(
    n_species_tested = n_distinct(accepted_scientific_name),
    n_species_positive_strict = n_distinct(accepted_scientific_name[positive_all_tested_methods]),
    .groups = "drop"
  ) %>%
  arrange(desc(n_species_tested), desc(n_species_positive_strict), combo)

print(combo_summary, n = Inf)

combo_summary_plot <- species_methods_wide %>%
  filter(!is.na(combo)) %>%
  transmute(
    accepted_scientific_name,
    combo,
    ELISA = tested_ELISA,
    MOL   = tested_MOL,
    VNT   = tested_VNT,
    positive_strict = positive_all_tested_methods
  ) %>%
  group_by(combo, ELISA, MOL, VNT) %>%
  summarise(
    n_species_tested = n_distinct(accepted_scientific_name),
    n_species_positive_strict = n_distinct(accepted_scientific_name[positive_strict]),
    .groups = "drop"
  ) %>%
  arrange(desc(n_species_tested), desc(n_species_positive_strict), combo) %>%
  mutate(
    x = seq_len(n()),
    combo_short = combo
  )

p_top <- ggplot(combo_summary_plot, aes(x = x)) +
  geom_col(aes(y = n_species_tested, fill = "Bird species tested"), width = 0.88) +
  geom_col(aes(y = n_species_positive_strict, fill = "Bird species positive for WNV"), width = 0.88) +
  geom_text(aes(y = n_species_tested + 8, label = n_species_tested), size = 5.5) +
  geom_text(aes(y = n_species_positive_strict / 2, label = n_species_positive_strict),
            color = "white", size = 4.8) +
  scale_fill_manual(
    values = c("Bird species tested" = "#c4814b",
               "Bird species positive for WNV" = "#c40000"),
    name = NULL
  ) +
  scale_x_continuous(breaks = combo_summary_plot$x, labels = rep("", nrow(combo_summary_plot))) +
  labs(y = "Species count", x = NULL) +
  theme_classic(base_size = 18) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = c(0.77, 0.87),
    legend.text = element_text(size = 14),
    plot.margin = margin(5, 20, 5, 20)
  )

matrix_long <- combo_summary_plot %>%
  select(x, ELISA, MOL, VNT) %>%
  pivot_longer(cols = c(ELISA, MOL, VNT), names_to = "method", values_to = "present") %>%
  mutate(
    y = case_when(
      method == "ELISA" ~ 1,
      method == "MOL"   ~ 2,
      method == "VNT"   ~ 3
    )
  )

segments_df <- combo_summary_plot %>%
  rowwise() %>%
  mutate(
    ymin = min(c(ifelse(ELISA, 1, NA), ifelse(MOL, 2, NA), ifelse(VNT, 3, NA)), na.rm = TRUE),
    ymax = max(c(ifelse(ELISA, 1, NA), ifelse(MOL, 2, NA), ifelse(VNT, 3, NA)), na.rm = TRUE)
  ) %>%
  ungroup()

p_matrix <- ggplot() +
  annotate("rect", xmin = 0.4, xmax = nrow(combo_summary_plot) + 0.6,
           ymin = 1.5, ymax = 2.5, fill = "#22a9e1", alpha = 0.85) +
  geom_segment(data = segments_df,
               aes(x = x, xend = x, y = ymin, yend = ymax), linewidth = 1) +
  geom_point(data = matrix_long,
             aes(x = x, y = y, fill = present),
             shape = 21, size = 6, color = "black", stroke = 0.3) +
  scale_fill_manual(values = c("TRUE" = "#0d5d83", "FALSE" = "#cfcfcf")) +
  scale_x_continuous(breaks = combo_summary_plot$x, labels = rep("", nrow(combo_summary_plot))) +
  scale_y_continuous(
    breaks = c(3, 2, 1),
    labels = c("serological.VNT", "molecular_detection", "serological.ELISA"),
    limits = c(0.4, 3.6)
  ) +
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = 16) +
  theme(
    legend.position = "none",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 15),
    plot.margin = margin(5, 20, 5, 20)
  )

p_bottom <- ggplot(combo_summary_plot, aes(x = x, y = 1, label = combo_short)) +
  geom_text(angle = 45, hjust = 1, vjust = 1, size = 5) +
  scale_x_continuous(limits = c(0.5, nrow(combo_summary_plot) + 0.5)) +
  coord_cartesian(ylim = c(0.5, 1.5), clip = "off") +
  labs(x = "WNV Detection Method", y = NULL) +
  theme_void(base_size = 16) +
  theme(
    axis.title.x = element_text(size = 20, face = "bold", margin = margin(t = 30)),
    plot.margin = margin(0, 20, 10, 20)
  )

final_plot <- p_top / p_matrix / p_bottom +
  plot_layout(heights = c(2.8, 1.5, 0.9))

print(final_plot)
ggsave("bird_method_combinations_strict_ordered.png", final_plot,
       width = 12, height = 10, dpi = 300)
