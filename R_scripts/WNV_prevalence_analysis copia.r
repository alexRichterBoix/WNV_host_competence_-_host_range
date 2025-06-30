# Load necessary libraries
library(dplyr)
library(readxl)
library(writexl)
library(metafor)
library(binom)
library(tidyr)
library(writexl)
library(ggplot2)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(patchwork)
library(superheat)
library(tidyverse)
library(ComplexUpset)
library(ape)
library(ggtree)
library(RColorBrewer)



#--------- Descriptive statistics------------------------------------

# Read in your data
wnv_data <- read_excel("WNV_Host_Susceptability_Review.xlsx")  # Adjust the path if necessary

# Clean the data (remove references of lists of species without "prevalence information" only species detection)
wnv_data_clean <- wnv_data %>%
  filter(!ref.id...2 %in% c(304, 222, 323))

# Remove rows with missing or literal "NA" in accepted_scientific_name
#    and convert total.tested & positive.individuals.m1 to numeric
wnv_data_clean <- wnv_data_clean %>%
  filter(!is.na(accepted_scientific_name) & accepted_scientific_name != "NA") %>%
  mutate(
    total.tested = as.numeric(total.tested),
    positive.individuals.m1 = as.numeric(positive.individuals.m1)
  )

library(dplyr)

# 1. Total number of unique studies in the whole dataset:
total_studies <- wnv_data_clean %>% 
  distinct(ref.id...2) %>% 
  nrow()
cat("Total unique studies:", total_studies, "\n")


#1.2 Total number of records in the whole dataset:
total_records <- wnv_data_clean %>% 
  distinct(numb) %>% 
  nrow()
cat("Total unique records:", total_records, "\n")


# 2. Number of unique studies for class = "Aves":
aves_studies <- wnv_data_clean %>% 
  filter(class == "Aves") %>% 
  distinct(ref.id...2) %>% 
  nrow()
cat("Unique studies for Aves:", aves_studies, "\n")


# 3. Number of unique studies for class = "Mammalia":
mammalia_studies <- wnv_data_clean %>% 
  filter(class == "Mammalia") %>% 
  distinct(ref.id...2) %>% 
  nrow()
cat("Unique studies for Mammalia:", mammalia_studies, "\n")


# 4. Count unique species per class:
unique_species_by_class <- wnv_data_clean %>%
   group_by(class) %>%
   summarize(unique_species = n_distinct(accepted_scientific_name), .groups = "drop")
print("Unique species by class:")
print(unique_species_by_class)


# Flag species as positive if at least one record has ratio.m1 > 0
species_status <- wnv_data_clean %>%
  group_by(class, accepted_scientific_name) %>%
  summarize(positive = any(ratio.m1 > 0, na.rm = TRUE), .groups = "drop") %>%
  mutate(status = if_else(positive, "Positive", "Negative"))

# Count unique species per class by status
unique_species_status <- species_status %>%
  group_by(class, status) %>%
  summarize(n_species = n_distinct(accepted_scientific_name), .groups = "drop")
print(unique_species_status)


# Flag each species as positive if at least one record has ratio.m1 > 0, grouping by order, class, and species.
species_status_order <- wnv_data_clean %>%
  group_by(order, class, accepted_scientific_name) %>%
  summarize(positive = any(ratio.m1 > 0, na.rm = TRUE), .groups = "drop") %>%
  mutate(status = if_else(positive, "Positive", "Negative"))

# Count unique species per order (with class information) by status.
unique_species_status_order <- species_status_order %>%
  group_by(order, class, status) %>%
  summarize(n_species = n_distinct(accepted_scientific_name), .groups = "drop")

# Calculate the total number of individuals tested per order
total_tested_order <- wnv_data_clean %>%
  group_by(order, class) %>%
  summarize(total_tested_individuals = sum(total.tested, na.rm = TRUE), .groups = "drop")

# Merge the species status by order and the total number of tested individuals
pivoted_results <- unique_species_status_order %>%
  pivot_wider(names_from = status, values_from = n_species) %>%
  left_join(total_tested_order, by = c("order", "class"))

# View the final table with the total number of individuals tested per order
print(pivoted_results, n = 55)


# Save the output as an Excel file
write_xlsx(pivoted_results, "unique_species_status_by_order_wide_with_total_tested.xlsx")



# 5. Count the number of unique "numb" values for class "Aves" and class "Mammalia"
# (Assuming 'numb' is a unique identifier for each row/study)
numb_counts <- wnv_data_clean %>%
   filter(class %in% c("Aves", "Mammalia")) %>%
   group_by(class) %>%
   summarize(numb_count = n_distinct(numb), .groups = "drop")
   
print("Unique 'numb' counts for Aves and Mammalia:")
print(numb_counts)


# 6. Number of unique studies per unique country:
studies_by_country <- wnv_data_clean %>%
  group_by(country) %>%
  summarize(unique_studies = n_distinct(ref.id...2), .groups = "drop") %>%
  arrange(desc(unique_studies))

print(studies_by_country, n=91)



###-----------------------------------------------------------------------

# Number of negative & positive species by family and zoorealms

# 1. Flag each species as positive if at least one record has ratio.m1 > 0
species_status_family_zoorealm <- wnv_data_clean %>%
  group_by(family, class, accepted_scientific_name, zoorealms) %>%
  summarize(
    positive = any(ratio.m1 > 0, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(status = if_else(positive, "Positive", "Negative"))

# 2. Count unique species per family, class, zoorealm, and status
unique_species_status_family_zoorealm <- species_status_family_zoorealm %>%
  group_by(family, class, zoorealms, status) %>%
  summarize(
    n_species = n_distinct(accepted_scientific_name),
    .groups = "drop"
  )

# 3. (Optional) Arrange the output nicely
unique_species_status_family_zoorealm <- unique_species_status_family_zoorealm %>%
  arrange(class, family, zoorealms, status)

# View the table
print(unique_species_status_family_zoorealm)

# Pivot wider: status becomes columns
unique_species_status_family_zoorealm_wide <- unique_species_status_family_zoorealm %>%
  pivot_wider(
    names_from = status,
    values_from = n_species,
    values_fill = 0  # Fill missing combinations with zero
  )

# View the result
print(unique_species_status_family_zoorealm_wide)

# Save the table as an Excel file
write_xlsx(unique_species_status_family_zoorealm_wide, "unique_species_status_family_zoorealm_wide.xlsx")



###-----------------------------------------------------------------------

# Number of negative & positive species by class and zoorealms

# 1. Flag each species as positive if at least one record has ratio.m1 > 0
species_status_class_zoorealm <- wnv_data_clean %>%
  group_by(class, accepted_scientific_name, zoorealms) %>%
  summarize(
    positive = any(ratio.m1 > 0, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(status = if_else(positive, "Positive", "Negative"))

# 2. Count unique species per family, class, zoorealm, and status
unique_species_status_class_zoorealm <- species_status_class_zoorealm %>%
  group_by(class, zoorealms, status) %>%
  summarize(
    n_species = n_distinct(accepted_scientific_name),
    .groups = "drop"
  )

# 3. (Optional) Arrange the output nicely
unique_species_status_class_zoorealm <- unique_species_status_class_zoorealm %>%
  arrange(class, zoorealms, status)

# View the table
print(unique_species_status_class_zoorealm)

# Pivot wider: status becomes columns
unique_species_status_class_zoorealm_wide <- unique_species_status_class_zoorealm %>%
  pivot_wider(
    names_from = status,
    values_from = n_species,
    values_fill = 0  # Fill missing combinations with zero
  )

# View the result
print(unique_species_status_class_zoorealm_wide)


# Calculate the percentage of positive species for each zoorealm
unique_species_status_class_zoorealm_wide <- unique_species_status_class_zoorealm_wide %>%
  mutate(percentage_positive = Positive / (Positive + Negative) * 100)


# Calculate species tested for each zoorealm
unique_species_status_class_zoorealm_wide <- unique_species_status_class_zoorealm_wide %>%
  mutate(total_tested = Positive + Negative)


# 3. Call the zooregions map
Sys.setenv(SHAPE_RESTORE_SHX = "YES")
library(sf)
zooregions <- st_read("realms_maps/newRealms.shp")

# Check the structure of the data to confirm the attribute columns are loaded
print(zooregions)

# Check the column names to ensure the relevant attribute (e.g., "zoorealms") is available
colnames(zooregions)

# Merge the zooregions shapefile with the species data
zooregions_data <- merge(zooregions, unique_species_status_class_zoorealm_wide, by.x = "Realm", by.y = "zoorealms", all.x = TRUE)

# Check the merged data
print(zooregions_data)

# Plot the map with percentage of positive species for WNV

# Filter data for 'Aves' class
aves_data <- subset(zooregions_data, class == "Aves")

# Filter data for 'Mammalia' class
mammalia_data <- subset(zooregions_data, class == "Mammalia")


# Create the plot for 'Aves' class
plot_aves_tested <- ggplot(data = aves_data) +
  geom_sf(aes(fill = total_tested)) +  # Fill regions based on the percentage
  scale_fill_gradientn(
    colors = brewer.pal(9, "Reds"),  # Use the regular color scale where higher values are darker
    name = "Species tested", 
    limits = c(),            # Set the color scale from 0 to 100 to reflect percentage
    na.value = "lightgrey"              # Show polygons without data in grey
  ) + 
  theme_classic() +
  labs(
    title = "Number of Species tested for WNV in Aves Zoorealm",
    subtitle = "",
    caption = ""
  ) +
  theme(
    legend.position = "right",    # Position the legend to the right
    legend.key.size = unit(0.5, "cm")  # Reduce the size of the legend keys
  )

# Create the plot for 'Aves' class
plot_aves <- ggplot(data = aves_data) +
  geom_sf(aes(fill = percentage_positive)) +  # Fill regions based on the percentage
  scale_fill_gradientn(
    colors = brewer.pal(9, "YlOrRd"),  # Use the regular color scale where higher values are darker
    name = "% of Positive Species", 
    limits = c(0, 100),            # Set the color scale from 0 to 100 to reflect percentage
    na.value = "lightgrey"              # Show polygons without data in grey
  ) + 
  theme_classic() +
  labs(
    title = "Percentage of Positive Species for WNV in Aves Zoorealm",
    subtitle = "",
    caption = ""
  ) +
  theme(
    legend.position = "right",    # Position the legend to the right
    legend.key.size = unit(0.5, "cm")  # Reduce the size of the legend keys
  )

# Combine both plots into a single layout using patchwork
combined_plot_1 <- plot_aves_tested / plot_aves

# Display the combined plot
print(combined_plot_1)

# Create the plot for 'Aves' class
plot_mammalia_tested <- ggplot(data = mammalia_data) +
  geom_sf(aes(fill = total_tested)) +  # Fill regions based on the percentage
  scale_fill_gradientn(
    colors = brewer.pal(9, "Blues"),  # Use the regular color scale where higher values are darker
    name = "Species tested", 
    limits = c(),            # Set the color scale from 0 to 100 to reflect percentage
    na.value = "lightgrey"              # Show polygons without data in grey
  ) + 
  theme_classic() +
  labs(
    title = "Number of Species tested for WNV in Mammals by Zoorealm",
    subtitle = "",
    caption = ""
  ) +
  theme(
    legend.position = "right",    # Position the legend to the right
    legend.key.size = unit(0.5, "cm")  # Reduce the size of the legend keys
  )

# Create the plot for 'Mammalia' class
plot_mammalia <- ggplot(data = mammalia_data) +
  geom_sf(aes(fill = percentage_positive)) +  # Fill regions based on the percentage
  scale_fill_gradientn(
    colors = brewer.pal(9, "YlOrRd"),  # Use the regular color scale where higher values are darker
    name = "% of Positive Species", 
    limits = c(0, 100),            # Set the color scale from 0 to 100 to reflect percentage
    na.value = "lightgrey"              # Show polygons without data in grey
  ) + 
  theme_classic() +
  labs(
    title = "Percentage of Positive Species for WNV in Mammalia Zoorealm",
    subtitle = "",
    caption = ""
  ) +
  theme(
    legend.position = "right",    # Position the legend to the right
    legend.key.size = unit(0.5, "cm")  # Reduce the size of the legend keys
  )

# Combine both plots into a single layout using patchwork
combined_plot <- plot_mammalia_tested / plot_mammalia

# Display the combined plot
print(combined_plot)

# Combine both plots into a single layout using patchwork
combined_plot_3 <- plot_aves_tested / plot_mammalia_tested

# Display the combined plot
print(combined_plot_3)


#-------------Same but instead of zoorealms, % by bird upgm
# Step 1: Ensure species status is computed correctly
species_status_class_zoorealm <- wnv_data_clean %>%
  group_by(class, accepted_scientific_name, bir_upgma_) %>%
  summarize(
    positive = any(ratio.m1 > 0, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(status = if_else(positive, "Positive", "Negative"))

# Step 2: Count unique species per family, class, bir_upgma, and status
unique_species_status_class_zoorealm <- species_status_class_zoorealm %>%
  group_by(class, bir_upgma_, status) %>%
  summarize(
    n_species = n_distinct(accepted_scientific_name),
    .groups = "drop"
  )

# Step 3: Pivot to wide format (status becomes columns)
unique_species_status_class_zoorealm_wide <- unique_species_status_class_zoorealm %>%
  pivot_wider(
    names_from = status,
    values_from = n_species,
    values_fill = 0  # Fill missing combinations with zero
  )

str(unique_species_status_class_zoorealm_wide)

# Step 4: Calculate percentage of positive species for each bir_upgma
unique_species_status_class_zoorealm_wide <- unique_species_status_class_zoorealm_wide %>%
  mutate(percentage_positive = Positive / (Positive + Negative) * 100)

# Step 5: Read the shapefile
zooregions <- st_read("realms_maps/bir.shp")
str(zooregions)

# Check the structure of zooregions to make sure the column exists
print(colnames(zooregions))  # Check if 'bir_upgma_' exists in zooregions

# If 'bir_upgma_' exists in zooregions, merge the datasets correctly
zooregions_data <- merge(zooregions, unique_species_status_class_zoorealm_wide, 
                         by.x = "bir_upgma_", by.y = "bir_upgma_", all.x = TRUE)

# Inspect the merged data to ensure it looks correct
print(zooregions_data)

# Step 7: Filter data for 'Aves' and 'Mammalia' classes
aves_data <- subset(zooregions_data, class == "Aves")
mammalia_data <- subset(zooregions_data, class == "Mammalia")

# Step 8: Plot for 'Aves' class
plot_aves <- ggplot(data = aves_data) +
  geom_sf(aes(fill = percentage_positive)) +  # Fill regions based on the percentage
  scale_fill_gradientn(
    colors = brewer.pal(9, "YlOrRd"),  # Use the regular color scale where higher values are darker
    name = "% of Positive Species", 
    limits = c(0, 100),            # Set the color scale from 0 to 100 to reflect percentage
    na.value = "lightgrey"              # Show polygons without data in grey
  ) + 
  theme_classic() +
  labs(
    title = "Percentage of Positive Species for WNV in Aves Zoorealm",
    subtitle = "",
    caption = ""
  ) +
  theme(
    legend.position = "right",    # Position the legend to the right
    legend.key.size = unit(0.5, "cm")  # Reduce the size of the legend keys
  )

# Display the plot for Aves
print(plot_aves)

# Step 9: Plot for 'Mammalia' class (similar to Aves)
plot_mammalia <- ggplot(data = mammalia_data) +
  geom_sf(aes(fill = percentage_positive)) +  # Fill regions based on the percentage
  scale_fill_gradientn(
    colors = brewer.pal(9, "YlOrRd"),  # Inverted color scale (darker colors for higher %)
    name = "% of Positive Species", 
    limits = c(0, 100),            # Set the color scale from 0 to 100 to reflect percentage
    na.value = "lightgrey"              # Show polygons without data in light grey
  ) + 
  theme_classic() +
  labs(
    title = "Percentage of Positive Species for WNV in Mammalia Zoorealm",
    subtitle = "",
    caption = ""
  ) +
  theme(
    legend.position = "right",    # Position the legend to the right
    legend.key.size = unit(0.5, "cm")  # Reduce the size of the legend keys
  )

# Display the plot for Mammalia
print(plot_mammalia)

# Step 10: Combine the plots for Aves and Mammalia
combined_plot <- plot_aves / plot_mammalia

# Display the combined plot
print(combined_plot)





## -----------------------------------------------------------------------------
### Plot a map with the number of studies by country

# Your studies_by_country data, already ordered:
studies_by_country <- wnv_data_clean %>%
  group_by(country) %>%
  summarize(unique_studies = n_distinct(ref.id...2), .groups = "drop") %>%
  arrange(desc(unique_studies))

# Load the world map (as an sf object) at medium scale
world <- ne_countries(scale = "medium", returnclass = "sf")

# Get the unique country names in your studies data
studies_countries <- unique(studies_by_country$country)
print("Country names in studies data:")
print(studies_countries)

# Get the unique country names in the world map data
world_countries <- unique(world$name)
print("Country names in world map:")
print(world_countries)

# Find the country names that are in your studies data but not in the world map
missing_names <- setdiff(studies_countries, world_countries)
cat("The following country names from your studies data are not found in the world map:\n")
print(missing_names)

# Create a named vector to map from problematic names in your data to the names in world_countries.

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
  "Republic of Trinidad and Tobago" = "Trinidad and Tobago",
  "Reunion" = "Réunion",         # adjust if your world map uses "Réunion"
  "Western Indian Ocean" = NA     # not a country; set to NA
)

# Create a new column 'country_new' that recodes names where needed.
wnv_data_clean <- wnv_data_clean %>%
  mutate(country_new = recode(country, !!!recode_country))

# Check which country names have been changed:
changed_names <- wnv_data_clean %>%
  filter(country != country_new) %>%
  distinct(country, country_new)

cat("The following country names were recoded:\n")
print(changed_names)

# Count unique studies per recoded country name
studies_by_country_new <- wnv_data_clean %>%
  group_by(country_new) %>%
  summarize(unique_studies = n_distinct(ref.id...2), .groups = "drop") %>%
  arrange(desc(unique_studies))
  
print(studies_by_country_new, n=91)

# Join studies data with the world map, matching world$name with country_new
world_data <- left_join(world, studies_by_country_new, by = c("name" = "country_new"))

# Create bins for unique_studies
world_data <- world_data %>%
  mutate(
    studies_cat = cut(unique_studies,
                      breaks = c(-Inf, 5, 10, 20, 50, Inf),
                      labels = c("0–5", "6–10", "11–20", "21–50", ">50"))
  )

# Plot the map with a categorical color scale using scale_fill_brewer
map_plot <- ggplot(world_data) +
  geom_sf(aes(fill = studies_cat)) +
  scale_fill_brewer(palette = "Reds", na.value = "grey95") +
  labs(fill = "Unique Studies",
       title = "Unique WNV Studies by Country") +
  theme_classic()

print(map_plot)


#-----------Create a map with number of studies, animal tested and prevalence (combined serology & molecular) for each unique site to plot maps-------

# 1. Identify the unique sites in the dataset using the latitude and longitude coordinates of the studies
wnv_data_clean <- wnv_data_clean %>%
  mutate(site = paste(round(lat, 4), round(long, 4), sep = "_"))

# 2. Identify the number of studies and total number of animals tested for each unique site
site_summary <- wnv_data_clean %>%
   # Keep only birds and mammals
   filter(class %in% c("Aves", "Mammalia")) %>%
   # Group by site, class, and keep lat/long for plotting
   group_by(site, class, lat, long) %>%
   summarize(
     unique_studies = n_distinct(`ref.id...2`),
     total_tested = sum(total.tested, na.rm = TRUE),
     .groups = "drop"
   )
 
print(site_summary)


# 3. Call the zooregions map
Sys.setenv(SHAPE_RESTORE_SHX = "YES")
library(sf)
zooregions <- st_read("realms_maps/newRealms.shp")
print(zooregions)

plot(zooregions)

# 4. Filter the site summary to only include birds and mammals
# ------------------ For Birds (Aves) ------------------
sites_aves <- site_summary %>%
  filter(class == "Aves")

sites_sf_aves <- st_as_sf(sites_aves, coords = c("long", "lat"), crs = st_crs(zooregions))

p1 <- ggplot() +
  geom_sf(data = zooregions, fill = "grey90", color = "black") +
  geom_sf(data = sites_sf_aves, aes(size = unique_studies),
          color = "red", fill = NA, lwd = 1.2, alpha = 0.3) +
  scale_size_continuous(name = "Unique Studies", range = c(2, 8)) +
  labs(title = "Bird WNV Studies by Site",
       subtitle = "Sites (class = Aves) overlayed on zooregions") +
  theme_minimal()

# ------------------ For Mammals (Mammalia) ------------------
sites_mam <- site_summary %>%
  filter(class == "Mammalia")

sites_sf_mam <- st_as_sf(sites_mam, coords = c("long", "lat"), crs = st_crs(zooregions))

p2 <- ggplot() +
  geom_sf(data = zooregions, fill = "grey90", color = "black") +
  geom_sf(data = sites_sf_mam, aes(size = unique_studies),
          color = "blue", fill = NA, lwd = 1.2, alpha = 0.3) +
  scale_size_continuous(name = "Unique Studies", range = c(1, 4)) +
  labs(title = "Mammalia WNV Studies by Site",
       subtitle = "Sites (class = Mammalia) overlayed on zooregions") +
  theme_minimal()

# ------------------ Combine the two plots vertically ------------------
combined_plot <- p1 / p2  # using patchwork: aves on top, mammals below

print(combined_plot)


# For Birds: Total Tested (using sites_sf_aves)
p1b <- ggplot() +
  geom_sf(data = zooregions, fill = "grey90", color = "black") +
  geom_sf(data = sites_sf_aves, aes(size = total_tested), 
          color = "#ab4901", fill = NA, lwd = 1.2, alpha = 0.4) +
  scale_size_continuous(
    name = "Total Tested",
    range = c(2, 18),
    breaks = c(200, 2000, 20000, 40000),
    labels = function(x) format(x, scientific = FALSE)
  ) +
  labs(title = "Birds: Total Tested",
       subtitle = "Sites (class = Aves) overlayed on zooregions") +
  theme_classic()

# For Mammals: Total Tested (using sites_sf_mam)
p2b <- ggplot() +
  geom_sf(data = zooregions, fill = "grey90", color = "black") +
  geom_sf(data = sites_sf_mam, aes(size = total_tested), 
          color = "#065471", fill = NA, lwd = 1.2, alpha = 0.4) +
  scale_size_continuous(
    name = "Total Tested",
    range = c(1, 5),
    breaks = c(200, 2000, 6000),
    labels = function(x) format(x, scientific = FALSE)
  ) +
  labs(title = "Mammals: Total Tested",
       subtitle = "Sites (class = Mammalia) overlayed on zooregions") +
  theme_classic()

# Combine the two plots vertically: Aves on top, Mammals below.
combined_plot2 <- p1b / p2b

print(combined_plot2)

# Plot a zoom of Europe: sampled animals per site
if(is.na(st_crs(zooregions))) {
  st_crs(zooregions) <- 4326
}

coord_sf(xlim = c(-25, 50), ylim = c(20, 80))
p_zoom1 <- ggplot() +
  geom_sf(data = zooregions, fill = "grey90", color = "black") +
  geom_sf(data = sites_sf_aves, aes(size = total_tested), 
          color = "red", fill = NA, lwd = 1.2, alpha = 0.3) +
  scale_size_continuous(name = "Total Tested", range = c(3, 20),
                        breaks = pretty(sites_sf_aves$total_tested, n = 5)) +
  labs(title = "Birds: Total Tested (Zoomed to Palearctic)",
       subtitle = "Sites (class = Aves) overlayed on zooregions") +
  coord_sf(xlim = c(-25, 70), ylim = c(20, 65)) +
  theme_minimal()
print(p_zoom1)

coord_sf(xlim = c(-150, -50), ylim = c(0, 60))
p_zoom2 <- ggplot() +
  geom_sf(data = zooregions, fill = "grey90", color = "black") +
  geom_sf(data = sites_sf_aves, aes(size = total_tested), 
          color = "red", fill = NA, lwd = 1.2, alpha = 0.3) +
  scale_size_continuous(name = "Total Tested", range = c(3, 20),
                        breaks = pretty(sites_sf_aves$total_tested, n = 5)) +
  labs(title = "Birds: Total Tested (Zoomed to Nearctic)",
       subtitle = "Sites (class = Aves) overlayed on zooregions") +
  coord_sf(xlim = c(-150, -50), ylim = c(0, 60)) +
  theme_minimal()
print(p_zoom2)

combined_plot3 <- p_zoom1 / p_zoom2

#-------------------------------------------------------------------------------

#----Plot the sites zoomed for the Palearctic and Nearctic regions for the zoorealms determined for birds
birdregions <- st_read("realms_maps/bir.shp")

library(sf)
library(ggplot2)
library(RColorBrewer)

# Ensure that birdregions has a defined CRS (if not, set it to WGS84)
if(is.na(st_crs(birdregions))) {
  st_crs(birdregions) <- 4326
}

# Ensure that sites_sf_aves has a defined CRS (set it to the same as birdregions)
if(is.na(st_crs(sites_sf_aves))) {
  st_crs(sites_sf_aves) <- st_crs(birdregions)
}

n_regions <- length(unique(birdregions$bir_upgma_))  # replace with actual column
birdregions$fill_color <- factor(birdregions$bir_upgma_)  # group variable
pastel_colors <- brewer.pal(min(8, n_regions), "Pastel1")

p_zoom1 <- ggplot() +
  geom_sf(data = birdregions, fill = "grey90", color = "black") +
  geom_sf(data = sites_sf_aves, aes(size = total_tested), 
          color = "#ab4901", fill = NA, lwd = 1.2, alpha = 0.4) +
  scale_size_continuous(name = "Total Tested", range = c(3, 20),
                        breaks = c(200, 2000, 20000, 40000),
                        labels = function(x) format(x, scientific = FALSE)) +
  labs(title = "Birds: Total Tested (Zoomed to Palearctic)",
       subtitle = "Sites (class = Aves) overlayed on bird regions") +
  coord_sf(xlim = c(-25, 70), ylim = c(20, 65)) +
  theme_classic()
print(p_zoom1)

p_zoom2 <- ggplot() +
  geom_sf(data = birdregions, fill = "grey90", color = "black") +
  geom_sf(data = sites_sf_aves, aes(size = total_tested), 
          color = "#ab4901", fill = NA, lwd = 1.2, alpha = 0.4) +
  scale_size_continuous(name = "Total Tested", range = c(3, 20),
                        breaks = c(200, 2000, 20000, 40000),
                        labels = function(x) format(x, scientific = FALSE)) +
  labs(title = "Birds: Total Tested (Zoomed to Nearctic)",
       subtitle = "Sites (class = Aves) overlayed on bird regions") +
  coord_sf(xlim = c(-150, -50), ylim = c(0, 60)) +
  theme_classic()
print(p_zoom2)

combined_plot3 <- p_zoom1 / p_zoom2
print(combined_plot3)

##----------------------------------------------------
# Plot the number of individuals tested by Class for each year
# Summarise total individuals tested per class
# Summarise total individuals tested per class
class_summary <- wnv_data_clean %>%
  group_by(class) %>%
  summarise(Total_Individuals_Tested = sum(total.tested, na.rm = TRUE)) %>%
  arrange(desc(Total_Individuals_Tested))

print(class_summary)

# Create custom color palette
class_colors <- c(
  "Aves" = "#ab4901",
  "Mammalia" = "#065471",
  "Amphibia" = "#cb4c78",
  "Sauropsida" = "#39b54a"
)

# Create cleaned dataset for plotting
class_data <- data.frame(
  class = c("Aves", "Mammalia", "Amphibia", "Sauropsida"),
  total_tested = c(431581, 102028, 96, 1101)
)

# Calculate percentages
class_data <- class_data %>%
  mutate(percent = 100 * total_tested / sum(total_tested),
         label = paste0(total_tested, " (", round(percent, 1), "%)"))

# Stacked Bar Chart (one bar, split by class)
# Plot
# Plot
library(scales) 
ggplot(class_data, aes(x = reorder(class, total_tested), y = total_tested, fill = class)) +
  geom_col(width = 0.6, alpha = 0.7) +
  geom_text(aes(label = label), hjust = -0.1, size = 4.5) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.15))) +  # no scientific, add margin
  scale_fill_manual(values = c("#cb4c78", "#ab4901", "#065471", "#39b54a")) +
  coord_flip() +
  theme_classic() +
  labs(title = "Total Number of Individuals Tested per Class (WNV Surveillance)",
       x = NULL, y = "Total Tested") +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 17),
    axis.text.y = element_text(size = 16),
    plot.title = element_text(size = 14, face = "bold")
  )

ggplot(class_data, aes(x = reorder(class, total_tested), y = total_tested, fill = class)) +
  geom_col(width = 0.6, alpha = 0.7) +
  geom_text(aes(label = label), hjust = -0.1, size = 3.5) +
  scale_fill_manual(values = c("#cb4c78", "#ab4901", "#065471", "#39b54a")) +
  coord_flip() +
  theme_classic() +
  labs(title = "Total Number of Individuals Tested per Class (WNV Surveillance)",
       x = NULL, y = "Total Tested") +
  theme(legend.position = "none")




# Summarize the total tested per year and class
wnv_summary <- wnv_data_clean %>%
  filter(!is.na(sampling.year), !is.na(class)) %>%
  group_by(sampling.year, class) %>%
  summarise(total_tested = sum(total.tested, na.rm = TRUE), .groups = 'drop')

# Define manual colors for classes
class_colors <- c(
  "Aves" = "#ab4901",        # red
  "Mammalia" = "#065471",    # blue
  "Amphibia" = "#cb4c78",    # purple
  "Sauropsida" = "#39b54a"   # green
)

# Define break points for every 10 years from 1950 to latest year
year_breaks <- seq(1950, max(wnv_summary$sampling.year, na.rm = TRUE), by = 10)

# Plot using ggplot
ggplot(wnv_summary, aes(x = sampling.year, y = total_tested, fill = class)) +
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
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12))

##----------------------------------------------------
# Plot the number of new species tested by year and the cummulative curve of positive species over time.
# Load necessary libraries
library(dplyr)
library(ggplot2)

# Step 1: Filter dataset for Aves and Mammalia & ensure valid sampling.year
wnv_filtered <- wnv_data_clean %>%
  filter(class %in% c("Aves", "Mammalia")) %>%   # Keep only Aves & Mammalia
  filter(!is.na(sampling.year)) %>%              # Remove missing years
  mutate(sampling.year = as.numeric(sampling.year)) %>%  # Ensure year is numeric
  filter(!is.na(sampling.year))  # Remove rows where conversion failed

# Step 2: Identify **first time a species was tested** (regardless of WNV status)
wnv_first_tested <- wnv_filtered %>%
  select(class, sampling.year, accepted_scientific_name) %>%
  group_by(class, accepted_scientific_name) %>%
  summarize(first_tested_year = min(sampling.year, na.rm = TRUE), .groups = "drop") # Get first year

# Step 3: Identify **first time a species tested positive for WNV** (ratio.m1 > 0)
wnv_first_positive <- wnv_filtered %>%
  filter(ratio.m1 > 0) %>%
  select(class, sampling.year, accepted_scientific_name) %>%
  group_by(class, accepted_scientific_name) %>%
  summarize(first_positive_year = min(sampling.year, na.rm = TRUE), .groups = "drop") # Get first positive year

# Step 4: Count cumulative species over time
cumulative_tested <- wnv_first_tested %>%
  count(class, first_tested_year) %>%
  arrange(class, first_tested_year) %>%
  group_by(class) %>%
  mutate(cumulative_species = cumsum(n)) %>%
  rename(sampling.year = first_tested_year) %>%
  mutate(status = "Tested")  # Add status for legend

cumulative_positive <- wnv_first_positive %>%
  count(class, first_positive_year) %>%
  arrange(class, first_positive_year) %>%
  group_by(class) %>%
  mutate(cumulative_species = cumsum(n)) %>%
  rename(sampling.year = first_positive_year) %>%
  mutate(status = "WNV-Positive")  # Add status for legend

# Step 5: Combine the two datasets
cumulative_combined <- bind_rows(cumulative_tested, cumulative_positive)

# Define break points for every 10 years from 1950 to latest year
year_breaks <- seq(1950, max(wnv_summary$sampling.year, na.rm = TRUE), by = 10)

# Step 6: Plot cumulative species curves for first-time tested and WNV-positive species
ggplot(cumulative_combined, aes(x = sampling.year, y = cumulative_species, color = class, linetype = status)) +
   geom_line(linewidth = 0.8) +  # Thinner lines
   geom_point(size = 2) +  # Adjust point size
   scale_color_manual(values = c("Aves" = "#ab4901", "Mammalia" = "#065471")) +  # Custom colors
   scale_linetype_manual(values = c("Tested" = "dotted", "WNV-Positive" = "solid")) +  # Different line types
   scale_x_continuous(breaks = year_breaks) +
   labs(title = "Cumulative First-Time Tested and WNV-Positive Species Over Time",
        subtitle = "Only counting first detection or test year for each species",
        x = "Sampling Year",
        y = "Cumulative Number of Species",
        color = "Class",
        linetype = "Species Status") +
   theme_minimal() +
   theme(
     text = element_text(size = 14),
     panel.grid = element_blank(),  # Remove all grid lines
     panel.border = element_blank(),  # Remove panel border
     axis.line = element_line(linewidth = 0.5),  # Keep only axis lines
     legend.position = "top"  # Move legend to top for better visibility
   )

#------------------------------------------------------------------------
# Load necessary libraries
library(dplyr)
library(ggplot2)
library(forcats)  # For manual factor ordering

# Step 1: Filter dataset for 'Aves' only and remove NA in Zoorealms
wnv_filtered <- wnv_data_clean %>% 
  filter(class == "Aves") %>%   # Keep only Aves
  filter(!is.na(sampling.year)) %>%  # Remove missing years
  mutate(sampling.year = as.numeric(sampling.year))  # Ensure year is numeric

# Step 2: Identify the first time a species was tested (regardless of WNV status)
wnv_first_tested <- wnv_filtered %>%
  select(class, sampling.year, accepted_scientific_name) %>%
  group_by(class, accepted_scientific_name) %>%
  summarize(first_tested_year = min(sampling.year, na.rm = TRUE), .groups = "drop")  # Get first test year

# Step 3: Identify the first time a species tested positive for WNV (ratio.m1 > 0)
wnv_first_positive <- wnv_filtered %>%
  filter(ratio.m1 > 0) %>%
  select(class, sampling.year, accepted_scientific_name) %>%
  group_by(class, accepted_scientific_name) %>%
  summarize(first_positive_year = min(sampling.year, na.rm = TRUE), .groups = "drop")  # Get first positive year

# Step 4: Count cumulative species over time (Grouped by year)
cumulative_tested <- wnv_first_tested %>%
  count(class, first_tested_year) %>%
  arrange(class, first_tested_year) %>%
  group_by(class) %>%
  mutate(cumulative_species = cumsum(n)) %>%
  rename(sampling.year = first_tested_year) %>%
  mutate(status = "Tested")  # Add status for legend

cumulative_positive <- wnv_first_positive %>%
  count(class, first_positive_year) %>%
  arrange(class, first_positive_year) %>%
  group_by(class) %>%
  mutate(cumulative_species = cumsum(n)) %>%
  rename(sampling.year = first_positive_year) %>%
  mutate(status = "WNV-Positive")  # Add status for legend

# Step 5: Combine the two datasets
cumulative_combined <- bind_rows(cumulative_tested, cumulative_positive)

# Step 6: Create the cumulative plot with bars for cumulative species tested and a line for cumulative positive species
p_aves <- ggplot(cumulative_combined, aes(x = sampling.year, y = cumulative_species, color = status)) +
  geom_bar(data = subset(cumulative_combined, status == "Tested"), 
           stat = "identity", position = "dodge", fill = "#ab4901", alpha = 0.5, color = NA) +  # Bars for tested species (no border)
  geom_line(data = subset(cumulative_combined, status == "WNV-Positive"), 
            aes(group = 1), size = 1, color = "#ab4901") +  # Line for positive species (fixed color code)
  geom_point(data = subset(cumulative_combined, status == "WNV-Positive"), 
             size = 3, color = "#ab4901") +  # Points for positive species
  labs(title = "Cumulative First-Time Tested and WNV-Positive Species Over Time",
       subtitle = "For Aves Class Only",
       x = "Sampling Year",
       y = "Cumulative Number of Species",
       color = "Species Status") +
  theme_minimal() +
  theme(
    text = element_text(size = 12),
    panel.grid = element_blank(),  # Remove grid lines
    panel.border = element_blank(),  # Remove panel border
    axis.line = element_line(size = 0.5),  # Keep only axis lines
    legend.position = "top"  # Position the legend at the top
  ) +
  scale_y_continuous(limits = c(0, 1400)) +  # Set Y-axis range
  scale_x_continuous(limits = c(1950, 2025), breaks = seq(1950, 2025, by = 10))  # Set X-axis range and breaks every 10 years

# Display the plot
print(p_aves)


## SAME FOR MAMMALS
# Step 1: Filter dataset for 'Mammals' only and remove NA in Zoorealms
wnv_filtered_m <- wnv_data_clean %>%  
  filter(class == "Mammalia") %>%   # Keep only Mammals
  filter(!is.na(sampling.year)) %>%  # Remove missing years
  mutate(sampling.year = as.numeric(sampling.year))  # Ensure year is numeric

# Step 2: Identify the first time a species was tested (regardless of WNV status)
wnv_first_tested_m <- wnv_filtered_m %>%
  select(class, sampling.year, accepted_scientific_name) %>%
  group_by(class, accepted_scientific_name) %>%
  summarize(first_tested_year = min(sampling.year, na.rm = TRUE), .groups = "drop")  # Get first test year

# Step 3: Identify the first time a species tested positive for WNV (ratio.m1 > 0)
wnv_first_positive_m <- wnv_filtered_m %>%
  filter(ratio.m1 > 0) %>%
  select(class, sampling.year, accepted_scientific_name) %>%
  group_by(class, accepted_scientific_name) %>%
  summarize(first_positive_year = min(sampling.year, na.rm = TRUE), .groups = "drop")  # Get first positive year

# Step 4: Count cumulative species over time (Grouped by year)
cumulative_tested_m <- wnv_first_tested_m %>%
  count(class, first_tested_year) %>%
  arrange(class, first_tested_year) %>%
  group_by(class) %>%
  mutate(cumulative_species = cumsum(n)) %>%
  rename(sampling.year = first_tested_year) %>%
  mutate(status = "Tested")  # Add status for legend

cumulative_positive_m <- wnv_first_positive_m %>%
  count(class, first_positive_year) %>%
  arrange(class, first_positive_year) %>%
  group_by(class) %>%
  mutate(cumulative_species = cumsum(n)) %>%
  rename(sampling.year = first_positive_year) %>%
  mutate(status = "WNV-Positive")  # Add status for legend

# Step 5: Combine the two datasets
cumulative_combined_m <- bind_rows(cumulative_tested_m, cumulative_positive_m)

# Step 6: Create the cumulative plot with bars for cumulative species tested and a line for cumulative positive species
p_mammals <- ggplot(cumulative_combined_m, aes(x = sampling.year, y = cumulative_species, color = status)) + 
  geom_bar(data = subset(cumulative_combined_m, status == "Tested"), 
           stat = "identity", position = "dodge", fill = "#065471", alpha = 0.5, color = NA) +  # Bars for tested species (no border) 
  geom_line(data = subset(cumulative_combined_m, status == "WNV-Positive"), 
            aes(group = 1), size = 1, color = "#065471") +  # Line for positive species (fixed color code)
  geom_point(data = subset(cumulative_combined_m, status == "WNV-Positive"), 
             size = 3, color = "#065471") +  # Points for positive species 
  labs(title = "Cumulative First-Time Tested and WNV-Positive Species Over Time",
       subtitle = "For Mammals Class Only",
       x = "Sampling Year",
       y = "Cumulative Number of Species",
       color = "Species Status") +
  theme_minimal() +
  theme(
    text = element_text(size = 12),
    panel.grid = element_blank(),  # Remove grid lines
    panel.border = element_blank(),  # Remove panel border
    axis.line = element_line(size = 0.5),  # Keep only axis lines
    legend.position = "top"  # Position the legend at the top
  ) +
  scale_y_continuous(limits = c(0, 700)) +  # Set Y-axis range
  scale_x_continuous(limits = c(1950, 2025), breaks = seq(1950, 2025, by = 10))  # Set X-axis range and breaks every 10 years

# Display the plot
print(p_mammals)


# Combine both plots (aves + mammals)
combined_plot_cumulative <- p_aves / p_mammals
print(combined_plot_cumulative)







# Step 1: Filter dataset for Aves and Mammalia & remove NA in Zoorealms
wnv_filtered <- wnv_data_clean %>%
  filter(class %in% c("Aves", "Mammalia")) %>%   # Keep only Aves & Mammalia
  filter(!is.na(sampling.year) & !is.na(zoorealms)) %>%  # Remove missing years & zoorealms
  mutate(sampling.year = as.numeric(sampling.year))  # Ensure year is numeric

# Step 2: Identify **first time a species was tested** (regardless of WNV status)
wnv_first_tested <- wnv_filtered %>%
  select(class, sampling.year, accepted_scientific_name, zoorealms) %>%
  group_by(class, accepted_scientific_name, zoorealms) %>%
  summarize(first_tested_year = min(sampling.year, na.rm = TRUE), .groups = "drop")  # Get first test year

# Step 3: Identify **first time a species tested positive for WNV** (ratio.m1 > 0)
wnv_first_positive <- wnv_filtered %>%
  filter(ratio.m1 > 0) %>%
  select(class, sampling.year, accepted_scientific_name, zoorealms) %>%
  group_by(class, accepted_scientific_name, zoorealms) %>%
  summarize(first_positive_year = min(sampling.year, na.rm = TRUE), .groups = "drop")  # Get first positive year

# Step 4: Count cumulative species over time **(Grouped by Zoorealm)**
cumulative_tested <- wnv_first_tested %>%
  count(class, zoorealms, first_tested_year) %>%
  arrange(class, zoorealms, first_tested_year) %>%
  group_by(class, zoorealms) %>%
  mutate(cumulative_species = cumsum(n)) %>%
  rename(sampling.year = first_tested_year) %>%
  mutate(status = "Tested")  # Add status for legend

cumulative_positive <- wnv_first_positive %>%
  count(class, zoorealms, first_positive_year) %>%
  arrange(class, zoorealms, first_positive_year) %>%
  group_by(class, zoorealms) %>%
  mutate(cumulative_species = cumsum(n)) %>%
  rename(sampling.year = first_positive_year) %>%
  mutate(status = "WNV-Positive")  # Add status for legend

# Step 5: Combine the two datasets
cumulative_combined <- bind_rows(cumulative_tested, cumulative_positive)

# Step 6: Manually order Zoorealms for facet arrangement
custom_zoorealm_order <- c(
  "Saharo-Arabian", "Afrotropical", "Palearctic", "Madagascan",  # Row 1
  "Nearctic", "Panamanian", "Neotropical",  # Row 2
  "Oriental", "Sino-Japanese", "Australian"  # Row 3
)

cumulative_combined <- cumulative_combined %>%
  mutate(zoorealms = factor(zoorealms, levels = custom_zoorealm_order)) %>%
  filter(!is.na(zoorealms))  # Ensure NA zoorealms are removed

# Step 7: Create the cumulative plot, faceting by Zoorealm in 3 rows
p <- ggplot(cumulative_combined, aes(x = sampling.year, y = cumulative_species, 
                                     color = class, linetype = status)) +
  geom_line(size = 0.8) +  # Thinner lines
  geom_point(size = 1.5) + # Adjust point size
  scale_color_manual(values = c("Aves" = "#e31a1c", "Mammalia" = "#1f78b4")) +  # Custom colors
  scale_linetype_manual(values = c("Tested" = "dotted", "WNV-Positive" = "solid")) +  # Different line types
  labs(title = "Cumulative First-Time Tested and WNV-Positive Species Over Time",
       subtitle = "Grouped by Zoorealm (Only first detection/test year counted)",
       x = "Sampling Year",
       y = "Cumulative Number of Species",
       color = "Class",
       linetype = "Species Status") +
  theme_minimal() +
  theme(
    text = element_text(size = 12),
    panel.grid = element_blank(),  # Remove all grid lines
    panel.border = element_blank(),  # Remove panel border
    axis.line = element_line(size = 0.5),  # Keep only axis lines
    legend.position = "top"  # Move legend to top for better visibility
  ) +
  facet_wrap(~zoorealms, nrow = 3, scales = "free_y")  # Facet by zoorealms with 3 rows

# Display the plot
print(p)





#### -------- CALCULATE HOW MANY SPECIES HAVE BEEN TESTED WITH DIFFERENT METHODS
## AND HOW MANY SPECIES HAVE BEEN DETECTED WITH THE DIFFRERENT METHODS

# Clean and ensure numeric conversion
wnv_data_clean <- wnv_data_clean %>%
  filter(!is.na(accepted_scientific_name), accepted_scientific_name != "NA") %>%
  mutate(
    total.tested = as.numeric(total.tested),
    positive.individuals.m1 = as.numeric(positive.individuals.m1),
    positive.individuals.m2 = as.numeric(positive.individuals.m2)
  )

# Pivot to long format to treat method.1 and method.2 separately
long_df <- wnv_data_clean %>%
  pivot_longer(cols = c(method.1, method.2), names_to = "method_type", values_to = "method") %>%
  mutate(
    positive = ifelse(method_type == "method.1", positive.individuals.m1, positive.individuals.m2)
  ) %>%
  filter(!is.na(method), method != "NA")

# Summarise per species and method
summary_species_method <- long_df %>%
  group_by(accepted_scientific_name, family, order, class, method) %>%
  summarise(
    total_tested = sum(total.tested, na.rm = TRUE),
    total_positive = sum(positive, na.rm = TRUE),
    .groups = "drop"
  )

# View result
head(summary_species_method)

# Summarize WNV testing results by method and class
summary_by_method_class <- summary_species_method %>%
  group_by(class, method) %>%
  summarise(
    total_species_tested = n_distinct(accepted_scientific_name),
    total_individuals_tested = sum(total_tested, na.rm = TRUE),
    total_species_positive = n_distinct(accepted_scientific_name[total_positive > 0]),
    .groups = "drop"
  ) %>%
  arrange(class, method)

# View the result
print(summary_by_method_class)

# Create a dataframe for plotting
plot_data <- summary_by_method_class %>%
  group_by(class) %>%
  summarise(
    tested = sum(total_species_tested),
    positive = sum(total_species_positive)
  ) %>%
  pivot_longer(cols = c(tested, positive), names_to = "status", values_to = "n_species") %>%
  mutate(status = factor(status, levels = c("tested", "positive")))

# Define custom colors
status_colors <- c("tested" = "#a6cee3", "positive" = "#1f78b4")

# Plot: stacked bar chart
ggplot(plot_data, aes(x = class, y = n_species, fill = status)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = status_colors,
                    labels = c("Total Tested", "Positive Species")) +
  labs(
    title = "WNV Testing Outcomes by Class",
    x = "Taxonomic Class",
    y = "Number of Species",
    fill = "Testing Outcome"
  ) +
  theme_classic(base_size = 13)

#-------------Plot number of species tested & postive with each method

# Reshape data to long format
plot_data <- summary_by_method_class %>%
  pivot_longer(cols = c(total_species_tested, total_species_positive),
               names_to = "status", values_to = "n_species") %>%
  mutate(status = recode(status,
                         total_species_tested = "tested",
                         total_species_positive = "positive"),
         status_class = paste(class, status, sep = "_"))

# Custom colors
custom_colors <- c(
  "Aves_tested" = "#ab9601", "Aves_positive" = "#ab4901",
  "Mammalia_tested" = "#069071", "Mammalia_positive" = "#065471",
  "Amphibia_tested" = "#cba078", "Amphibia_positive" = "#cb4c78",
  "Sauropsida_tested" = "#8cc63f", "Sauropsida_positive" = "#39b54a"
)

# Plot
ggplot(plot_data, aes(x = method, y = n_species, fill = status_class)) +
  geom_col(position = "stack") +
  scale_fill_manual(values = custom_colors) +
  labs(title = "Number of Species Tested and Positive for WNV by Method and Class",
       x = "Diagnostic Method",
       y = "Number of Species",
       fill = "Class / Status") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#-------------

# Create dataset with negative species count
summary_with_negatives <- summary_by_method_class %>%
  filter(class %in% c("Aves", "Mammalia")) %>%
  mutate(total_species_negative = total_species_tested - total_species_positive)

# Reshape to long format for plotting
long_data <- summary_with_negatives %>%
  select(method, class, total_species_positive, total_species_negative) %>%
  pivot_longer(cols = c(total_species_positive, total_species_negative),
               names_to = "status", values_to = "count") %>%
  complete(method = unique(summary_by_method_class$method),
           class = c("Aves", "Mammalia"),
           status = c("total_species_positive", "total_species_negative"),
           fill = list(count = 0))

# Order factors
long_data$method <- factor(long_data$method, levels = unique(summary_by_method_class$method))
long_data$status <- factor(long_data$status,
                           levels = c("total_species_negative", "total_species_positive"))

# Custom colors for negative and positive species
fill_colors <- c("total_species_negative" = "#1f78b4",  # light grey
                 "total_species_positive" = "#ab0e01")  # blue

# Plot for Aves
plot_aves <- ggplot(long_data %>% filter(class == "Aves"),
                    aes(x = method, y = count, fill = status)) +
  geom_col(position = "stack") +
  scale_fill_manual(values = fill_colors,
                    labels = c("Negative species", "Positive species")) +
  labs(title = "Bird Species Positive and Negative by Method",
       x = "", y = "Number of Species", fill = "") +
  theme_classic(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Plot for Mammalia
plot_mammals <- ggplot(long_data %>% filter(class == "Mammalia"),
                       aes(x = method, y = count, fill = status)) +
  geom_col(position = "stack") +
  scale_fill_manual(values = fill_colors,
                    labels = c("Negative species", "Positive species")) +
  labs(title = "Mammal Species Positive and Negative by Method",
       x = "Diagnostic Method", y = "Number of Species", fill = "") +
  theme_classic(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Combine both plots
plot_aves / plot_mammals + plot_layout(ncol = 1)

###-----------------------------------------------------------------------
library(tidyverse)
library(patchwork)

# Define serological group
serological_methods <- c("HI", "IgC", "IgM", "MIA", "serological", "serological.ELISA")

# Recode methods into grouped categories
summary_grouped <- summary_by_method_class %>%
  filter(class %in% c("Aves", "Mammalia")) %>%
  mutate(method_grouped = case_when(
    method %in% serological_methods ~ "serological",
    method == "serological.VNT" ~ "serological.VNT",
    method == "molecular_detection" ~ "molecular_detection",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(method_grouped)) %>%
  group_by(class, method = method_grouped) %>%
  summarise(
    total_species_tested = sum(total_species_tested, na.rm = TRUE),
    total_species_positive = sum(total_species_positive, na.rm = TRUE),
    total_species_negative = total_species_tested - total_species_positive,
    .groups = "drop"
  )

# Reshape for plotting
long_data <- summary_grouped %>%
  select(method, class, total_species_positive, total_species_negative) %>%
  pivot_longer(cols = c("total_species_positive", "total_species_negative"),
               names_to = "status", values_to = "count") %>%
  complete(method = c("molecular_detection", "serological.VNT", "serological"),
           class = c("Aves", "Mammalia"),
           status = c("total_species_positive", "total_species_negative"),
           fill = list(count = 0))

# Set factor order
long_data$method <- factor(long_data$method,
                           levels = c("molecular_detection", "serological.VNT", "serological"))
long_data$status <- factor(long_data$status,
                           levels = c("total_species_negative", "total_species_positive"))

# Custom fill colors
fill_colors <- c("total_species_negative" = "#1f78b4",  # blue
                 "total_species_positive" = "#ab0e01")  # red

# Plot for Aves
plot_aves <- ggplot(long_data %>% filter(class == "Aves"),
                    aes(x = method, y = count, fill = status)) +
  geom_col(position = "stack", width = 0.6) +
  scale_fill_manual(values = fill_colors,
                    labels = c("Negative species", "Positive species")) +
  labs(title = "Bird Species Tested for WNV by Method",
       x = "", y = "Number of Species", fill = "") +
  theme_classic(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Plot for Mammalia
plot_mammals <- ggplot(long_data %>% filter(class == "Mammalia"),
                       aes(x = method, y = count, fill = status)) +
  geom_col(position = "stack", width = 0.6) +
  scale_fill_manual(values = fill_colors,
                    labels = c("Negative species", "Positive species")) +
  labs(title = "Mammal Species Tested for WNV by Method",
       x = "Diagnostic Method", y = "Number of Species", fill = "") +
  theme_classic(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Combine plots
plot_aves / plot_mammals + plot_layout(ncol = 1)


##---------------------------------------------------------
# Filter for birds and selected methods
methods_of_interest <- c("molecular_detection", "serological.VNT", "serological")

bird_subset <- summary_species_method %>%
  filter(class == "Aves", method %in% methods_of_interest) %>%
  mutate(present = total_tested > 0)

# Create presence/absence matrix
presence_matrix <- bird_subset %>%
  select(accepted_scientific_name, method, present) %>%
  pivot_wider(names_from = method, values_from = present, values_fill = FALSE)

# Count how many species fall into each method combination
summary_combinations <- presence_matrix %>%
  group_by(across(all_of(methods_of_interest))) %>%
  summarise(n_species = n(), .groups = "drop")

print(summary_combinations)

# Define methods of interest
methods_of_interest <- c("molecular_detection", "serological.VNT", "serological")

# Prepare binary presence matrix
presence_matrix <- summary_species_method %>%
  filter(class == "Aves", method %in% methods_of_interest) %>%
  mutate(present = total_tested > 0) %>%
  select(accepted_scientific_name, method, present) %>%
  pivot_wider(names_from = method, values_from = present, values_fill = FALSE) %>%
  mutate(across(all_of(methods_of_interest), as.integer))

# Generate UpSet plot
upset(
  presence_matrix,
  intersect = methods_of_interest,
  name = "Tested by Method",
  base_annotations = list(
    'Species Count' = intersection_size()
  ),
  width_ratio = 0.1
)




#--------------
# Step 1: Clean and prepare data
methods_of_interest <- c("molecular_detection", "serological.VNT", 
                          "serological.ELISA")

wnv_data_clean <- wnv_data_clean %>%
  filter(!is.na(accepted_scientific_name), accepted_scientific_name != "NA") %>%
  mutate(
    total.tested = as.numeric(total.tested),
    positive.individuals.m1 = as.numeric(positive.individuals.m1),
    positive.individuals.m2 = as.numeric(positive.individuals.m2)
  )

# Step 2: Long format for both methods
long_df <- wnv_data_clean %>%
  pivot_longer(cols = c(method.1, method.2),
               names_to = "method_type", values_to = "method") %>%
  mutate(
    positive = ifelse(method_type == "method.1", positive.individuals.m1, positive.individuals.m2)
  ) %>%
  filter(!is.na(method), method %in% methods_of_interest, class == "Aves")

# Step 3: Summarize tested & positive per species × method
summary_matrix <- long_df %>%
  group_by(accepted_scientific_name, method) %>%
  summarise(
    total_tested = sum(total.tested, na.rm = TRUE),
    total_positive = sum(positive, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    tested_flag = total_tested > 0,
    positive_flag = total_positive > 0 & tested_flag  # Positive only if also tested
  )

# Step 4: Pivot to tested/positive presence matrices
get_presence_matrix <- function(df, flag_col) {
  df %>%
    filter(.data[[flag_col]]) %>%
    mutate(present = 1L) %>%
    select(accepted_scientific_name, method, present) %>%
    pivot_wider(names_from = method, values_from = present, values_fill = list(present = 0)) %>%
    arrange(accepted_scientific_name)
}

tested_presence <- get_presence_matrix(summary_matrix, "tested_flag")
positive_presence <- get_presence_matrix(summary_matrix, "positive_flag")

# Step 5: Create combo labels
combine_methods <- function(df) {
  df %>%
    mutate(combo = apply(select(., -accepted_scientific_name), 1, paste, collapse = "_")) %>%
    count(combo, name = "n")
}

tested_combos   <- combine_methods(tested_presence)
positive_combos <- combine_methods(positive_presence)

# Step 6: Compare combo counts
comparison <- full_join(tested_combos, positive_combos, by = "combo", suffix = c("_tested", "_positive")) %>%
  replace_na(list(n_tested = 0, n_positive = 0)) %>%
  mutate(diff = n_positive - n_tested)

# Step 7: Sanity checks
cat("\n--- Sanity Checks ---\n")
cat("Species tested per method:\n")
print(colSums(select(tested_presence, -accepted_scientific_name)))

cat("\nSpecies tested by how many methods:\n")
tested_presence %>%
  mutate(n_methods = rowSums(select(., -accepted_scientific_name))) %>%
  count(n_methods, name = "n_species") %>%
  print()

cat("\n✅ Final Comparison Table (no more false positives!):\n")
print(comparison %>% arrange(desc(n_tested)), n = 50)

# Define method order
method_order <- c("molecular_detection", "serological.VNT", "serological.ELISA")

# Reorder the columns if necessary
tested_presence <- tested_presence %>%
  select(accepted_scientific_name, all_of(method_order))

# Make sure each column is integer
tested_presence[method_order] <- lapply(tested_presence[method_order], as.integer)

# Tested species combinations
upset(
  tested_presence,
  intersect = method_order,
  name = "WNV Detection Method",
  base_annotations = list(
    'Species count' = intersection_size()
  ),
  width_ratio = 0.1
) +
  theme_classic(base_size = 16) +
  labs(title = "Species Tested for WNV by Method Combination")

upset(
  positive_presence,
  intersect = methods_of_interest,
  name = "WNV Detection Method",
  base_annotations = list(
    'Species count' = intersection_size()
  ),
  width_ratio = 0.1
) +
  theme_minimal(base_size = 14) +
  labs(title = "Species Positive for WNV by Method Combination")


#----- Species & individuals tested with different method by countries

# Clean and ensure numeric conversion
wnv_data_clean <- wnv_data_clean %>%
  filter(!is.na(accepted_scientific_name), accepted_scientific_name != "NA") %>%
  mutate(
    total.tested = as.numeric(total.tested),
    positive.individuals.m1 = as.numeric(positive.individuals.m1),
    positive.individuals.m2 = as.numeric(positive.individuals.m2)
  )

# Pivot to long format to treat method.1 and method.2 separately
long_df <- wnv_data_clean %>%
  pivot_longer(cols = c(method.1, method.2), names_to = "method_type", values_to = "method") %>%
  mutate(
    positive = ifelse(method_type == "method.1", positive.individuals.m1, positive.individuals.m2)
  ) %>%
  filter(!is.na(method), method != "NA")

# Summarize per species and method
summary_species_method <- long_df %>%
  group_by(accepted_scientific_name, family, order, class, country, method) %>%
  summarise(
    total_tested = sum(total.tested, na.rm = TRUE),
    total_positive = sum(positive, na.rm = TRUE),
    .groups = "drop"
  )

# Now, calculate the total species tested and total individuals tested by class and method, for Aves and Mammals separately.

# For class Aves
summary_aves <- summary_species_method %>%
  filter(class == "Aves") %>%
  group_by(country, method) %>%
  summarise(
    num_species = n_distinct(accepted_scientific_name),
    total_individuals_tested = sum(total_tested, na.rm = TRUE),
    .groups = "drop"
  )

# For class Mammalia
summary_mammals <- summary_species_method %>%
  filter(class == "Mammalia") %>%
  group_by(country, method) %>%
  summarise(
    num_species = n_distinct(accepted_scientific_name),
    total_individuals_tested = sum(total_tested, na.rm = TRUE),
    .groups = "drop"
  )

# View summary for Aves
print(summary_aves)

# View summary for Mammals
print(summary_mammals)

# Write the results to an Excel file for further analysis
write_xlsx(list(Aves = summary_aves, Mammals = summary_mammals), "species_testing_summary_by_method.xlsx")

library(gridExtra)

# Reorganize 'method' column: group specific methods into 'others'
summary_species_country <- summary_aves %>%
  mutate(method = case_when(
    method %in% c("molecular_detection", "serological.ELISA", "serological.VNT") ~ method,  # Keep these methods as is
    method %in% c("HI", "IgC", "IgM", "IHC", "MIA", "mortality", "serological") ~ "others",  # Group these into 'others'
    TRUE ~ method  # Leave any other methods unchanged (if present)
  ))

# Summarize the number of species tested by country and method
summary_species_country_method <- summary_species_country %>%
  group_by(country, method) %>%
  summarise(num_species = sum(num_species, na.rm = TRUE), .groups = "drop")

# Calculate the total number of species tested per country
country_totals <- summary_species_country_method %>%
  group_by(country) %>%
  summarise(total_species = sum(num_species, na.rm = TRUE)) %>%
  arrange(total_species)  # Arrange countries by total species tested, descending order

# Reorder countries based on total species tested
summary_species_country_method$country <- factor(summary_species_country_method$country, 
                                                 levels = country_totals$country)

# Plot the data: Countries on the Y-axis, stacked bars for methods
ggplot(summary_species_country_method, aes(x = num_species, y = country, fill = method)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +  # Stacked bars
  labs(
    title = "Number of Species Tested by Country and Method (Aves)",
    x = "Number of Species Tested",
    y = "Country",
    fill = "Method"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 10),  # Adjust the size of X-axis text
    axis.text.y = element_text(size = 8),  # Adjust the size of Y-axis text for readability
    plot.title = element_text(hjust = 0.5),
    legend.position = "right"
  )


# Reorganize 'method' column: group specific methods into 'others'
summary_aves_individuals <- summary_aves %>%
  mutate(method = case_when(
    method %in% c("molecular_detection", "serological.ELISA", "serological.VNT") ~ method,  # Keep these methods as is
    method %in% c("HI", "IgC", "IgM", "IHC", "MIA", "mortality", "serological") ~ "others",  # Group these into 'others'
    TRUE ~ method  # Leave any other methods unchanged (if present)
  ))

# Summarize the total number of individuals tested by country and method
summary_species_country_method_individuals <- summary_aves_individuals %>%
  group_by(country, method) %>%
  summarise(total_individuals = sum(total_individuals_tested, na.rm = TRUE), .groups = "drop")

# Calculate the total number of individuals tested per country
country_totals_individuals <- summary_species_country_method_individuals %>%
  group_by(country) %>%
  summarise(total_individuals_tested = sum(total_individuals, na.rm = TRUE)) %>%
  arrange(total_individuals_tested)  # Arrange countries by total individuals tested, ascending order

# Reorder countries based on total individuals tested
summary_species_country_method_individuals$country <- factor(summary_species_country_method_individuals$country, 
                                                             levels = country_totals_individuals$country)

# Plot the data: Countries on the Y-axis, stacked bars for methods
ggplot(summary_species_country_method_individuals, aes(x = total_individuals, y = country, fill = method)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +  # Stacked bars
  labs(
    title = "Total Individuals Tested by Country and Method (Aves)",
    x = "Total Individuals Tested",
    y = "Country",
    fill = "Method"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 10),  # Adjust the size of X-axis text
    axis.text.y = element_text(size = 8),  # Adjust the size of Y-axis text for readability
    plot.title = element_text(hjust = 0.5),
    legend.position = "right"
  )

