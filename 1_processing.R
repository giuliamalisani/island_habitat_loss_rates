# =============================================================================
#Content: all habitat-loss (HL) and habitat-loss-rate (HLr) calculations and calibration of global and local HumHL. 
#See 2_visualization.R for stats and figures.

# Steps in this workflow:
#1. Intact habitat extraction (Kennedy et al. HMc raster -> per-island area) + island_information
#2. Global HumHL & Global HumHLr
#3. Local HumHL & Local HumHLr
#4. SeaHL & SeaHLr
#5. VolHL & VolHLr
#6. Calibration (global vs. local HumHL, per threshold)
#7. Join all HL sources + HumHLr/BaseHLr ratios
# =============================================================================

library(here)
library(readxl)
library(readr)
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(writexl)
library(sf)
library(terra)
library(exactextractr)

input_data_dir <- here("input_data")
output_data_dir <- here("output_data")

arrival_years <- read_xlsx(file.path(input_data_dir, "arrival_years.xlsx")) %>%
  distinct(island, archipelago, arrival_year)

#correct island and archipelago spelling
name_lookup <- read_xlsx(file.path(input_data_dir, "island_name_lookup.xlsx"))

harmonise <- function(df) {
  df %>%
    left_join(name_lookup, by = c("island" = "raw_name")) %>%
    mutate(island = coalesce(correct_name, island)) %>%
    select(-correct_name)
}

#----1.INTACT HABITAT EXTRACTION----
hmc <- rast(file.path(input_data_dir, "Kennedy_HMc.tif"))
island_boundaries_folder <- file.path(input_data_dir, "island_boundaries")
island_boundaries <- list.files(island_boundaries_folder, pattern = "\\.gpkg$", full.names = TRUE)

natural_habitat <- data.frame()
island_areas <- data.frame()

for (island_boundaries in island_boundaries) {
  island_name <- tools::file_path_sans_ext(basename(island_boundaries))
  island_boundary <- st_read(island_boundaries)
  
  #extract values
  hmc_values_list <- exact_extract(hmc, island_boundary)
  hmc_values <- do.call(rbind, hmc_values_list)
  hmc_values <- hmc_values[!is.na(hmc_values$value), ]
  
  #apply threshold
  area_threshold_0 <- sum(hmc_values$coverage_fraction[hmc_values$value <= 0.00])
  area_threshold_0.1 <- sum(hmc_values$coverage_fraction[hmc_values$value <= 0.10])
  area_threshold_0.15 <- sum(hmc_values$coverage_fraction[hmc_values$value <= 0.15])
  area_threshold_0.2 <- sum(hmc_values$coverage_fraction[hmc_values$value <= 0.20])
  area_threshold_0.3 <- sum(hmc_values$coverage_fraction[hmc_values$value <= 0.30])
  area_threshold_0.4 <- sum(hmc_values$coverage_fraction[hmc_values$value <= 0.40])
  
  #total island area
  total_area <- as.numeric(st_area(island_boundary)) %>% sum() / 1e6
  
  natural_habitat <- rbind(natural_habitat, data.frame(
    island = island_name,
    area_threshold_0 = area_threshold_0,
    area_threshold_0.10 = area_threshold_0.1,
    area_threshold_0.15 = area_threshold_0.15,
    area_threshold_0.20 = area_threshold_0.2,
    area_threshold_0.30 = area_threshold_0.3,
    area_threshold_0.40 = area_threshold_0.4))
  
  island_areas <- rbind(island_areas, data.frame(
    island = island_name,
    total_area_km2 = total_area))
}

#natural habitat
natural_habitat_long <- natural_habitat %>%
  pivot_longer(
    cols = starts_with("area_threshold"),
    names_to = "threshold",
    values_to = "hmc_natural_habitat_km2") %>%
  mutate(threshold = gsub("area_threshold_", "", threshold)) %>%
  harmonise()

write_xlsx(natural_habitat_long, file.path(output_data_dir, "global_humhl", "hmc_natural_habitat.xlsx"))

#island_information: names + present-day area + arrival data
island_information <- arrival_years %>%
  left_join(island_areas %>% harmonise(), by = "island")

write_xlsx(island_information, file.path(output_data_dir, "island_information.xlsx"))

#----2. Global HumHL & HumHLr----
Ghumhl_data <- natural_habitat_long %>%
  left_join(island_information, by = "island") %>%
  mutate(
    years_since_arrival = 2014 - arrival_year,
    global_humHL = total_area_km2 - hmc_natural_habitat_km2,
    global_humHLr = (total_area_km2 - hmc_natural_habitat_km2) / years_since_arrival)

write_xlsx(Ghumhl_data, file.path(output_data_dir, "global_humhl", "Ghumhl(R).xlsx"))

#----3. Local HumHL & Localc HumHLr----
Lhumhl_data <- read_xlsx(file.path(input_data_dir, "local_natural_covers.xlsx")) %>%
  select(-archipelago) %>%
  harmonise() %>%
  left_join(island_information, by = "island") %>%
  mutate(
    natural_habitat_km2 = coalesce(
      natural_cover_km2,
      natural_cover_pctg * total_area_km2),
    local_humHL = total_area_km2 - natural_habitat_km2,
    local_years_since_arrival = year_dataset - arrival_year,
    local_humHLr = local_humHL / local_years_since_arrival)

write_xlsx(Lhumhl_rates, file.path(output_data_dir, "Lhumhl(R).xlsx"))

#----4. SeaHL & SeaHLr----
parent_dir <- file.path(input_data_dir, "lgm_areas")
lgm_col <- "BP0021000"
lgm_years <- 21000

read_island <- function(folder) {
  island_name <- basename(folder)
  
  df <- read_delim(file.path(folder, "recarea.csv"), delim = ";",
                   show_col_types = FALSE, name_repair = "minimal") %>%
    rename(row_idx = 1)
  
  island_row <- df %>% filter(str_to_lower(recname) == str_to_lower(island_name))
  
  tibble(island = island_name, area_lgm_m2 = as.numeric(island_row[[lgm_col]]))
}

seahl_data <- map_dfr(list.dirs(parent_dir, recursive = FALSE), read_island) %>%
  harmonise() %>%
  left_join(island_information, by = "island") %>%
  mutate(
    area_lgm_km2 = area_lgm_m2 / 1e6,
    seaHL = area_lgm_km2 - total_area_km2,
    seaHLr = seaHL / lgm_years) %>%
  select(island, area_present_km2 = total_area_km2, area_lgm_km2, seaHL, seaHLr)

write_xlsx(seahl_data, file.path(output_data_dir, "seahl(R).xlsx"))

#----4. VolHL & VolHLr----
volcov <- read_xlsx(file.path(input_data_dir, "volcanic_covers.xlsx"))

volhl_data <- volcov %>%
  harmonise() %>%
  select(island, vol_area_pctg) %>%
  left_join(island_information, by = "island") %>%
  mutate(
    vol_area_pctg = as.numeric(vol_area_pctg),
    volHL = vol_area_pctg / 100 * total_area_km2,
    volHLr = volHL / 11700) %>%
  select(island, archipelago, total_area_km2, vol_area_pctg, volHL, volHLr)

write_xlsx(volhl_data, file.path(output_data_dir, "volhl(R).xlsx"))

#----5. Calibration — global (per threshold) vs. local HumHL----
cal_data <- Ghumhl_data %>%
  mutate(threshold = as.numeric(threshold)) %>%
  select(island, threshold, global_humHL) %>%
  left_join(Lhumhl_rates %>% select(island, local_humHL), by = "island") %>%
  mutate(diff = global_humHL - local_humHL)

calibration_stats <- cal_data %>%
  group_by(threshold) %>%
  summarise(
    MAE = mean(abs(diff), na.rm = TRUE),
    RMSE = sqrt(mean(diff^2, na.rm = TRUE)),
    MedAE = median(abs(diff), na.rm = TRUE),
    Median_diff = median(diff, na.rm = TRUE),
    .groups = "drop")

write_xlsx(calibration_stats, file.path(output_data_dir, "calibration_stats.xlsx"))

#select threshold (after reviewing calibration_stats)
selected_threshold <- 0.20

#----7. Join all HL sources + HumHLr/BaseHLr ratios (selected threshold only)----
hl_joined <- Ghumhl_data %>%
  mutate(threshold = as.numeric(threshold)) %>%
  filter(threshold == selected_threshold) %>%
  select(island, archipelago, arrival_year, years_since_arrival, total_area_km2,
         global_humHL, global_humHLr) %>%
  left_join(Lhumhl_data %>% select(island, local_humHL, local_humHLr), by = "island") %>%
  left_join(seahl_data %>% select(island, seaHL, seaHLr), by = "island") %>%
  left_join(volhl_data %>% select(island, volHL, volHLr), by = "island") %>%
  mutate(
    baseHLr = seaHLr + volHLr,
    global_ratio = if_else(baseHLr > 0, global_humHLr / baseHLr, NA_real_),
    local_ratio = if_else(baseHLr > 0, local_humHLr  / baseHLr, NA_real_))

write_xlsx(hl_joined, file.path(output_data_dir, "habitat_losses(R).xlsx"))
