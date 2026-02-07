
# LST analysis

# ============================================================
# PART 0: LIBRARIES AND PATH
# ============================================================

# ============================================================
# Load Required Libraries 
# ============================================================
library(terra)       # raster handling
library(sf)          # spatial vector data
library(sp)          # spatial objects for Gi*
library(spdep)       # spatial autocorrelation / Gi*
library(ggplot2)     # plotting
library(tmap)        # thematic mapping (optional)
library(dplyr)       # data manipulation
library(tidyr)       # data reshaping
library(lubridate)   # date handling
library(zoo)         # rolling mean / smoothing
library(data.table)  # run-length encoding
library(scales)      # formatting scales
library(e1071)       # skewness / stats
library(tidyverse)   # includes ggplot2, dplyr, tidyr (optional)

# ============================================================
# PART 1: Multi-City Multi-Year Spatial Distribution of LST (2021–2024)
# ============================================================

# Define cities, years, and base path
cities <- c("Addis", "Adama", "Jimma", "Harar")
years  <- 2021:2024
base_path <- "/Users/drdesalewmoges/Documents/Papers/RS_Heat_Exposure/Data/LST/"


# Read and combine rasters (all cities & years)
df_all <- expand.grid(City = cities, Year = years) %>%
  mutate(
    file_path = paste0(base_path, City, "_LST_", Year, ".tif")
  ) %>%
  rowwise() %>%
  mutate(
    data = list(
      rast(file_path) %>%
        as.data.frame(xy = TRUE, na.rm = TRUE) %>%
        rename(LST = 3) %>%
        mutate(City = City, Year = Year)
    )
  ) %>%
  ungroup() %>%
  select(data) %>%
  bind_rows()


# Custom LST palette
palette_lst <- c(
  "#235cb1", "#2171b5", "#6baed6", "#86e26f",
  "#fff705", "#ffb613", "#ff500d", "#CD2626"
)


# Multi-city × multi-year plot
p <- ggplot(df_all, aes(x = x, y = y, fill = LST)) +
  geom_raster() +
  facet_grid(City ~ Year) +
  scale_fill_gradientn(colours = palette_lst, name = "LST (°C)") +
  coord_equal() +
  theme_bw() +
  theme(
    strip.text = element_text(face = "bold"),
    axis.text  = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  labs(title = "Land Surface Temperature (2021–2024)")

print(p)


# ============================================================
# PART 2: Daily Land Surface Temperature Trends in Selected Ethiopian Cities (2021–2024)
# ============================================================

# Define cities and base path
cities <- c("Addis", "Adama", "Jimma", "Harar")
base_path <- "/Users/drdesalewmoges/Documents/Papers/RS_Heat_Exposure/Data/LST/"

# Read and preprocess LST (all cities)
lst <- lapply(cities, function(cty) {
  read.csv(paste0(base_path, cty, "_Daily_LST_2021-2024.csv")) %>%
    mutate(
      City = cty,
      date = as.Date(date)
    )
}) %>%
  bind_rows() %>%
  arrange(City, date) %>%
  group_by(City) %>%
  mutate(
    LST_Day   = na.locf(LST_Day, na.rm = FALSE),
    LST_Night = na.locf(LST_Night, na.rm = FALSE),
    Year      = year(date),
    Month     = month(date, label = TRUE, abbr = TRUE)
  ) %>%
  ungroup() %>%
  filter(LST_Night >= 0)

# Build compressed seasonal index (per city)
lst_season <- lst %>%
  group_by(City, Year) %>%
  mutate(season_day = row_number()) %>%
  ungroup() %>%
  group_by(City) %>%
  mutate(
    season_index = season_day + (Year - min(Year)) * max(season_day)
  ) %>%
  ungroup()

# Compute tick positions (month labels per city)
ticks <- lst_season %>%
  group_by(City, Year, Month) %>%
  summarize(
    breaks = mean(season_index),
    labels = as.character(Month),
    .groups = "drop"
  )

# Plot: smooth LST time series by city
ts_plot <- ggplot(lst_season, aes(x = season_index)) +
  geom_smooth(aes(y = LST_Day),
              method = "loess", span = 0.01,
              size = 0.5, se = FALSE, color = "#CD2626") +
  geom_smooth(aes(y = LST_Night),
              method = "loess", span = 0.01,
              size = 0.5, se = FALSE, color = "#EE9A00") +
  facet_wrap(~ City, ncol = 1, scales = "free_x") +
  scale_x_continuous(
    breaks = ticks$breaks,
    labels = ticks$labels
  ) +
  labs(
    x = "Month",
    y = "LST (°C)"
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(size = 10, face = "bold"),
    axis.text.y = element_text(size = 10, face = "bold"),
    axis.title  = element_text(size = 10, face = "bold"),
    legend.position = "bottom"
  )

ts_plot

# ============================================================
# PART 3: Seasonal Distribution of Daytime and Nighttime LST Across Ethiopian Cities (2021–2024)
# ============================================================

# Define cities and base path
cities <- c("Addis", "Adama", "Jimma", "Harar")
base_path <- "/Users/drdesalewmoges/Documents/Papers/RS_Heat_Exposure/Data/LST/"

# Read and preprocess LST (all cities)
lst <- lapply(cities, function(cty) {
  read.csv(paste0(base_path, cty, "_Daily_LST_2021-2024.csv")) %>%
    mutate(
      City = cty,
      date = as.Date(date)
    )
}) %>%
  bind_rows() %>%
  arrange(City, date) %>%
  group_by(City) %>%
  mutate(
    LST_Day   = na.locf(LST_Day, na.rm = FALSE),
    LST_Night = na.locf(LST_Night, na.rm = FALSE),
    Year      = year(date),
    Month     = month(date, label = TRUE, abbr = TRUE)
  ) %>%
  ungroup() %>%
  filter(LST_Night >= 0)

# Pivot longer (Day / Night)
lst_long <- lst %>%
  pivot_longer(
    cols = c(LST_Day, LST_Night),
    names_to  = "Time_of_Day",
    values_to = "Temperature"
  ) %>%
  mutate(
    Time_of_Day = recode(
      Time_of_Day,
      "LST_Day"   = "Day",
      "LST_Night" = "Night"
    ),
    Month = factor(
      Month,
      levels = c("Feb", "Mar", "Apr", "May")
    ),
    Year = factor(Year)
  )

# Boxplot: City × Year facets
bplot <- ggplot(
  lst_long,
  aes(x = Month, y = Temperature, fill = Time_of_Day)
) +
  geom_boxplot(
    outlier.size = 0.5,
    outlier.colour = "#CD0000",
    width = 0.75,
    position = position_dodge(0.75)
  ) +
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 23,
    size = 1.5,
    color = "#545454",
    position = position_dodge(0.75)
  ) +
  facet_grid(City ~ Year) +
  scale_fill_manual(
    values = c(
      "Day"   = "#CD2626",
      "Night" = "#EE9A00"
    ),
    name = "Time of Day"
  ) +
  labs(
    x = "Month",
    y = "LST (°C)"
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(face = "bold", size = 10),
    axis.text.x  = element_text(size = 10, face = "bold"),
    axis.text.y  = element_text(size = 10, face = "bold"),
    axis.title.x = element_text(size = 10, face = "bold"),
    axis.title.y = element_text(size = 10, face = "bold"),
    legend.position = "bottom"
  )

bplot

# ============================================================
# PART 4: LST Getis–Ord Gi* Hotspot Analysis for Multiple Cities
# ============================================================

# Global parameters
base_path  <- "/Users/drdesalewmoges/Documents/Papers/RS_Heat_Exposure/Data/LST/"
output_dir <- "/Users/drdesalewmoges/Documents/Papers/RS_Heat_Exposure/Output/Hotspot"
k_neighbors <- 6
utm_crs <- "+proj=utm +zone=37 +datum=WGS84 +units=m +no_defs"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Gi* Hotspot processing function
process_LST_hotspot <- function(lst_raster, city, k, utm_crs) {
  
  # Raster → data frame
  lst_df <- as.data.frame(lst_raster, xy = TRUE, na.rm = TRUE)
  colnames(lst_df)[3] <- "LST"
  
  # Data frame → SpatialPoints
  coordinates(lst_df) <- ~ x + y
  proj4string(lst_df) <- CRS("+proj=longlat +datum=WGS84")
  
  # Reproject to UTM
  lst_proj <- spTransform(lst_df, CRS(utm_crs))
  coords <- coordinates(lst_proj)
  
  # Spatial weights (k-NN)
  nb <- knn2nb(knearneigh(coords, k = k))
  lw <- nb2listw(nb, style = "W")
  
  # Gi* statistic
  lst_proj$GiZScore <- as.numeric(localG(lst_proj$LST, lw))
  
  # Convert to sf
  lst_sf <- st_as_sf(lst_proj)
  
  # Hotspot classification
  lst_sf <- lst_sf %>%
    mutate(
      HotspotClass = case_when(
        GiZScore <= -2.58 ~ "Strong Coldspot (99%)",
        GiZScore <= -1.96 ~ "Moderate Coldspot (95%)",
        GiZScore <= -1.65 ~ "Weak Coldspot (90%)",
        GiZScore <=  1.65 ~ "Not Significant",
        GiZScore <=  1.96 ~ "Weak Hotspot (90%)",
        GiZScore <=  2.58 ~ "Moderate Hotspot (95%)",
        TRUE             ~ "Strong Hotspot (99%)"
      ),
      HotspotClass = factor(
        HotspotClass,
        levels = c(
          "Strong Coldspot (99%)",
          "Moderate Coldspot (95%)",
          "Weak Coldspot (90%)",
          "Not Significant",
          "Weak Hotspot (90%)",
          "Moderate Hotspot (95%)",
          "Strong Hotspot (99%)"
        )
      )
    )
  
  # Color palette
  hotspot_cols <- c(
    "Strong Coldspot (99%)"   = "#053061",
    "Moderate Coldspot (95%)" = "#2166AC",
    "Weak Coldspot (90%)"     = "#4393C3",
    "Not Significant"         = "#EEEE00",
    "Weak Hotspot (90%)"      = "#FDAE61",
    "Moderate Hotspot (95%)"  = "#D6604D",
    "Strong Hotspot (99%)"    = "#B2182B"
  )
  
  # --- Hotspot map
  map_plot <- ggplot(lst_sf) +
    geom_sf(aes(color = HotspotClass), size = 2.2, alpha = 0.85) +
    scale_color_manual(values = hotspot_cols) +
    labs(
      title = paste(city, "LST Hotspot Analysis (Getis–Ord Gi*)")
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.position = "right",
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      panel.border = element_blank()
    )
  
  # --- Gi* Z-score histogram
  hist_plot <- ggplot(lst_sf, aes(GiZScore)) +
    geom_histogram(binwidth = 0.5, fill = "skyblue", color = "white") +
    geom_vline(xintercept = c(-1.96, 1.96), linetype = "dashed", color = "red") +
    labs(
      title = paste(city, "Gi* Z-Score Distribution"),
      x = "Z-Score", y = "Frequency"
    ) +
    theme_bw() +
    theme(panel.grid = element_blank())
  
  # --- Area percentage
  area_df <- lst_sf %>%
    st_drop_geometry() %>%
    count(HotspotClass) %>%
    mutate(percent = 100 * n / sum(n))
  
  bar_plot <- ggplot(area_df, aes(HotspotClass, percent, fill = HotspotClass)) +
    geom_col() +
    scale_fill_manual(values = hotspot_cols) +
    labs(
      title = paste(city, "Hotspot Class Area Percentage"),
      x = "Hotspot Class", y = "Area (%)"
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none",
      plot.title = element_text(face = "bold", hjust = 0.5)
    )
  
  list(
    hotspot_map = map_plot,
    histogram   = hist_plot,
    bar_chart   = bar_plot,
    z_scores    = lst_sf$GiZScore
  )
}

# Load rasters (all cities)
cities <- c("Addis", "Adama", "Jimma", "Harar")

lst_files <- setNames(
  lapply(cities, function(cty) {
    rast(paste0(base_path, cty, "_LST_2024.tif"))
  }),
  cities
)

# Run analysis for all cities
city_results <- lapply(
  names(lst_files),
  function(city) process_LST_hotspot(
    lst_raster = lst_files[[city]],
    city       = city,
    k          = k_neighbors,
    utm_crs    = utm_crs
  )
)
names(city_results) <- cities

# Export outputs
for (city in cities) {
  
  ggsave(paste0(output_dir, "/", city, "_HotspotMap.png"),
         city_results[[city]]$hotspot_map,
         width = 10, height = 8, dpi = 300)
  
  ggsave(paste0(output_dir, "/", city, "_GiZScoreHistogram.png"),
         city_results[[city]]$histogram,
         width = 8, height = 6, dpi = 300)
  
  ggsave(paste0(output_dir, "/", city, "_HotspotBarChart.png"),
         city_results[[city]]$bar_chart,
         width = 8, height = 6, dpi = 300)
  
  message("Saved outputs for ", city)
}

# Gi* Z-score skewness summary
city_skew_df <- data.frame(
  City = cities,
  GiZScore_Skewness = sapply(cities, function(cty) {
    skewness(city_results[[cty]]$z_scores, na.rm = TRUE)
  })
)

city_skew_df

# ============================================================
# PART 5: LST-Based Heatwave Detection (Daytime & Nighttime)
# Multiple Cities | CTX90pct | 2021–2024
# ============================================================

# Generic heatwave detection function
detect_heatwaves <- function(file_path,
                             city_name,
                             time_period = c("Daytime", "Nighttime"),
                             lst_column,
                             output_plot_path) {
  
  time_period <- match.arg(time_period)
  
  # --- Load & clean data
  df <- read.csv(file_path) %>%
    mutate(
      date = as.Date(system.time_start, format = "%b %d, %Y"),
      lst  = .data[[lst_column]],
      doy  = yday(date),
      year = year(date)
    ) %>%
    select(date, lst, doy, year) %>%
    arrange(date) %>%
    drop_na()
  
  # --- CTX90pct threshold by DOY (smoothed)
  doy_thresholds <- df %>%
    group_by(doy) %>%
    summarize(
      raw_p90 = quantile(lst, 0.9, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      smoothed_p90 = rollapply(raw_p90, 15, mean,
                               align = "center",
                               fill = NA, na.rm = TRUE),
      smoothed_p90 = na.approx(smoothed_p90, rule = 2)
    )
  
  df <- df %>%
    left_join(doy_thresholds, by = "doy") %>%
    mutate(hot_day = lst > smoothed_p90)
  
  # --- Identify heatwave days (≥3 consecutive days)
  df <- df %>%
    group_by(year) %>%
    mutate(
      run_id = rleid(hot_day),
      run_length = sequence(rle(hot_day)$lengths),
      heatwave_day = hot_day & run_length >= 3
    ) %>%
    ungroup()
  
  # --- Extract heatwave events
  events <- df %>%
    filter(heatwave_day) %>%
    group_by(year, run_id) %>%
    summarize(
      start_date = min(date),
      end_date   = max(date),
      duration   = as.numeric(end_date - start_date) + 1,
      mean_temp  = mean(lst),
      max_temp   = max(lst),
      mean_excess = mean(lst - smoothed_p90),
      .groups = "drop"
    ) %>%
    filter(duration >= 3) %>%
    arrange(start_date)
  
  # --- Plot
  heatwave_plot <- ggplot(df, aes(x = date)) +
    geom_line(aes(y = lst, color = paste(time_period, "LST")),
              alpha = 0.7, linewidth = 0.7) +
    geom_line(aes(y = smoothed_p90, color = "Seasonal Threshold"),
              linetype = "dashed", linewidth = 1) +
    geom_point(
      data = filter(df, heatwave_day),
      aes(y = lst, color = "Heatwave Day"),
      size = 3
    ) +
    scale_color_manual(
      values = c(
        paste(time_period, "LST") = "#003366",
        "Seasonal Threshold" = "#A020F0",
        "Heatwave Day" = "#CD0000"
      ),
      name = ""
    ) +
    facet_wrap(~ year, scales = "free_x", ncol = 4) +
    labs(
      title = paste(
        time_period,
        "LST-Derived Heatwave Events in",
        city_name,
        "(CTX90pct, 2021–2024)"
      ),
      x = "Date",
      y = "LST (°C)"
    ) +
    theme_bw(base_size = 18) +
    theme(
      legend.position = "bottom",
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.title = element_text(face = "bold"),
      strip.text = element_text(face = "bold", size = 18)
    )
  
  ggsave(output_plot_path, heatwave_plot,
         width = 12, height = 5, dpi = 300)
  
  # --- Summary table
  summary_tbl <- events %>%
    mutate(
      City = city_name,
      Period = time_period,
      Date_Range = paste(format(start_date, "%b %d"),
                         "-", format(end_date, "%b %d")),
      Intensity_C = round(mean_excess, 1)
    ) %>%
    select(
      City, Period, Year = year,
      Duration, Date_Range, Intensity_C
    )
  
  # --- Console summary
  cat("\n", toupper(time_period), "HEATWAVE SUMMARY for", city_name, "\n")
  cat("--------------------------------------------------\n")
  cat("Total heatwave days :", sum(df$heatwave_day), "\n")
  cat("Total events        :", nrow(events), "\n")
  cat("Mean duration (days):", round(mean(events$duration), 1), "\n")
  cat("Max duration (days) :", max(events$duration), "\n")
  
  invisible(list(
    city = city_name,
    period = time_period,
    data = df,
    events = events,
    summary = summary_tbl,
    plot = heatwave_plot
  ))
}

# City configuration (Day & Night)
cities <- tribble(
  ~city, ~day_file, ~night_file,
  "Addis Ababa", "addis_day.csv", "addis_night.csv",
  "Adama",       "adama_day.csv", "adama_night.csv",
  "Harar",       "harar_day.csv", "harar_night.csv",
  "Jimma",       "jimma_day.csv", "jimma_night.csv"
)

data_dir   <- "/Users/drdesalewmoges/Documents/Papers/RS_Heat_Exposure/Data/LST/"
output_dir <- "/Users/drdesalewmoges/Documents/Papers/RS_Heat_Exposure/Output/"

# Run heatwave detection for all cities & periods
results <- pmap(
  cities,
  function(city, day_file, night_file) {
    
    list(
      Daytime = detect_heatwaves(
        file_path = paste0(data_dir, day_file),
        city_name = city,
        time_period = "Daytime",
        lst_column = "LST_Day_C",
        output_plot_path = paste0(output_dir, "daytime_heatwave_", tolower(gsub(" ", "_", city)), ".png")
      ),
      Nighttime = detect_heatwaves(
        file_path = paste0(data_dir, night_file),
        city_name = city,
        time_period = "Nighttime",
        lst_column = "LST_Night_C",
        output_plot_path = paste0(output_dir, "nighttime_heatwave_", tolower(gsub(" ", "_", city)), ".png")
      )
    )
  }
)

# Combined summary table (all cities & periods)
heatwave_summary_all <- bind_rows(
  lapply(results, function(x) {
    bind_rows(x$Daytime$summary, x$Nighttime$summary)
  })
)

heatwave_summary_all

# ============================================================
# PART 6: Mean Monthly Land Surface Temperature Patterns Across Major Ethiopian Cities (2021–2024)
# ============================================================

# Define cities and paths
cities <- c("Addis", "Adama", "Jimma", "Harar")

data_dir   <- "/Volumes/Data/Papers/Heat_Exposure/Data/LST/"
output_dir <- "/Volumes/Data/Papers/Heat_Exposure/Output/"

# Common color palette and breaks
breaks_common <- round(seq(5, 50, length.out = 7))

palette_lst <- colorRampPalette(c(
  "#0502e6", "#0602ff", "#235cb1", "#307ef3", "#269db1", "#30c8e2", "#32d3ef",
  "#3be285", "#3ff38f", "#86e26f", "#3ae237", "#b5e22e", "#d6e21f",
  "#fff705", "#ffd611", "#ffb613", "#ff8b13", "#ff6e08",
  "#ff500d", "#ff0000", "#de0101", "#c21301"
))(length(breaks_common) - 1)

# Loop through cities
for (city in cities) {
  
  # Load monthly rasters
  monthly_lst <- rast(
    paste0(data_dir, "Mean_LST_", city, "_Month_", 1:12, ".tif")
  )
  
  names(monthly_lst) <- paste0("Month_", 1:12)
  
  # Raster → data frame
  lst_df <- as.data.frame(monthly_lst, xy = TRUE) %>%
    pivot_longer(
      cols = starts_with("Month_"),
      names_to = "Month",
      values_to = "LST"
    ) %>%
    mutate(
      Month = factor(
        Month,
        levels = paste0("Month_", 1:12),
        labels = month.name
      )
    )
  
  # Plot
  p <- ggplot(lst_df) +
    geom_raster(aes(x = x, y = y, fill = LST)) +
    scale_fill_gradientn(
      colors = palette_lst,
      limits = c(5, 50),
      breaks = breaks_common,
      name = "LST (°C)"
    ) +
    facet_wrap(~ Month, ncol = 4) +
    coord_equal() +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.title = element_blank(),
      axis.text  = element_blank(),
      axis.ticks = element_blank(),
      strip.text = element_text(face = "bold", size = 14),
      legend.position = "bottom"
    ) +
    labs(
      title = paste("Mean Monthly Land Surface Temperature –", city)
    )
  
  # Save plot
  ggsave(
    filename = paste0(output_dir, "monthly_lst_", tolower(city), ".png"),
    plot = p,
    width = 10,
    height = 6,
    dpi = 300
  )
  
  cat("Saved monthly LST map for", city, "\n")
}



