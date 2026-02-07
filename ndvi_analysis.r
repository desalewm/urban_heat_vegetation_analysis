
##### NDVI Analysis 2021–2024 #####

# ============================================================
# PART 0: Libraries and Paths
# ============================================================
library(terra)
library(sf)
library(sp)
library(spdep)
library(ggplot2)
library(tmap)
library(dplyr)
library(tidyr)
library(lubridate)
library(zoo)
library(data.table)
library(scales)
library(tidyverse)
library(Kendall)  # Mann-Kendall trend

data_dir   <- "/Users/drdesalewmoges/Documents/Papers/RS_Heat_Exposure/Data/NDVI/"
output_dir <- "/Users/drdesalewmoges/Documents/Papers/RS_Heat_Exposure/Output/NDVI/"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

cities <- c("Addis", "Adama", "Harar", "Jimma")
years  <- 2021:2024

# NDVI color palette & breaks
breaks_ndvi <- seq(0, 1, by = 0.1)
palette_ndvi <- colorRampPalette(c(
  "#CD0000","#df923d","#FFFF00","#66a000",
  "#529400","#3e8601","#207401","#056201","#004c00"
))(length(breaks_ndvi) - 1)

# ============================================================
# PART 1: Multi-City × Multi-Year NDVI Spatial Distribution
# ============================================================
ndvi_df_all <- expand.grid(City = cities, Year = years) %>%
  mutate(
    file_path = paste0(data_dir, "NDVI_", City, "_FebMay_", Year, ".tif")
  ) %>%
  rowwise() %>%
  mutate(
    data = list(
      rast(file_path) %>%
        as.data.frame(xy = TRUE, na.rm = TRUE) %>%
        rename(NDVI = 3) %>%
        mutate(City = City, Year = Year)
    )
  ) %>%
  ungroup() %>%
  select(data) %>%
  bind_rows()

p_ndvi <- ggplot(ndvi_df_all, aes(x = x, y = y, fill = NDVI)) +
  geom_raster() +
  facet_grid(City ~ Year) +
  scale_fill_gradientn(
    colours = palette_ndvi,
    limits  = c(0, 1),
    breaks  = breaks_ndvi,
    name    = "NDVI"
  ) +
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
  labs(title = "NDVI Spatial Distribution (Feb–May, 2021–2024)")

print(p_ndvi)
ggsave(paste0(output_dir, "NDVI_MultiCity_2021_2024.png"), p_ndvi,
       width = 12, height = 8, dpi = 300)

# ============================================================
# PART 2: Daily NDVI Time Series
# ============================================================
# Example: Addis Ababa (repeat for other cities)
ndvi_daily <- read.csv(paste0(data_dir, "NDVI_Daily_Addis.csv")) %>%
  mutate(
    date  = as.Date(date),
    Year  = year(date),
    Month = month(date, label = TRUE, abbr = TRUE),
    NDVI  = na.locf(NDVI, na.rm = FALSE)
  ) %>%
  filter(month(date) %in% 2:5) %>%
  arrange(Year, Month, date) %>%
  mutate(month_index = match(Month, month.abb[2:5])) %>%
  group_by(Year, Month) %>%
  mutate(
    day_in_month = row_number(),
    total_days_in_month = n(),
    month_frac = (day_in_month - 1) / total_days_in_month
  ) %>%
  ungroup() %>%
  mutate(seq_month = (Year - min(Year)) * 4 + month_index - 1 + month_frac)

month_ticks <- ndvi_daily %>%
  group_by(Year, Month) %>%
  summarize(x = mean(seq_month), label = unique(Month), .groups = "drop")

ndvi_ts <- ggplot(ndvi_daily, aes(x = seq_month, y = NDVI)) +
  geom_smooth(method = "loess", span = 0.05, size = 0.75, se = FALSE, color = "#228B22") +
  scale_y_continuous(name = "NDVI") +
  scale_x_continuous(breaks = month_ticks$x, labels = month_ticks$label) +
  labs(x = "", y = "NDVI") +
  theme_bw() +
  theme(
    axis.title.y = element_text(face = "bold", size = 10),
    axis.text.x  = element_text(face = "bold", size = 10),
    axis.text.y  = element_text(face = "bold", size = 10)
  )

print(ndvi_ts)
ggsave(paste0(output_dir, "NDVI_Daily_TimeSeries_Addis.png"), ndvi_ts,
       width = 8, height = 4, dpi = 300)

# ============================================================
# PART 3: Monthly NDVI Boxplots
# ============================================================
ndvi_bp <- ndvi_daily %>%
  mutate(Month = factor(Month, levels = month.abb[2:5])) %>%
  ggplot(aes(x = Month, y = NDVI)) +
  geom_boxplot(
    fill = "#00CD66",
    outlier.size = 2,
    outlier.colour = "#008B45",
    width = 0.8
  ) +
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 23,
    size = 1.5,
    fill = "white",
    color = "#1A1A1A"
  ) +
  facet_wrap(~ Year, nrow = 1) +
  labs(x = "", y = "NDVI") +
  theme_bw() +
  theme(
    axis.text  = element_text(face = "bold", size = 10),
    strip.text = element_text(face = "bold", size = 10),
    legend.position = "none"
  )

print(ndvi_bp)
ggsave(paste0(output_dir, "NDVI_Boxplot_FebMay.png"), ndvi_bp,
       width = 8, height = 4, dpi = 300)

# ============================================================
# PART 4: NDVI Trend Analysis (Mann-Kendall)
# ============================================================
mk_ndvi <- Kendall::MannKendall(ndvi_daily$NDVI)
print(mk_ndvi)

# ============================================================
# PART 5: Mean Monthly NDVI Maps
# ============================================================
ndvi_monthly <- rast(paste0("/Volumes/Data/Papers/Heat_Exposure/Data/NDVI/MeanMonthly_NDVI_Addis_Month_", 1:12, ".tif"))
names(ndvi_monthly) <- paste0("Month_", 1:12)

ndvi_df <- as.data.frame(ndvi_monthly, xy = TRUE) %>%
  pivot_longer(
    cols = starts_with("Month_"),
    names_to = "Month",
    values_to = "NDVI"
  ) %>%
  mutate(Month = factor(Month, levels = paste0("Month_", 1:12), labels = month.name))

breaks_common <- seq(0, 1, length.out = 5)
palette_ndvi <- colorRampPalette(c(
  "#ce7e45","#df923d","#f1b555","#fcd163","#99b718","#74a901",
  "#66a000","#529400","#3e8601","#207401","#056201","#004c00"
))(100)

p_monthly <- ggplot(ndvi_df) +
  geom_raster(aes(x = x, y = y, fill = NDVI)) +
  scale_fill_gradientn(
    colors = palette_ndvi,
    limits = c(0, 1),
    breaks = breaks_common,
    name = "NDVI"
  ) +
  facet_wrap(~ Month, ncol = 4) +
  coord_equal() +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    strip.text = element_text(face = "bold", size = 14),
    legend.position = "bottom"
  )

print(p_monthly)
ggsave(paste0(output_dir, "Monthly_NDVI_Addis.png"), p_monthly,
       width = 10, height = 6, dpi = 300)
