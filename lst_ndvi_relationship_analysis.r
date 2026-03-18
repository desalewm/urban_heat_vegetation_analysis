
# =====================================================
# Title: Remote Sensing of Urban Heat Dynamics and the Cooling Effect of Urban Green Spaces in Ethiopian Cities
# LST–NDVI Relationship Analysis
# =====================================================

# PART 0. Load libraries and define paths
library(terra)
library(sf)
library(sp)
library(spgwr)
library(ggplot2)
library(viridis)
library(dplyr)
library(tidyr)
library(gridExtra)
library(patchwork)
library(units)

# Paths
lst_base_path   <- "Data/LST/"
ndvi_base_path  <- "Data/NDVI/"
residual_path   <- "Residuals/"
plot_path       <- "Residual_Maps/"
output_gwr_path <- "GWR_Results/"

dir.create(residual_path, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_path, recursive = TRUE, showWarnings = FALSE)
dir.create(output_gwr_path, recursive = TRUE, showWarnings = FALSE)

cities_orig <- c("Addis", "Adama", "Harar", "Jimma")
cities_plot <- c("Addis Ababa", "Adama", "Harar", "Jimma")
years <- 2021:2024
target_crs <- "EPSG:32637"

# Residual color palette
palette_resid <- colorRampPalette(c(
  "#053061","#2166ac","#4393c3","#92c5de","#f7f7f7",
  "#f4a582","#d6604d","#b2182b","#67001f"
))(100)

# Lists for storing plots and results
plots_NDVI <- list()
plots_R2 <- list()
results <- list()

# PART 1: Global Regression Residual Maps
for(i in seq_along(cities_orig)){
  city_orig <- cities_orig[i]
  city_plot <- cities_plot[i]
  
  lst_stack <- rast(paste0(lst_base_path, 
                           city_orig, "_LST_",
                           years, ".tif")) %>% app(mean, na.rm=TRUE) %>% project(target_crs)
  ndvi_stack <- rast(paste0(ndvi_base_path, "NDVI_", 
                            city_orig, "_FebMay_", 
                            years, ".tif")) %>% app(mean, na.rm=TRUE) %>% project(target_crs) %>% resample(lst_stack)
  
  df <- as.data.frame(cbind(values(ndvi_stack), values(lst_stack)), xy=TRUE, na.rm=TRUE)
  colnames(df) <- c("NDVI","LST","x","y")
  
  # Fit global regression
  fit <- lm(LST ~ NDVI, data=df)
  lst_pred <- df$NDVI * coef(fit)[2] + coef(fit)[1]
  residual <- df$LST - lst_pred
  resid_df <- cbind(df[,c("x","y")], residual=residual)
  
  # Symmetric limits for color scale
  lim <- max(abs(resid_df$residual), na.rm=TRUE)
  
  map <- ggplot(resid_df) +
    geom_raster(aes(x=x, y=y, fill=residual)) +
    scale_fill_gradientn(colours=palette_resid, limits=c(-lim,lim), name="Residuals (°C)") +
    coord_equal() +
    labs(title=paste(city_plot, "Global Regression Residuals")) +
    theme_void(base_size=12) +
    theme(legend.position="bottom", legend.title=element_text(face="bold", size=14))

}

# PART 2: GWR Maps – Local NDVI Effect & Local R²
for(i in seq_along(cities_orig)){
  city_orig <- cities_orig[i]
  city_plot <- cities_plot[i]
  
  lst_stack <- rast(paste0(lst_base_path, 
                           city_orig, "_LST_", years, ".tif")) %>% app(mean, na.rm=TRUE) %>% project(target_crs)
  ndvi_stack <- rast(paste0(ndvi_base_path, "NDVI_", 
                            city_orig, "_FebMay_", years, ".tif")) %>% app(mean, na.rm=TRUE) %>% project(target_crs) %>% resample(lst_stack)
  
  df <- as.data.frame(cbind(values(ndvi_stack), values(lst_stack)), xy=TRUE, na.rm=TRUE)
  colnames(df) <- c("NDVI","LST","x","y")
  
  # Spatial conversion for GWR
  df_sf <- st_as_sf(df, coords=c("x","y"), crs=target_crs)
  df_sp <- as(df_sf, "Spatial")
  
  # GWR
  bw <- gwr.sel(LST ~ NDVI, data=df_sp, adapt=TRUE, verbose=FALSE)
  gwr_model <- gwr(LST ~ NDVI, data=df_sp, adapt=bw, hatmatrix=TRUE, se.fit=TRUE)
  
  coef_df <- as.data.frame(gwr_model$SDF)
  coords <- coordinates(df_sp)
  coef_df$x <- coords[,1]; coef_df$y <- coords[,2]
  
  pts <- vect(cbind(coords, coef_df[,c("NDVI","localR2")]), geom=c("x","y"), crs=crs(lst_stack))
  r_NDVI <- rasterize(pts, lst_stack, field="NDVI")
  r2 <- rasterize(pts, lst_stack, field="localR2")
  
  ndvi_df <- as.data.frame(r_NDVI, xy=TRUE) %>% rename(NDVI_coef=NDVI)
  p_ndvi <- ggplot(ndvi_df, aes(x=x, y=y, fill=NDVI_coef)) +
    geom_tile() + scale_fill_gradient2(low="blue", mid="white", high="red", midpoint=0) +
    coord_sf(crs=target_crs) + labs(title=paste(city_plot, "Local NDVI Effect")) + theme_bw()
  
  r2_df <- as.data.frame(r2, xy=TRUE) %>% rename(R2=localR2)
  p_r2 <- ggplot(r2_df, aes(x=x, y=y, fill=R2)) +
    geom_tile() + scale_fill_viridis(option="magma") +
    coord_sf(crs=target_crs) + labs(title=paste(city_plot, "Local R²")) + theme_bw()
  
  ggsave(paste0(output_gwr_path, "/", city_plot, "_NDVI_coef.png"), p_ndvi, width=8, height=6)
  ggsave(paste0(output_gwr_path, "/", city_plot, "_R2.png"), p_r2, width=8, height=6)
  
  plots_NDVI[[city_plot]] <- p_ndvi
  plots_R2[[city_plot]] <- p_r2
  
  # Store results
  results[[city_plot]] <- list(global_corr=cor(df$NDVI, df$LST), gwr_bw=bw, n_points=nrow(df))
}

# Combined GWR plots
wrap_plots(plots_NDVI, ncol=2) + plot_annotation(title="Local NDVI Effects on LST")
wrap_plots(plots_R2, ncol=2) + plot_annotation(title="Local R² Values")

# PART 3: Scatter Plots – LST vs NDVI
df_all <- do.call(rbind, lapply(seq_along(cities_orig), function(i){
  city_orig <- cities_orig[i]; city_plot <- cities_plot[i]
  df_city <- bind_rows(lapply(years, function(y){
    lst <- rast(paste0(lst_base_path, city_orig, "_LST_", y, ".tif"))
    ndvi <- rast(paste0(ndvi_base_path, "NDVI_", city_orig, "_FebMay_", y, ".tif"))
    ndvi <- project(ndvi, lst) %>% resample(lst) %>% clamp(0,0.8)
    lst <- aggregate(lst, fact=10, fun=mean, na.rm=TRUE)
    ndvi <- aggregate(ndvi, fact=10, fun=mean, na.rm=TRUE)
    df <- as.data.frame(cbind(values(ndvi), values(lst)), na.rm=TRUE)
    colnames(df) <- c("NDVI","LST"); df$Year <- y; df$City <- city_plot
    df
  }))
  df_city
}))

# Regression labels
reg_labels <- df_all %>%
  group_by(City, Year) %>%
  group_modify(~ tibble(label=sprintf("LST=%.2f*NDVI+%.2f", 
                                      coef(lm(LST~NDVI,data=.x))[2], 
                                      coef(lm(LST~NDVI,data=.x))[1]),
                        NDVI=0.05, LST=20))

ggplot(df_all, aes(NDVI,LST)) +
  geom_point(alpha=0.7, size=0.5, color="darkblue") +
  geom_smooth(method="lm", se=TRUE, color="red") +
  geom_text(data=reg_labels, aes(x=NDVI,y=LST,label=label), inherit.aes=FALSE) +
  facet_grid(City ~ Year) +
  labs(x="NDVI", y="LST (°C)") + theme_bw(base_size=12)


#//

