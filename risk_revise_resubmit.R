# The risk degraded and intact forests face to climate and land-use change
# Author: ERNEST FRIMPONG ASAMOAH, Ph.D
# Research fellow

library(terra)
forest_integrity <- rast("./data/raw_data/pre_processed/forIntegrity.tif")
plot(forest_integrity)

# Import Risk Maps and Resample from 24k to 5km
climateVelocity <- rast(list.files("./data/raw_data/vocc_gcms", pattern = ".tif$", full.names =TRUE))
climateVelocity <- project(climateVelocity, forest_integrity, method = "near")
plot(climateVelocity[[1]])

# Write to drive as NC file
writeRaster(climateVelocity, 
            file.path("./data/raw_data/pre_processed", paste0(filename=names(rr),".tif")),
            overwrite = TRUE
)


# primf, primn, secdf, secdn, urban, c3ann, c4ann, c3per, c4per, c3nfx, pastr, range, secmb, secma

# SSP5-RCP8.5
land_use_states <- "E:/LUH/LUH2_v2 - 2015-2100/rcp8.5/ssp585_states_2015-2100.nc"
vlst <- c("primf", "secdf")

land_use_2050 <- list()
for(var in vlst) {
  rr <- rast(land_use_states, subds = var)
  # Get the time values
  t <- time(rr)
  
  land_use_2050[[var]] <- app(rr[[which(t >= 2045 & t <= 2055)]], "mean", na.rm = TRUE)
}

land_use_2050 <- rast(land_use_2050)
land_use_2050 <- app(land_use_2050, "sum", na.rm = TRUE)
plot(land_use_2050)

land_use_2050_proj <- project(land_use_2050, forest_integrity, method = "near")
names(land_use_2050_proj) <- "forest_cover_2050_ssp585"

writeRaster(land_use_2050_proj, 
            file.path("./data/raw_data/pre_processed", paste0(names(land_use_2050_proj), ".tif")), 
            overwrite = TRUE
)

# RCP2.6
land_use_states <- "E:/LUH/LUH2_v2 - 2015-2100/rcp2.6/ssp126_states_2015-2100.nc"
vlst <- c("primf", "secdf")

land_use_rcp26 <- list()
for(var in vlst) {
  rr <- rast(land_use_states, subds = var)
  # Get the time values
  t <- time(rr)
  
  land_use_rcp26[[var]] <- app(rr[[which(t >= 2045 & t <= 2055)]], "mean", na.rm = TRUE)
}

land_use_rcp26 <- rast(land_use_rcp26)
land_use_rcp26 <- app(land_use_rcp26, "sum", na.rm = TRUE)
plot(land_use_rcp26)

land_use_rcp26_proj <- project(land_use_rcp26, forest_integrity, method = "near")
names(land_use_rcp26_proj) <- "forest_cover_2050_ssp126"

writeRaster(land_use_rcp26_proj, 
            file.path("./data/raw_data/pre_processed", paste0(names(land_use_rcp26_proj), ".tif")), 
            overwrite = TRUE
)
# Historical 
land_use_states <- "E:/LUH/LUH2 v2h 850-2015 AD/states.nc"
land_use_current <- list()

for(var in vlst) {
  rr <- rast(land_use_states, subds = var)
  
  # Get the time values
  t <- time(rr)
  
  land_use_current[[var]] <- app(rr[[which(t >= 2005 & t <= 2015)]], "mean", na.rm = TRUE)
}

land_use_current <- rast(land_use_current)
land_use_current <- app(land_use_current, "sum", na.rm = TRUE)
plot(land_use_current)

land_use_current_proj <- project(land_use_current, forest_integrity, method = "near")

names(land_use_current_proj) <- "forest_cover_2015"
writeRaster(land_use_current_proj, 
            file.path("./data/raw_data/pre_processed", paste0(names(land_use_current_proj), ".tif")), 
            overwrite = TRUE
)

# List all files
all_data_r <- rast(list.files("./data/raw_data/pre_processed", "\\.tif$", full.names = TRUE))
all_data_r <- as.data.frame(all_data_r, xy = TRUE)
all_data_r <- dplyr::filter(all_data_r, !is.na(eco_id))

forest_data <- dplyr::filter(all_data_r, !is.na(forIntegrity))
colnames(forest_data)[colnames(forest_data) == "eco_id"] <- "ECO_ID"

# Import Ecoregions Attributes
ecoregions <- read.csv("./data/ecoregion_attributes.csv")
forest_data_final <- merge(forest_data, ecoregions, by = "ECO_ID", all.x = TRUE)

saveRDS(forest_data_final, "./data/ForestRiskDec2025.rds")

# install.packages

if(!requireNamespace("terra", quietly=TRUE))  install.packages("terra", quiet=TRUE, dependencies=TRUE)
if(!requireNamespace("ggplot2", quietly=TRUE)) install.packages("ggplot2", quiet=TRUE, dependencies=TRUE)
if(!requireNamespace("sf", quietly=TRUE)) install.packages("sf", quiet=TRUE, dependencies=TRUE)
if(!requireNamespace("dplyr", quietly=TRUE)) install.packages("dplyr", quiet=TRUE, dependencies=TRUE)
if(!requireNamespace("spatialEco", quietly=TRUE)) install.packages("spatialEco", quiet=TRUE, dependencies=TRUE)
if(!requireNamespace("rnaturalearth", quietly=TRUE)) install.packages("rnaturalearth", quiet=TRUE, dependencies=TRUE)

rm(list=ls())

# Load required package
library(rnaturalearth)
library(sf)
library(dplyr)
library(ggplot2)
library(terra)
library(tidyverse)
library(RColorBrewer)

# Import file
all.data <- readRDS("./data/ForestRiskDec2025.rds")
head(all.data)

land_ne <- ne_countries(
  scale = 50, returnclass = "sf"
) |> st_transform(crs = "ESRI:54009")
land_ne <- land_ne[,c("admin", "adm0_a3")]

bbox_ne <- ne_download(
  scale = 110, 
  type = "wgs84_bounding_box", 
  category = "physical",
  returnclass = "sf"
) |> st_transform(crs = "ESRI:54009")

###########################################################################
# CLIMATE EXPOSURE
###########################################################################
head(all.data)

all.data$vocc_hist = apply(all.data[colnames(all.data)[grep("current_2000", colnames(all.data))]], 1, mean)
all.data$vocc_rcp26 = apply(all.data[colnames(all.data)[grep("rcp26_2050", colnames(all.data))]], 1, mean)
all.data$vocc_rcp85 = apply(all.data[colnames(all.data)[grep("rcp85_2050", colnames(all.data))]], 1, mean)

# ====================================
# Climate Velocity Transformation
# ====================================

# Biological reference: poleward migration rate for birds, insects, mammals
bio_ref <- 1.69  # km/year

# Calculate transformation coefficient
# Target: bio_ref velocity transforms to 0.632 (1 - 1/e)
k <- 1 / bio_ref  # = 0.5917159763313609

# Apply exponential transformation to normalize velocity distribution
# and reflect asymptotic exposure-velocity relationship
all.data$exposure_rcp85 <- 1 - exp(-k * all.data$vocc_rcp85)
all.data$exposure_rcp26 <- 1 - exp(-k * all.data$vocc_rcp26)
all.data$exposure_hist <- 1 - exp(-k * all.data$vocc_hist)


all.data  <- all.data %>% group_by(ECO_ID) %>% mutate(
  forest_exposure_rcp85 = mean(exposure_rcp85, na.rm = TRUE),
  forest_exposure_rcp26 = mean(exposure_rcp26, na.rm = TRUE),
  forest_exposure_hist = mean(exposure_hist, na.rm = TRUE)
)

# Verification: check transformation at biological reference
cat("At biological reference (", bio_ref, " km/yr):\n")
cat("Transformed value:", 1 - exp(-k * bio_ref), "\n")
cat("Expected: 0.632\n")

#################################################################################################################
# LAND DEGRADATION THREAT (LDR)
#################################################################################################################

ggplot() + 
  geom_tile(data = all.data, aes(x, y, fill = (forest_cover_2050_ssp585 - forest_cover_2015))) +
  theme_void() + coord_fixed()+ scale_fill_viridis_c()

# Calculate proportional change (bounded) for SSP5-8.5
# Assuming LC values are percentages or areas
all.data <- all.data %>% mutate(
  forest_change_ssp585 = (forest_cover_2050_ssp585 - forest_cover_2015) / forest_cover_2015,
  forest_change_ssp126 = (forest_cover_2050_ssp126 - forest_cover_2015) / forest_cover_2015
)

# Check the distribution
hist(all.data$forest_change_ssp585, breaks = 30, 
     main = "Land Use Change Distribution (SSP5-8.5)",
     xlab = "Change in Land Cover")

# Handle edge cases
all.data$forest_change_ssp585[is.infinite(all.data$forest_change_ssp585)] <- 0
all.data$forest_change_ssp585[is.na(all.data$forest_change_ssp585)] <- 0

all.data$forest_change_ssp126[is.infinite(all.data$forest_change_ssp126)] <- 0
all.data$forest_change_ssp126[is.na(all.data$forest_change_ssp126)] <- 0

# Calculate integrity with bounds checking
# If forest loss: integrity decreases and if forest gain: integrity stays same or increases slightly
all.data <- all.data %>% mutate(
  integrity_ssp585 = forIntegrity * pmax(0, pmin(1, 1 + forest_change_ssp585)),
  integrity_ssp126 = forIntegrity * pmax(0, pmin(1, 1 + forest_change_ssp126))
)

# Estimate AC
all.data  <- all.data %>% 
  group_by(ECO_ID) %>% 
  mutate(
    total_forest_pixels = n(), 
    
    current_integrity = sum(forIntegrity, na.rm = TRUE) / (10000 * total_forest_pixels),
    future_integrity_ssp585 = sum(integrity_ssp585, na.rm = TRUE) / (10000 * total_forest_pixels),
    future_integrity_ssp126 = sum(integrity_ssp126, na.rm = TRUE) / (10000 * total_forest_pixels),
    
    # RISK METRIC
    future_risk_ssp585 = sqrt(forest_exposure_rcp85^2 + (1-future_integrity_ssp585)^2),
    future_risk_ssp126 = sqrt(forest_exposure_rcp26^2 + (1-future_integrity_ssp126)^2),
    current_risk = sqrt(forest_exposure_hist^2 + (1-current_integrity)^2),
    
    # Risk without degradation
    prime_risk = sqrt(forest_exposure_rcp85^2 + (1-current_integrity)^2)
  ) %>% ungroup()

# Create a copy
with(all.data, plot(future_risk_ssp585, current_risk))
head(all.data)

# Include lowest value in first bin
# Thresholds
breaks <- c(
  -Inf,       # Minimum
  0.63,       # Biological reference with intact forests
  0.85,       # High exposure OR high degradation
  1.10,       # Both factors elevated (85% of maximum)
  Inf         # Maximum
)

labels <- c("Low", "Moderate", "High", "Very High")

risk_colors <- c(
  "Low" = "#fee08b",        # Pale yellow
  "Moderate" = "#fdae61",   # Orange
  "High" = "#f46d43",       # Red-orange
  "Very High" = "#d73027"   # Dark red
)

all.data$risk_category_ssp585 <- cut(all.data$future_risk_ssp585, breaks = breaks, labels = labels, include.lowest = TRUE)
all.data$risk_category_ssp126 <- cut(all.data$future_risk_ssp126, breaks = breaks, labels = labels, include.lowest = TRUE)
all.data$risk_category_current <- cut(all.data$current_risk, breaks = breaks, labels = labels, include.lowest = TRUE)

# CURRENT RISK MAP
p1 <- ggplot() + 
  geom_sf(data = bbox_ne, fill = NA)+
  geom_tile(data = filter(all.data, !is.na(risk_category_current)), aes(x, y, fill = (risk_category_current)))+
  scale_fill_manual(name ="", values = risk_colors) +
  # geom_sf(data = land_ne, fill = NA)+
  theme_void() + 
  theme(legend.position = "none")

ggsave(
  plot = p1,
  filename = paste0("./Figs/Dec2025/risk_category_current.png"), 
  dpi = 1200, 
  width = 17155, 
  height = 10356,
  units = "px"
)

# SOCIOECONOMIC DEVELOPMENT PATHWAY -SSP585
p2 <- ggplot() + 
  geom_sf(data = bbox_ne, fill = NA)+
  geom_tile(data = filter(all.data, !is.na(risk_category_ssp585)), aes(x, y, fill = factor(risk_category_ssp585)))+
  scale_fill_manual(name ="", values = risk_colors) + 
  # geom_sf(data = land_ne, fill = NA) +
  theme_void()+ 
  theme(legend.position = "none")

ggsave(
  plot = p2,
  filename = paste0("./Figs/Dec2025/risk_category_ssp585.png"), 
  dpi = 1200, 
  width = 17155, 
  height = 10356,
  units = "px"
)

# SOCIOECONOMIC DEVELOPMENT PATHWAY -SSP585
p3 <- ggplot() + 
  geom_sf(data = bbox_ne, fill = NA)+
  geom_tile(data = filter(all.data, !is.na(risk_category_ssp126)), aes(x, y, fill = factor(risk_category_ssp126)))+
  scale_fill_manual(name = "", values = risk_colors) + 
  # geom_sf(data = land_ne, fill = NA) +
  theme_void()+ 
  theme(legend.position = "none")

ggsave(
  plot = p3,
  filename = paste0("./Figs/Dec2025/risk_category_ssp126.png"), 
  dpi = 1200, 
  width = 17155, 
  height = 10356,
  units = "px"
)

# (rsk.inset1 <- ggplot(all.data, aes(x = forest_exposure_rcp85, after_stat(count))) + geom_density(n = 10, colour = "red4", lwd = 1) +
#     theme_classic(base_size = 18) + theme(panel.grid.major = element_blank()) + labs(x = "Degree of exposure", y = "Number of forest types"))
# ggsave(plot = rsk.inset1, "./Figs/exposureline.png", dpi = 1200, width = 4, height = 4)


# Isolate KEY variables
testing <- aggregate(
  cbind(current_risk, future_risk_ssp585, future_risk_ssp126,
        forest_exposure_hist, current_integrity,
        forest_exposure_rcp85, future_integrity_ssp585, 
        forest_exposure_rcp26, future_integrity_ssp126,
        prime_risk) ~ 
    ECO_ID + BIOME_NUM + ECO_NAME,
  data = all.data,
  FUN = mean,
  na.rm = TRUE
)

# add the total pixel counts
forest_pixel_count <- all.data %>% group_by(ECO_ID) %>% summarise(
  N = n()
)

testing_forests <- merge(testing, forest_pixel_count, by = "ECO_ID", all.x = TRUE)
testing_forests <- filter(testing_forests, !is.na(N) & N > 10)

# testing_forests <- filter(testing_forests, !BIOME_NAME %in% c("Deserts & Xeric Shrublands","Flooded Grasslands & Savannas"))
nrow(testing_forests)
# [1] 663

head(testing_forests)
testing_forests <- testing_forests %>% 
  mutate(
    risk_category_ssp585 = cut(future_risk_ssp585, breaks = breaks, labels = labels, include.lowest = TRUE),
    risk_category_ssp126 = cut(future_risk_ssp126, breaks = breaks, labels = labels, include.lowest = TRUE),
    risk_category_current = cut(current_risk, breaks = breaks, labels = labels, include.lowest = TRUE)
  )
# GLOBAL SUMMARY
table(testing_forests$risk_category_ssp585)
prop.table(table(testing_forests$risk_category_ssp585)) * 100

table(testing_forests$risk_category_ssp126)
prop.table(table(testing_forests$risk_category_ssp126)) * 100

table(testing_forests$risk_category_current)
prop.table(table(testing_forests$risk_category_current)) * 100

with(testing_forests, 
     table(risk_category_current, risk_category_ssp585)
)

# Testing Hypothesis
# Check distribution
with(testing_forests, cor.test(future_risk, current_risk, method = "spearman"))
with(testing_forests, plot(future_risk - current_risk, future_risk))

with(testing_forests, hist(future_risk))
with(testing_forests, hist(current_risk))



library(effectsize)
cohens_d(testing_forests$prime_risk, testing_forests$future_risk)

# Effect size for Mann–Whitney U
rank_biserial(testing_forests$prime_risk, testing_forests$future_risk)



# Simplified version using only biome codes
testing_forests <- testing_forests %>%
  mutate(
    forest_domain = case_when(
      BIOME_NUM %in% c(1, 2, 3) ~ "Tropical",
      BIOME_NUM %in% c(4, 5) ~ "Temperate",
      BIOME_NUM == 6 ~ "Boreal",
      BIOME_NUM == 12 ~ "Mediterranean",
      BIOME_NUM == 14 ~ "Mangroves",
      BIOME_NUM %in% c(7, 8, 9, 10, 11, 13) ~ "Other forest types",
      TRUE ~ NA_character_
    )
  )

# Future Bars
future_risk_bin <- testing_forests %>% group_by(forest_domain, risk_category_future) %>% summarise(
  N = n() ) %>% 
  ungroup() %>% 
  group_by(forest_domain) %>% 
  mutate(grp_sum = sum(N), 
         prop_bins = N/grp_sum)

future_risk_bin$forest_domain <- factor(
  future_risk_bin$forest_domain, 
  levels = c(
    "Other forest types", "Boreal", "Mangroves", "Mediterranean", "Temperate", "Tropical" 
  )
)
futureRiskBiome <- ggplot(future_risk_bin, 
                          aes(x = forest_domain, 
                              y = prop_bins, 
                              fill = risk_category_future))+
  geom_col() + coord_flip()+
  scale_fill_manual(name ="", values = risk_colors)+
  scale_y_continuous(name = "", expand = c(0,0), 
                     labels = scales::percent)+
  scale_x_discrete(name = "") +   
  theme_classic(base_size = 20) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        axis.title = element_text(face = "bold"),
        axis.text = element_text(size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "top",
        plot.margin = margin(20, 60, 10, 10))

print(futureRiskBiome)
ggsave(plot = futureRiskBiome, 
       "./Figs/Dec2025/futureRiskBiomes.png", 
       dpi = 1200, width = 8, height = 6)

# Current Risk Bars
current_risk_bin <- testing_forests %>% 
  group_by(forest_domain, risk_category_current) %>% 
  summarise(
    N = n() ) %>% ungroup() %>% 
  group_by(forest_domain) %>% 
  mutate(grp_sum = sum(N), 
         prop_bins = N/grp_sum)

current_risk_bin$forest_domain <- factor(
  current_risk_bin$forest_domain, 
  levels = c(
    "Other forest types", "Boreal", "Mangroves", "Mediterranean", "Temperate", "Tropical" 
  )
)
CurrentRiskBiome <- ggplot(current_risk_bin, 
                           aes(x = forest_domain, 
                               y = prop_bins, 
                               fill = risk_category_current))+
  geom_col() + coord_flip()+
  scale_fill_manual(name ="", values = risk_colors)+
  scale_y_continuous(name = "", expand = c(0,0), 
                     labels = scales::percent)+
  scale_x_discrete(name = "") +   
  theme_classic(base_size = 20) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        axis.title = element_text(face = "bold"),
        axis.text = element_text(size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "top",
        plot.margin = margin(20, 60, 10, 10))

print(CurrentRiskBiome)
ggsave(plot = CurrentRiskBiome, 
       "./Figs/Dec2025/CurrentRiskBiome.png", 
       dpi = 1200, width = 8, height = 6)




library(ggplot2)

# Calculate risk change
testing_forests$risk_change <- testing_forests$future_risk - testing_forests$current_risk

# Calculate medians for reference lines
median_change <- median(testing_forests$risk_change, na.rm = TRUE)
median_future <- median(testing_forests$future_risk, na.rm = TRUE)

# Get plot limits
x_range <- range(testing_forests$risk_change, na.rm = TRUE)
y_range <- range(testing_forests$future_risk, na.rm = TRUE)
y_range[1] <- 0

# Create the plot
p <- ggplot(testing_forests, aes(x = risk_change, y = future_risk)) +
  
  # Top header bars
  annotate("rect", xmin = x_range[1], xmax = median_change, 
           ymin = y_range[2], ymax = y_range[2] + (y_range[2] - y_range[1]) * 0.05, 
           fill = "#5DADE2", alpha = 1) +
  annotate("rect", xmin = median_change, xmax = x_range[2], 
           ymin = y_range[2], ymax = y_range[2] + (y_range[2] - y_range[1]) * 0.05, 
           fill = "#E67E22", alpha = 1) +
  annotate("text", x = (x_range[1] + median_change) / 2, 
           y = y_range[2] + (y_range[2] - y_range[1]) * 0.025,
           label = "slower increase", color = "white", fontface = "bold", size = 4) +
  annotate("text", x = (median_change + x_range[2]) / 2, 
           y = y_range[2] + (y_range[2] - y_range[1]) * 0.025,
           label = "faster increase", color = "white", fontface = "bold", size = 4) +
  
  # Side label bar
  annotate("rect", xmin = x_range[2], xmax = x_range[2] + (x_range[2] - x_range[1]) * 0.05,
           ymin = median_future, ymax = y_range[2], 
           fill = "#E67E22", alpha = 1) +
  annotate("rect", xmin = x_range[2], xmax = x_range[2] + (x_range[2] - x_range[1]) * 0.05,
           ymin = y_range[1], ymax = median_future, 
           fill = "#5DADE2", alpha = 1) +
  annotate("text", x = x_range[2] + (x_range[2] - x_range[1]) * 0.025, 
           y = (median_future + y_range[2]) / 2,
           label = "higher risk", angle = 270, color = "white", 
           fontface = "bold", size = 4) +
  annotate("text", x = x_range[2] + (x_range[2] - x_range[1]) * 0.025, 
           y = (y_range[1] + median_future) / 2,
           label = "lower risk", angle = 270, color = "white", 
           fontface = "bold", size = 4) +
  
  # Reference lines
  geom_vline(xintercept = median_change, 
             linetype = "dashed", color = "orange", linewidth = 1) +
  geom_hline(yintercept = median_future, 
             linetype = "dashed", color = "orange", linewidth = 1) +
  
  # Points
  geom_point(color = "gray70", alpha = 0.6, size = 2.5) +
  
  # Scales
  scale_x_continuous(breaks = c(-0.1, 0, 0.2, 0.4, 0.6, 0.8),
                     expand = expansion(mult = c(0.02, 0.08))) + 
  scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4),
                     expand = expansion(mult = c(0.02, 0.08)),
                     limits = c(0, 1.4)) + 
  
  # Labels and theme
  labs(x = "Change in risk",
       y = "Future risk",
       title = "") +
  theme_classic(base_size = 20) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        axis.title = element_text(face = "bold"),
        axis.text = element_text(size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right",
        plot.margin = margin(20, 60, 10, 10)) +
  coord_cartesian(clip = "off")

print(p)
ggsave(plot = p, "./Figs/Dec2025/forest_risk_change_analysis.png", dpi = 1200, width = 8, height = 8)


# install.packages("remotes")
# remotes::install_github("davidsjoberg/ggsankey")

library(ggsankey)

d <- data.frame(testing_forests[,c("risk_category_current", "risk_category_future")]) %>% na.omit()
names(d) <- c('Current', 'Future')

TotalCount = nrow(d)
df <- d %>%  make_long(Current, Future)
dagg <- df %>% dplyr::group_by(x,node)%>% tally() %>%  dplyr::mutate(pct = n/TotalCount)
df2 <- merge(df, dagg, by.x = c("x","node"), by.y = c("x","node"), all.x = TRUE)
df2$node <- factor(df2$node, levels = c("Low","Moderate","High","Very High"))
df2$type <- "Global"

(pl <- ggplot(df2, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = node, 
                       label = paste0(round(pct* 100,0),"%","\n",'(',n,')') )) +
    geom_sankey(flow.alpha = 0.5,  color = NA, show.legend = TRUE, lwd = 0.01) + 
    geom_sankey_text(size = 5, color = "black", width = 0.1) +
    theme_classic(base_size = 20) +
    theme(legend.position = "top", 
          plot.title = element_text(face = "bold", hjust = 0.5),
          axis.text.x = element_text(size = 20, face = "bold"),
          axis.text.y = element_blank(), 
          axis.ticks = element_blank(), 
          panel.grid = element_blank(),
          panel.border = element_blank(),
          axis.line = element_blank(),
          axis.title = element_blank(),
          strip.background = element_rect(color=NA, fill="grey95", linewidth=1.5, linetype="solid"),
          legend.key = element_blank(),
          plot.margin = margin(20, 60, 10, 10)
    )+
    scale_x_discrete(expand = c(0,0))+
    scale_fill_manual(
      name = "", 
      values = risk_colors
    )
)

ggsave(plot = pl, "./Figs/Dec2025/risk_sankey.png", dpi = 1200, width = 8, height = 8)

