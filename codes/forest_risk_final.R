#The risk degraded and intact forests face to climate and land-use change
#Author: ERNEST FRIMPONG ASAMOAH, Ph.D
#Research fellow

#install.packages
if(!requireNamespace("terra", quietly=TRUE))  install.packages("terra", quiet=TRUE, dependencies=TRUE)
if(!requireNamespace("ggplot2", quietly=TRUE)) install.packages("ggplot2", quiet=TRUE, dependencies=TRUE)
if(!requireNamespace("sf", quietly=TRUE)) install.packages("sf", quiet=TRUE, dependencies=TRUE)
if(!requireNamespace("dplyr", quietly=TRUE)) install.packages("dplyr", quiet=TRUE, dependencies=TRUE)
if(!requireNamespace("spatialEco", quietly=TRUE)) install.packages("spatialEco", quiet=TRUE, dependencies=TRUE)
if(!requireNamespace("rnaturalearth", quietly=TRUE)) install.packages("rnaturalearth", quiet=TRUE, dependencies=TRUE)

rm(list=ls())

# Load required package
library(rnaturalearth);library(sf);library(dplyr);library(ggplot2);library(terra)#;library(tidyverse)

#Import file
all.data <- readRDS("data/riskData.rds")

# Rename some new variables
colnames(all.data)[colnames(all.data) == "eco_id"] <- "ECO_ID"
colnames(all.data)[colnames(all.data) == "NatObjIDs"] <- "CID"

# Import admins attributes
gadmAttr <- as.data.frame(readxl::read_excel("data/gadm_attributes.xlsx"))[1:3]
colnames(gadmAttr)[colnames(gadmAttr) == "OBJECTID"] <- "CID"

data_info <- as.data.frame(readr::read_csv("data/ecoregion_attributes.csv")) #Merge all ecoregional data

# Merge but keep all risk dataset
all.data <- merge(all.data, gadmAttr, by = "CID", all.x = TRUE)
all.data <- merge(all.data, data_info, by = "ECO_ID", all.x = TRUE)

###########################################################################
# CLIMATE EXPOSURE
###########################################################################

# import geometric mean function
# https://stackoverflow.com/questions/2602583/geometric-mean-is-there-a-built-in
source("./codes/gmean.R")
# source("./codes/var_transformation.R")

# Store some few columns to be used for plotting later
# df_st <- all.data[,c("ECO_ID","CID","x","y","GID_0","NAME_0")]
# write.csv(df_st, "data/risk_coord.csv", row.names = F)

# INDIVIDUAL THRESHOLDS
source("./codes/exposure_func.R")
variable = c("temp","prec")

all.data$baseline = apply(all.data[colnames(all.data)[grep("current_2000",colnames(all.data))]],1,mean)
(exposureNormed <- threshold.based.exposure(data = all.data, scenario = c("hadgem","mpiesm","noresm"), id = "ECO_ID", th = .90, baseline = "baseline"))

# Calculate the total forest extent within each forest type
df1 <- all.data %>% group_by(ECO_ID) %>% summarise(N = n()); colnames(df1) <- c("ID", "N")
exposureNormed <- merge(exposureNormed, df1, by = "ID")

# Multimodal averaging and find exposure as percentage of the extent projected to be stressed

exposureNormed$vocc_temp_rcp26 = apply(exposureNormed[colnames(exposureNormed)[grep("temp_rcp26_",colnames(exposureNormed))]],1,gm_mean)/exposureNormed[,"N"]
exposureNormed$vocc_temp_rcp85 = apply(exposureNormed[colnames(exposureNormed)[grep("temp_rcp85_",colnames(exposureNormed))]],1,gm_mean)/exposureNormed[,"N"]
exposureNormed$vocc_prec_rcp26 = apply(exposureNormed[colnames(exposureNormed)[grep("prec_rcp26_",colnames(exposureNormed))]],1,gm_mean)/exposureNormed[,"N"]
exposureNormed$vocc_prec_rcp85 = apply(exposureNormed[colnames(exposureNormed)[grep("prec_rcp85_",colnames(exposureNormed))]],1,gm_mean)/exposureNormed[,"N"]

exposureNormed$vocc_rcp85 = apply(exposureNormed[c("vocc_prec_rcp85","vocc_temp_rcp85")],1,max)
with(exposureNormed, plot(vocc_temp_rcp85, vocc_prec_rcp85))

with(exposureNormed, hist(sqrt(vocc_temp_rcp85), breaks= 30))
with(exposureNormed, hist(sqrt(vocc_prec_rcp85), breaks= 30))
with(exposureNormed, hist((vocc_rcp85), breaks= 30))

summary(exposureNormed$vocc_rcp85*100)
nrow(subset(exposureNormed, vocc_prec_rcp85>0.66))/nrow(exposureNormed)

##############################################################################################################################
# FOREST INTEGRITY GRADIENT
###################################################################################

all.data$int_binned <- ifelse(all.data$FLII >= 9.6,1,NA)
areal_sum_high <- aggregate(int_binned ~ ECO_ID, data = all.data, sum); colnames(areal_sum_high) <- c("ID", "a_sum_high_int")
avg_integrity <- aggregate(FLII ~ ECO_ID, data = all.data, gm_mean); colnames(avg_integrity) <- c("ID", "avg_integrity")

intactness <- merge(avg_integrity, areal_sum_high, all = TRUE)
intactness$a_sum_high_int[is.na(intactness$a_sum_high_int)] <- 0

#################################################################################################################
#LAND DEGRADATION THREAT (LDR)
#################################################################################################################

#all.data$ssp5_agric <- apply( all.data[,c("ssp5_c3ann","ssp5_c3nfx","ssp5_c3per","ssp5_c4ann")],1,mean)
#https://www.nature.com/articles/s41586-024-07264-9


all.data <- all.data %>% mutate( 
  ssp5_c3ann_nd = (ifelse(ssp5_c3ann <=0, 0, ssp5_c3ann)),
  ssp5_c3nfx_nd = (ifelse(ssp5_c3nfx <=0, 0, ssp5_c3nfx)),
  ssp5_c3per_nd = (ifelse(ssp5_c3per <=0, 0, ssp5_c3per)),
  ssp5_c4ann_nd = (ifelse(ssp5_c4ann <=0, 0, ssp5_c4ann)),
  ssp5_c4per_nd = (ifelse(ssp5_c4per <=0, 0, ssp5_c4per)),
  ssp5_pastr_nd = (ifelse(ssp5_pastr <=0, 0, ssp5_pastr)),
  ssp5_urban_nd = (ifelse(ssp5_urban <=0, 0, ssp5_urban)))

all.data$ssp5_luse = apply(all.data[,c("ssp5_c3ann_nd","ssp5_c3nfx_nd","ssp5_c3per_nd","ssp5_c4ann_nd","ssp5_c4per_nd","ssp5_pastr_nd","ssp5_urban_nd")],1,function(x)(sum(x, na.rm = T)))
all.data$luse_bin <- all.data[,"ssp5_luse"]#*all.data[,"int_binned"]
summary(all.data$luse_bin)

# Aggregate to for each forest major group
lu_risk <- aggregate(luse_bin ~ ECO_ID, data = all.data, mean)
colnames(lu_risk)[colnames(lu_risk)=="ECO_ID"] <- "ID"
summary(lu_risk)

#Merge all ecoregional data together
all.data.fnal <- merge(exposureNormed, intactness, by = "ID", all = T)
all.data.fnal$intactness <- ifelse(all.data.fnal[,"a_sum_high_int"]>0, all.data.fnal[,"a_sum_high_int"]/all.data.fnal[,"N"], 0)

all.data.fnal <- merge(all.data.fnal, lu_risk, by = "ID", all = T)
# all.data.fnal$lu_risk <- ifelse(all.data.fnal[,"luse_bin"]>0, all.data.fnal[,"luse_bin"]/all.data.fnal[,"a_sum_high_int"], 0)
hist(all.data.fnal$luse_bin)

#Save to disk
write.csv(all.data.fnal, "data/risk_data_agg.csv", row.names = F)

#House tidying
rm(all.data.fnal, areal_sum_high, data_info, df1, lu_risk, intactness, exposureNormed, gadmAttr, avg_integrity)

######################################################################################################
#CLIMATE RISK INDEX FOR EARTH'S FOREST
######################################################################################################
z_score_brk <- function(x) {
  x = sign(x)*log10(abs(x)+0.0001)
  x = scale(x)[,1]
  xbin <- as.numeric(cut(x, breaks = c(-Inf, -1.282, -0.675, 0.000, 0.675, 1.282, Inf), labels = seq(0,5,1)))
  return(xbin)
}


source("./codes/gmean.R")

# Load required package
library(sf);library(RColorBrewer);library(rnaturalearth);library(dplyr);library(ggplot2);library(terra)#;library(tidyverse)

# Import coordinates for mapping
risk_cood <- read.csv("data/risk_coord.csv")

# DOWNLOAD SOME MAP ELEMENTS
bbox <- rnaturalearth::ne_download(scale=50, type="wgs84_bounding_box",category="physical", returnclass = "sf")%>%st_transform(crs = "+proj=moll")
land <- rnaturalearth::ne_download(scale=10, type="countries",category="cultural", returnclass = "sf") %>% st_transform(crs = "+proj=moll")
land <- filter(land, scalerank =="0") %>% group_by(CONTINENT) %>%summarise(N = n())
plot(land["CONTINENT"])

# Import ecoregional aggregates information
testing <- as.data.frame(readr::read_csv("data/risk_data_agg.csv"))
colnames(testing)[colnames(testing) == "ID"] = "ECO_ID"

#House Cleaning
data_info <- as.data.frame(readr::read_csv("data/ecoregion_attributes.csv")) #Merge all ecoregional data
testing <- merge(testing, data_info, by = "ECO_ID", all.x = TRUE)

testing <- filter(testing, !BIOME_NAME %in% c("Deserts & Xeric Shrublands","Flooded Grasslands & Savannas"))
testing <- filter(testing, !is.na(N)); testing <- filter(testing, N>10)
nrow(testing)
# [1] 605

# Degree of association between RCPS
with(testing, cor.test(vocc_temp_rcp85, vocc_temp_rcp26, method = "spearman"))
with(testing, plot(vocc_temp_rcp85, vocc_temp_rcp26))

# Compute t-test
m_data <- reshape::melt(testing[c("ECO_ID","vocc_temp_rcp85","vocc_prec_rcp85")], id.var = "ECO_ID")
(res <- t.test(value ~ variable, data = m_data, var.equal = F))

m_data <- reshape::melt(testing[c("ECO_ID","vocc_temp_rcp26","vocc_prec_rcp26")], id.var = "ECO_ID")
(res <- t.test(value ~ variable, data = m_data, var.equal = F))

#############################################################################################################
#TOPSIS APPROACH TO VULNERABILITY (HIGH EMISSION CLIMATE CHANGE SCENARIO)
############################################################################################################

# Function to perform TOPSIS
source("codes/topsis_func.R")
#source("./codes/var_transformation.R")

# Decision matrix (alternatives in rows, criteria in columns)
X_proj <- as.matrix(
  cbind( (testing[,"vocc_rcp85"]),
         (testing[,"intactness"]))
  )

# Weights for each criterion
W <- c(1, 1) #the final model was given equal weights

# Benefit/cost indicator for each criterion
is_benefit <- c(FALSE, TRUE)

testing$risk.rcp85 <- topsis(X_proj, W, is_benefit, std = F)
quantile(testing$risk.rcp85,probs = .90)

hist(testing$risk.rcp85, breaks = 30)
with(subset(testing, luse_bin>0.0001), plot(log(luse_bin), risk.rcp85))
with(testing, cor.test(log(luse_bin+0.0001), risk.rcp85, method = "pearson"))

filter(testing, risk.rcp85>=0.66) %>% group_by(BIOME_NAME) %>% summarise(N = n())
filter(testing, risk.rcp85>=0.66) %>% summarise(N_sum = sum(N), N = n())
filter(testing, risk.rcp85>=0.66) %>% summarise(N = median(vocc_rcp85))

summary(filter(testing, risk.rcp85<=.33) )

nrow(filter(testing, risk.rcp85<=0.33))/nrow(testing)

# impacts of mitigation
# testing$mitigation_impacts <- 100*(testing$climate_rcp85-testing$climate_rcp26)
# with(testing, cor.test(mitigation_impacts, risk.clim, method = "pearson"))
# with(testing, plot(mitigation_impacts, risk.clim))
# 
# testing$alp <- ifelse(testing$mitigation_impacts>0,1,.5)
# (mitigate_imp <- ggplot(data = testing) + 
#     geom_point(aes(risk.clim.wd, mitigation_impacts, alpha = alp), colour = "orange4")+
#     geom_hline(yintercept = 0, lty = 2)+
#     #scale_colour_viridis_c(name = "", option = "A")+
#     scale_x_continuous(name = "Vulnerability of forest with no mitigation", expand = c(0,0), limits = c(0,1), breaks = seq(0,1,0.2))+
#     scale_y_continuous(name = "Effects of strong, relative to weak mitigation \non climate exposure", expand = c(0,0), limits = c(-30, 90))+
#     theme_bw(base_size = 15)+ theme(legend.position = "none", panel.grid = element_blank()))
# 
# ggsave(plot = mitigate_luse, "figs/risk_mitigate.png", dpi = 1200, width = 5, height = 5)
# (zb <- rgeoda::natural_breaks(k=4, as.data.frame(testing[,"luse_bin"])))
# (za <- rgeoda::natural_breaks(k=4, as.data.frame(testing[,"risk.rcp85"])))

z_score_brk(testing$luse_bin)
(scale(testing$luse_bin))
testing$luse_bin_sc <- scales::rescale(log10(round(testing$luse_bin, 4)+0.001))

testing <- testing %>% 
  mutate( 
    luse_ = scale(log10(round(luse_bin, 4)+0.001), center = T, scale=T)[,1],
    luse_brk = factor(ifelse(luse_ >= 1.2816, "Very high", ifelse(luse_ >= 0.4125, "High", ifelse(luse_ >= -0.4399, "Moderate", "Low")))),
    vul_brk = factor(ifelse(risk.rcp85 >= .90, "Very high", ifelse(risk.rcp85 >= .66, "High", ifelse(risk.rcp85 >= .33, "Moderate", "Low"))))  )

(df <- testing %>% group_by(vul_brk, luse_brk) %>% summarise(N = n()) %>% ungroup())
# write.csv(df, "data.csv")


df1 <- testing[,c("ECO_ID","vul_brk","luse_brk")]
(df1 <- reshape::melt(as.data.frame(df1), id.var = "ECO_ID") %>% group_by(variable,value) %>% summarise(N = n()) %>% ungroup() %>% na.omit())
colnames(df1) <- c("variable","cat","N")

df1$cat <- factor(df1$cat, levels = c("Very high","High","Moderate","Low"))
levels(df1$variable) <- c("CV", "LCR"); df1$variable <- factor(df1$variable, levels = c("CV", "LCR"))
(riskBars <- ggplot() +
  geom_col(data = df1, aes(x = variable, N, fill = cat), width = .5)+
  scale_fill_manual(name = "", values = c("Very high" = "darkred","High" = "orange4", "Moderate" = "dodgerblue4", "Low" = "grey"))+
  scale_y_continuous(name = "Number of forest types", expand = c(0,0)) + scale_x_discrete(name = "") + theme_classic(base_size = 18) + 
  theme(legend.position = "top", axis.text.x = element_text(angle = 45, hjust = 1), legend.text = element_text(size = 10))
)

#ggsave(plot = riskBars, "figs/risk_bars.png", dpi = 1200, width = 3.5, height = 5)

# Libraries
# install.packages("remotes")
# remotes::install_github("davidsjoberg/ggsankey")
library(ggsankey)

d <- data.frame(testing[,c("vul_brk", "luse_brk")]) %>% na.omit()
names(d) <- c('FV', 'LC')

TotalCount = nrow(d)
df <- d %>%  make_long("FV", "LC")
dagg <- df %>% dplyr::group_by(x,node)%>% tally() %>%  dplyr::mutate(pct = n/TotalCount)
df2 <- merge(df, dagg, by.x = c("x","node"), by.y = c("x","node"), all.x = TRUE)
df2$node <- factor(df2$node, levels = c("Low","Moderate","High","Very high"))

df2$type <- "Global"

(pl <- ggplot(df2, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = node, 
                      # label = paste0(node," n=", n, "\n",'(',  round(pct* 100,0), '%)' )
                      # label = paste0(round(pct* 100,0),"%","\n",'('," n=", n,' )' ),
                      label = paste0(round(pct* 100,0),"%","\n",'(', n,')' )
                      )) +
  geom_sankey(flow.alpha = 0.5,  color = "gray40", show.legend = TRUE, lwd = 0.1) + 
  geom_sankey_label(size = 3, color = "black", fill= NA, 
                    label.padding = unit(0, "lines"), # Remove padding
                    label.size = NA # Remove the border completely
                    )+
  theme_bw() +
  theme(legend.position = "bottom", 
        axis.title = element_blank(), 
        axis.text.y = element_blank(), 
        axis.ticks = element_blank(), 
        axis.text.x = element_text(size = 18),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        legend.spacing.x = unit(0, 'cm'),
        legend.key.width = unit(.25, 'cm'),
        legend.key.height = unit(.25, 'cm'),
        strip.background = element_rect(color=NA, fill="grey95", linewidth=1.5, linetype="solid")
        ) + 
    scale_x_discrete(expand = c(0.01,0.01))+#scale_y_continuous(expand = c(0,0)) + 
    scale_fill_manual(name = "", values = c("Very high" = "darkred","High" = "orange4", "Moderate" = "dodgerblue4", "Low" = "grey") ))

ggsave(plot = pl, "figs/risk_sankey_new.png", dpi = 1200, width = 4, height = 4)


# levels(factor(testing$BIOME_NAME))
# (df <- testing %>% dplyr::select(BIOME_NAME, vul.bin, luse.bin) %>%
#   mutate(
#   BIOME_NAME = ifelse(BIOME_NAME == "Boreal Forests/Taiga", "BF/T", BIOME_NAME),
#   BIOME_NAME = ifelse(BIOME_NAME == "Mangroves", "MGV", BIOME_NAME),
#   BIOME_NAME = ifelse(BIOME_NAME == "Mediterranean Forests, Woodlands & Scrub", "MFWS", BIOME_NAME),
#   BIOME_NAME = ifelse(BIOME_NAME == "Montane Grasslands & Shrublands", "MGS", BIOME_NAME),
#   BIOME_NAME = ifelse(BIOME_NAME == "Temperate Broadleaf & Mixed Forests", "TBMF", BIOME_NAME),
#   BIOME_NAME = ifelse(BIOME_NAME == "Temperate Conifer Forests", "TCF", BIOME_NAME),
#   BIOME_NAME = ifelse(BIOME_NAME == "Temperate Grasslands, Savannas & Shrublands", "TGSS", BIOME_NAME),
#   BIOME_NAME = ifelse(BIOME_NAME == "Tropical & Subtropical Coniferous Forests", "TSCF", BIOME_NAME),
#   BIOME_NAME = ifelse(BIOME_NAME == "Tropical & Subtropical Dry Broadleaf Forests", "TSDBF", BIOME_NAME),
#   BIOME_NAME = ifelse(BIOME_NAME == "Tropical & Subtropical Grasslands, Savannas & Shrublands", "TSGSS", BIOME_NAME),
#   BIOME_NAME = ifelse(BIOME_NAME == "Tropical & Subtropical Moist Broadleaf Forests", "TSMBF", BIOME_NAME),
#   BIOME_NAME = ifelse(BIOME_NAME == "Tundra", "TDRA", BIOME_NAME))
# )

# d <- df %>% na.omit()
# names(d) <- c("Biome", 'CCV', 'LUR')
# biome <- unique(d$Biome)
# 
# plot_list <- lapply(1:length(biome), function(x) {
#   data = d %>% dplyr::filter(Biome %in% biome[x])
#   
#   N = nrow(data)
#   df <- data %>%  make_long(CCV, LUR)
#   dagg <- df %>% dplyr::group_by(x,node) %>% tally() %>%  dplyr::mutate(pct = n/N)
#   df2 <- merge(df, dagg, by.x = c("x","node"), by.y = c("x","node"), all.x = TRUE)
#   df2$node <- factor(df2$node, levels = c("Low","Moderate","High","Very high"))
#   
#   df2$type <- paste0(biome[x])
#   
#   pl <- ggplot(df2, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = node, 
#                          label = paste0(round(pct* 100,0),"%","\n",'('," n=", n,' )' )  )) +
#       geom_sankey(flow.alpha = 0.5, node.color = 1,  color = "gray40", show.legend = TRUE, lwd = 0.1) + 
#       geom_sankey_label(size = 1, color = "black", fill= "white",width = 0.0001,
#                         inherit.aes = TRUE)+
#     theme_bw() +facet_grid(~type)+
#       theme(legend.position = "none",
#             axis.title = element_blank(), 
#             axis.text.y = element_blank(), 
#             axis.ticks = element_blank(), 
#             panel.grid = element_blank(),
#             panel.border = element_blank(),
#             strip.background = element_rect(color=NA, fill="grey95", size=1.5, linetype="solid")      )+
#       scale_x_discrete(expand = c(0,0))+#scale_y_continuous(expand = c(0,0)) + 
#       scale_fill_manual(name = "", values = c("Very high" = "darkred","High" = "orange4", "Moderate" = "dodgerblue4", "Low" = "grey") )
# })
# (biome_plt <- gridExtra::grid.arrange(grobs = plot_list))
# ggsave(plot = biome_plt, "figs/biome_sankey_4b.png", dpi = 1200, width = 8, height = 8)


#import bivariate codes
source("codes/colmatrix.R")

# Define the number of breaks
nBreaks <- 99

# Create the colour matrix
#upperleft = "#34C21B", upperright = "#FFFFFF", bottomleft = "#595757",  bottomright = "#A874B8",
#upperleft = "#0096EB", upperright = "#820050", bottomleft = "#BEBEBE", bottomright = "#FFE60F"
col.matVul <- colmat(nbreaks = nBreaks, 
                     breakstyle = "quantile",
                     xlab = "X", ylab = "Y",
                     upperleft = "dodgerblue2",
                     upperright = "darkred",
                     bottomleft = "#BEBEBE",
                     bottomright = "#FFE60F",
                     saveLeg = FALSE, plotLeg = TRUE
)
# Retrieve bivariate colour pallet data
lgdBiv <- BivLegend$data; names(lgdBiv) <- c("binY", "binX", "BivCol", "UID")

#Climate risk and FLII overlaps
bivar.data <- cbind.data.frame( ECO_ID = testing$ECO_ID,
                                binY = ntile(testing$risk.rcp85, 99),
                                binX = ntile(testing$luse_, 99)) %>% inner_join(y = lgdBiv, by = c("binY", "binX"))

baseline_shp <- merge(risk_cood, bivar.data, by = "ECO_ID", all.y = TRUE)
(risk.map.c1 <- ggplot()+
    geom_sf(data=bbox, fill=NA,lwd=.1,colour="grey50")+
    geom_sf(data=land, fill=NA, colour = "grey20", lwd=0.1)+
    geom_raster(data=baseline_shp, aes(x,y,fill=BivCol))+
    #labs(title = "(a)")+
    scale_fill_identity() + theme_void(base_size = 14)+theme(legend.position = "none"))
ggsave(plot = risk.map.c1,"figs/ssp585x-N.png", width = 5.5, height = 3)

#Generate data for legend
testing$luse_4d_sc <- scales::rescale(log10(testing$luse_4d))*10

xR <- range(testing$luse_4d_sc);yR <- range(testing$risk.rcp85)
(legCI <- expand.grid(
  y = seq(0,1, diff(yR)/1000),
  x = seq(0,10, diff(xR)/1000)
  ) %>% 
    mutate(
      bnd = "10% Quantiles",
      binY = ntile(y,10),
      binX = ntile(x,10),
      #binY = ifelse(y>.66, 3, ifelse(y >.33, 2, 1)),
      #binX = ifelse(x>.66, 3, ifelse(x >.33, 2, 1))
    ) %>% inner_join(y = lgdBiv, by = c("binY", "binX")) %>% 
    
    ggplot() +
    geom_raster(aes(x = x, y = y, fill = BivCol)) + scale_fill_identity() +
    geom_point(data=testing, aes(x=round(luse_4d_sc,2), y=round(risk.rcp85,2),shape = "Future"), stroke = 0.5, size=.25) +
    #geom_smooth(data=testing, aes(x=round(landuse,2), y=round(risk.rcp85,2)), method = "lm",formula = y~poly(x,1), se = F) +
    theme_minimal(base_size = 25)+
    scale_x_discrete(expand = c(0,0))+
    scale_y_discrete(expand = c(0,0))+
    scale_shape_manual(name="", values = c("Future" = 19))+
    facet_wrap(bnd~.)+
    theme(legend.position = "none",
          axis.text.y = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank(),
          axis.title.y = element_text(angle = 90))+
    
    labs(x = "", y = ""))
ggsave(plot = legCI, "figs/bicolxx-N.png", dpi = 1200, width = 5, height = 5)


summary(testing$climate)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.01983 0.26654 0.33730 0.35363 0.42087 0.86875 
nrow(subset(testing, climate>=.5))/nrow(testing)
View(subset(testing, climate>=.5))



qbrks1 <- classInt::classIntervals(testing$vocc_rcp85, style = "quantile")




# (qbrks1 =  quantile((testing$vocc_rcp85), probs = seq(0,1,0.1)))
# testing$exp.bin <- cut(testing$vocc_rcp85, breaks =qbrks1, include.lowest = TRUE)
# levels(testing$exp.bin) <- c("<0.02", round(qbrks1[-c(1,2,11)], 2),"<1")

qbrks2 <- classInt::classIntervals(testing$intactness, style = "quantile")
# (qbrks2 =  quantile(abs(testing$integrity), probs = seq(0,1,0.1)))
testing$ac.bin <- cut(testing$intactness, breaks = qbrks2$brks, include.lowest = TRUE)
levels(testing$ac.bin) <- c(".48", round(qbrks2[-c(1,2,11)],2),"10")

qbrks3 <- classInt::classIntervals(abs(testing$luse_bin_sc), n = 10, style = "quantile")
# (qbrks3 <- quantile(testing$landuse, probs = seq(0,1,0.1)))
testing$deg.bin <- cut(testing$luse_bin_sc, breaks = qbrks3$brks, include.lowest = T)
levels(testing$deg.bin) <- c("<.22", round(qbrks3$brks[-c(1,2)], 2))

#PLOT EXPOSURE TO FOREST LOSS
risk.exp <- ggplot() +
  geom_sf(data=bbox,fill=NA,lwd=.05,colour="grey50")+
  geom_raster(data=merge(risk_cood, testing[,c("ECO_ID","vocc_rcp85")],by ="ECO_ID" ,all.y = TRUE),
              aes(x,y,fill=vocc_rcp85))+
  geom_sf(data=land, fill=NA, colour = "grey50", lwd=0.01)+
  # scale_fill_manual(name ="", values = colorRampPalette(brewer.pal(9, "Reds"))(10))+
  scale_fill_gradient(name ="", breaks = qbrks1$brks, label = round(nbrks1$brks, 2), low = "white", high = "red2", limits = c(0, 1) )+
  theme_void(base_size = 8)+labs(title = "A")+
  theme(legend.position = "right",
        legend.spacing.x = unit(0, 'cm'),
        legend.key.width = unit(.1, 'cm'),
        legend.key.height = unit(.25, 'cm'))
rsk.inset1 <- ggplot(testing, aes(x = vocc_rcp85, after_stat(count))) + geom_density(n=10, colour = "red4", lwd = .3)+
  theme_classic(base_size = 3)+theme(panel.grid.major = element_blank())+labs(x = "Degree of exposure", y = "Number of forest types")
risk.exp <- risk.exp + annotation_custom(ggplotGrob(rsk.inset1), ymin=-.2e6, ymax = -8e6, xmin = -17e6, xmax = -8.5e6)

#PLOT ADAPTIVE CAPACITY
risk.ac <- ggplot() +
  geom_sf(data=bbox,fill=NA,lwd=.05,colour="grey50")+
  geom_raster(data=merge(risk_cood,testing[,c("ECO_ID","intactness")],by ="ECO_ID" ,all.y = TRUE), 
              aes(x,y,fill=intactness))+
  geom_sf(data=land, fill=NA, colour = "grey50", lwd=0.01)+
  # scale_fill_manual(name ="", values = colorRampPalette(brewer.pal(9, "Greens"))(10))+
  scale_fill_gradient(name ="", breaks = qbrks2$brks, label = round(qbrks2$brks, 2), low = "white", high = "green4", limits = c(0, 1) )+
  theme_void(base_size = 8)+labs(title = "B")+
  theme(legend.position = "right",
        legend.spacing.x = unit(0, 'cm'),
        legend.key.width = unit(.1, 'cm'),
        legend.key.height = unit(.25, 'cm'))
rsk.inset3 <- ggplot(testing, aes(x = intactness, after_stat(count))) + geom_density(n=10, colour = "green4", lwd = 1)+
  theme_classic(base_size = 8)+theme(panel.grid.major = element_blank())+labs(x = "Degree of forest integrity", y = "Number of forest types")
ggsave(plot = rsk.inset3, "./figs/intrigtyline.png", dpi = 1200, width = 4, height = 4)
# risk.ac <- risk.ac + annotation_custom(ggplotGrob(rsk.inset3), ymin=-.2e6, ymax = -8e6, xmin = -17e6, xmax = -8.5e6)

#PLOT EXPOSURE TO FOREST LOSS
risk.degr <- ggplot() +
  geom_sf(data = bbox, fill = NA, lwd = .05, colour = "grey50")+
  geom_tile(data = merge(risk_cood,testing[,c("ECO_ID", "deg.bin")], by = "ECO_ID", all.y = TRUE), 
            aes(x, y, fill = deg.bin))+
  geom_sf(data = land, fill = NA, colour = "grey50", lwd=0.01)+
  # scale_fill_viridis_c(name ="", option = "A",trans="sqrt", limits = c(0,1))+
  scale_fill_manual(name ="", values = colorRampPalette(brewer.pal(9, "Blues"))(10))+
  theme_void(base_size = 8)+
  # labs(title = "C")+
  theme(legend.position = "right",
        legend.spacing.x = unit(0, 'cm'),
        legend.key.width = unit(.1, 'cm'),
        legend.key.height = unit(.25, 'cm'))
ggsave(plot = risk.degr, "./figs/degMap.png", dpi = 1200, width = 6.5, height = 4)


rsk.inset3 <- ggplot(testing, aes(x = luse_bin_sc, after_stat(count))) + geom_density(n=10, colour = "dodgerblue4", lwd = 1)+
  theme_classic(base_size = 18)+theme(panel.grid.major = element_blank())+labs(x = "Degree of integrity loss", y = "Number of forest types")
# risk.degr <- risk.degr + annotation_custom(ggplotGrob(rsk.inset3), ymin=-.2e6, ymax = -8e6, xmin = -17e6, xmax = -8.5e6)
ggsave(plot = rsk.inset3, "./figs/luseline.png", dpi = 1200, width = 4, height = 4)

library(patchwork)
figs <- (risk.exp + risk.ac + risk.degr + plot_layout(ncol = 1))
ggsave(plot = figs, here::here("./figs/cmbXX-n.png"), dpi = 1200, width = 3.5, height = 5)

# y_rng <- range(ExpRisk_shp$y)
# ExpRisk_shp$latitude <- as.numeric(as.character(
#   cut(ExpRisk_shp$y, breaks = seq(y_rng[1],y_rng[2], 55500), include.lowest = TRUE,
#       labels = seq(y_rng[1],y_rng[2], 55500)[-1])
#   ))
# dat.latitude <- ExpRisk_shp %>% group_by(latitude) %>% 
#   summarise(ssp585.clim.exp = median(ssp585.clim.exp, na.rm = TRUE),
#             adaptive_capacity = median(adaptive_capacity, na.rm = TRUE),
#             ssp585.land.exp = median(ssp585.land.exp, na.rm = TRUE))
# ggplot(na.omit(dat.latitude))+
#   geom_smooth(aes(x = latitude, y =ssp585.clim.exp, colour = "Exposure"),se=F, method = "loess")+
#   geom_smooth(aes(x = latitude, y =adaptive_capacity, colour = "Adaptive Capacity"),se=F, method = "loess")+
#   geom_smooth(aes(x = latitude, y =ssp585.land.exp, colour = "Degradation"),se=F, method = "loess")+
#   coord_flip()

(qbrks4 <- quantile(testing$vul.high.clim, probs = seq(0,1,0.1)))
testing$vul.bin <- cut(testing$vul.high.clim,breaks=qbrks4, include.lowest = TRUE)
levels(testing$vul.bin) <- c(">.017", round(qbrks4[-c(1,2,11)], 2),"<1")

#PLOT EXPOSURE TO FOREST LOSS
risk.exp <- ggplot()+
  geom_sf(data=bbox,fill=NA,lwd=.1,colour="grey50")+
  geom_raster(data=merge(dat_coods,testing[,c("ECO_ID","vul.bin")],by ="ECO_ID" ,all.y = TRUE),
              aes(x,y,fill=vul.bin))+
  geom_sf(data=land, fill=NA, colour = "grey20", lwd=0.1)+
  scale_fill_manual(name ="", values = colorRampPalette(brewer.pal(9, "Reds"))(10))+
  theme_void(base_size = 8)+
  theme(legend.position = "bottom",
        legend.spacing.x = unit(0, 'cm'),
        legend.key.width = unit(1, 'cm'),
        legend.key.height = unit(.25, 'cm'))+guides(fill=guide_legend(nrow=1, title.position = "top", title.hjust = 0.5, label.position = "bottom"))
rsk.inset1 <- ggplot(testing, aes(x = vul.high.clim, after_stat(count))) + geom_density(n=10, colour = "red4", lwd = .3)+
  theme_classic(base_size = 3)+theme(panel.grid.major = element_blank())+labs(x = "Degree of exposure", y = "Number of forest types")
ov.vuln <- risk.exp + annotation_custom(ggplotGrob(rsk.inset1), ymin=-.2e6, ymax = -8e6, xmin = -17e6, xmax = -8.5e6)

ggsave(plot = ov.vuln,"figs/ov_vuln.png", width = 5.5, height = 3)
# END HERE
#########################################################################################################################################
