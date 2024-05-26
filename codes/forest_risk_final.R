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

#DOWNLOAD SOME MAP ELEMENTS
bbox <- rnaturalearth::ne_download(scale=50, type="wgs84_bounding_box",category="physical", returnclass = "sf")%>%st_transform(crs = "+proj=moll")
land <- rnaturalearth::ne_download(scale=10, type="countries",category="cultural", returnclass = "sf") %>% st_transform(crs = "+proj=moll")
land <- filter(land, scalerank =="0") %>% group_by(CONTINENT) %>%summarise(N = n())
plot(land["CONTINENT"])


#Import file
all.data <- readRDS("data/riskData.rds")
#Rename some new variables
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
#source("./codes/var_transformation.R")

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
namelist <- colnames(exposureNormed)

exposureNormed$vocc_temp_rcp26 = apply(exposureNormed[namelist[grep("temp_rcp26_",namelist)]],1,gm_mean)/exposureNormed[,"N"]
exposureNormed$vocc_temp_rcp85 = apply(exposureNormed[namelist[grep("temp_rcp85_",namelist)]],1,gm_mean)/exposureNormed[,"N"]
exposureNormed$vocc_prec_rcp26 = apply(exposureNormed[namelist[grep("prec_rcp26_",namelist)]],1,gm_mean)/exposureNormed[,"N"]
exposureNormed$vocc_prec_rcp85 = apply(exposureNormed[namelist[grep("prec_rcp85_",namelist)]],1,gm_mean)/exposureNormed[,"N"]

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

# z_score_brk <- function(x) {
#   x = sign(x)*log10(abs(x)+0.0001)
#   x = scale(x)[,1]
#   xbin <- as.numeric(cut(x, breaks = c(-Inf, -1.282, -0.675, 0.000, 0.675, 1.282, Inf), labels = seq(0,5,1)))
#   return(xbin)
# }

all.data <- all.data %>% mutate( 
  ssp5_c3ann_nd = (ifelse(ssp5_c3ann <=0, 0, ssp5_c3ann)),
  ssp5_c3nfx_nd = (ifelse(ssp5_c3nfx <=0, 0, ssp5_c3nfx)),
  ssp5_c3per_nd = (ifelse(ssp5_c3per <=0, 0, ssp5_c3per)),
  ssp5_c4ann_nd = (ifelse(ssp5_c4ann <=0, 0, ssp5_c4ann)),
  ssp5_c4per_nd = (ifelse(ssp5_c4per <=0, 0, ssp5_c4per)),
  ssp5_pastr_nd = (ifelse(ssp5_pastr <=0, 0, ssp5_pastr)),
  ssp5_urban_nd = (ifelse(ssp5_urban <=0, 0, ssp5_urban)))

all.data$ssp5_luse = apply(all.data[,c("ssp5_c3ann_nd","ssp5_c3nfx_nd","ssp5_c3per_nd","ssp5_c4ann_nd","ssp5_c4per_nd","ssp5_pastr_nd","ssp5_urban_nd")],1,function(x)(sum(x, na.rm = T)))

# all.data$luse_bin <- ifelse(all.data$ssp5_luse>0.0021, 1, 0)*all.data$int_binned
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
#all.data.fnal$lu_risk <- ifelse(all.data.fnal[,"luse_bin"]>0, all.data.fnal[,"luse_bin"]/all.data.fnal[,"a_sum_high_int"], 0)
hist(all.data.fnal$luse_bin)

#Save to disk
write.csv(all.data.fnal, "data/risk_data_agg.csv", row.names = F)



######################################################################################################
#CLIMATE RISK INDEX FOR EARTH'S FOREST
######################################################################################################

source("./codes/gmean.R")

# Load required package
library(sf);library(RColorBrewer);library(rnaturalearth);library(dplyr);library(ggplot2);library(terra)#;library(tidyverse)

bbox <- rnaturalearth::ne_download(scale=50, type="wgs84_bounding_box",category="physical", returnclass = "sf")%>%st_transform(crs = "+proj=moll")
land <- rnaturalearth::ne_download(scale=10, type="coastline",category="physical", returnclass = "sf") %>% st_transform(crs = "+proj=moll")
land <- filter(land, min_zoom %in% "0", scalerank =="0")
plot(land["min_zoom"])
#eco_shp <- st_read("E:/LARGE FILE STORAGE (LFS)/CURATED SHPs/Ecoregions2017/Ecoregions2017.shp")

# Import coordinates for mapping
risk_cood <- read.csv("data/risk_coord.csv")

#Import ecoregional aggregates information
testing <- as.data.frame(readr::read_csv("data/risk_data_agg.csv"))
colnames(testing)[colnames(testing) == "ID"] = "ECO_ID"

#House Cleaning
data_info <- as.data.frame(readr::read_csv("data/ecoregion_attributes.csv")) #Merge all ecoregional data
testing <- merge(testing, data_info, by = "ECO_ID", all.x = TRUE)

testing <- filter(testing, !BIOME_NAME %in% c("Deserts & Xeric Shrublands","Flooded Grasslands & Savannas"))
testing <- filter(testing, !is.na(N)); testing <- filter(testing, N>100 &  !is.na(integrity))
nrow(testing)
#[1] 605

#Degree of association between RCPS
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

# testing$exposure.high.ssp <- exp_transform(apply(testing[,c("vocc_temp_rcp85","vocc_prec_rcp85")], 1, gm_mean))
# testing$exposure.low.ssp <- exp_transform(apply(testing[,c("vocc_temp_rcp26","vocc_prec_rcp26")], 1, gm_mean))

testing$exposure.high.ssp <- (apply(testing[,c("vocc_temp_rcp85","vocc_prec_rcp85")], 1, max))
testing$exposure.low.ssp <- (apply(testing[,c("vocc_temp_rcp26","vocc_prec_rcp26")], 1, max))

with(testing, cor.test(exposure.high.ssp, exposure.low.ssp))
with(testing, plot(exposure.high.ssp, exposure.low.ssp))

with(testing, hist(exposure.high.ssp, breaks = 30))

# Function to perform TOPSIS
source("codes/topsis_func.R")
#source("./codes/var_transformation.R")

# Decision matrix (alternatives in rows, criteria in columns)
X_proj <- as.matrix(
  cbind( (testing[,"vocc_rcp85"]),
         (testing[,"integrity"]/10))
  )

# Weights for each criterion
W <- c(1, 1) #the final model was given equal weights

# Benefit/cost indicator for each criterion
is_benefit <- c(FALSE, TRUE)

testing$risk.rcp85 <- topsis(X_proj, W, is_benefit, std = F)
quantile(testing$risk.rcp85,probs = .90)

hist(testing$risk.rcp85, breaks = 30)
with(testing, plot(log(luse_bin), risk.rcp85))
with(testing, cor.test(luse_bin, risk.rcp85, method = "spearman"))

#impacts of mitigation
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
(zb <- rgeoda::natural_breaks(k=4, as.data.frame(testing[,"std_landuse"])))
(za <- rgeoda::natural_breaks(k=4, as.data.frame(testing[,"risk.rcp85"])))

testing <- testing %>% 
  mutate( luse.bin = factor(ifelse(std_landuse >= zb[[3]], "Very high", ifelse(std_landuse >= zb[[2]], "High", ifelse(std_landuse >= zb[[1]], "Moderate", "Low")))),
          vul.bin = factor(ifelse(risk.rcp85 >= za[[3]], "Very high", ifelse(risk.rcp85 >= za[[2]], "High", ifelse(risk.rcp85 >= za[[1]], "Moderate", "Low"))))  )
#
df1 <- testing[,c("ECO_ID","vul.bin","luse.bin")]
(df1 <-reshape::melt(as.data.frame(df1), id.var = "ECO_ID") %>% group_by(variable,value) %>% summarise(N = n()) %>% ungroup() %>% na.omit())
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

d <- data.frame(testing[,c("vul.bin", "luse.bin")]) %>% na.omit()
names(d) <- c('Vulnerability', 'Landuse')

TotalCount = nrow(d)
df <- d %>%  make_long(Vulnerability, Landuse)
dagg <- df %>% dplyr::group_by(x,node)%>% tally() %>%  dplyr::mutate(pct = n/TotalCount)
df2 <- merge(df, dagg, by.x = c("x","node"), by.y = c("x","node"), all.x = TRUE)
df2$node <- factor(df2$node, levels = c("Low","Moderate","High","Very high"))

df2$type <- "Global"

(pl <- ggplot(df2, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = node, 
                      #label = paste0(node," n=", n, "\n",'(',  round(pct* 100,0), '%)' )
                      label = paste0(round(pct* 100,0),"%","\n",'('," n=", n,' )' )
                      #label = n
                      )) +
  geom_sankey(flow.alpha = 0.5,  color = "gray40", show.legend = TRUE, lwd = 0.1) + 
  geom_sankey_label(size = 1.5, color = "black", fill= "white")+
  theme_bw() +facet_grid(~type)+
  theme(legend.position = "bottom", axis.title = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(),panel.border = element_blank(),
        strip.background = element_rect(color=NA, fill="grey95", linewidth=1.5, linetype="solid")
        )+
  scale_x_discrete(expand = c(0,0))+#scale_y_continuous(expand = c(0,0)) + 
    scale_fill_manual(name = "", values = c("Very high" = "darkred","High" = "orange4", "Moderate" = "dodgerblue4", "Low" = "grey") ))

#ggsave(plot = pl, "figs/risk_sankey_4b.png", dpi = 1200, width = 4, height = 4)


levels(factor(testing$BIOME_NAME))
(df <- testing %>% dplyr::select(BIOME_NAME, vul.bin, luse.bin) %>%
  mutate(
  BIOME_NAME = ifelse(BIOME_NAME == "Boreal Forests/Taiga", "BF/T", BIOME_NAME),
  BIOME_NAME = ifelse(BIOME_NAME == "Mangroves", "MGV", BIOME_NAME),
  BIOME_NAME = ifelse(BIOME_NAME == "Mediterranean Forests, Woodlands & Scrub", "MFWS", BIOME_NAME),
  BIOME_NAME = ifelse(BIOME_NAME == "Montane Grasslands & Shrublands", "MGS", BIOME_NAME),
  BIOME_NAME = ifelse(BIOME_NAME == "Temperate Broadleaf & Mixed Forests", "TBMF", BIOME_NAME),
  BIOME_NAME = ifelse(BIOME_NAME == "Temperate Conifer Forests", "TCF", BIOME_NAME),
  BIOME_NAME = ifelse(BIOME_NAME == "Temperate Grasslands, Savannas & Shrublands", "TGSS", BIOME_NAME),
  BIOME_NAME = ifelse(BIOME_NAME == "Tropical & Subtropical Coniferous Forests", "TSCF", BIOME_NAME),
  BIOME_NAME = ifelse(BIOME_NAME == "Tropical & Subtropical Dry Broadleaf Forests", "TSDBF", BIOME_NAME),
  BIOME_NAME = ifelse(BIOME_NAME == "Tropical & Subtropical Grasslands, Savannas & Shrublands", "TSGSS", BIOME_NAME),
  BIOME_NAME = ifelse(BIOME_NAME == "Tropical & Subtropical Moist Broadleaf Forests", "TSMBF", BIOME_NAME),
  BIOME_NAME = ifelse(BIOME_NAME == "Tundra", "TDRA", BIOME_NAME))
)

d <- df %>% na.omit()
names(d) <- c("Biome", 'CCV', 'LUR')
biome <- unique(d$Biome)

plot_list <- lapply(1:length(biome), function(x) {
  data = d %>% dplyr::filter(Biome %in% biome[x])
  
  N = nrow(data)
  df <- data %>%  make_long(CCV, LUR)
  dagg <- df %>% dplyr::group_by(x,node) %>% tally() %>%  dplyr::mutate(pct = n/N)
  df2 <- merge(df, dagg, by.x = c("x","node"), by.y = c("x","node"), all.x = TRUE)
  df2$node <- factor(df2$node, levels = c("Low","Moderate","High","Very high"))
  
  df2$type <- paste0(biome[x])
  
  pl <- ggplot(df2, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = node, 
                         label = paste0(round(pct* 100,0),"%","\n",'('," n=", n,' )' )  )) +
      geom_sankey(flow.alpha = 0.5, node.color = 1,  color = "gray40", show.legend = TRUE, lwd = 0.1) + 
      geom_sankey_label(size = 1, color = "black", fill= "white",width = 0.0001,
                        inherit.aes = TRUE)+
    theme_bw() +facet_grid(~type)+
      theme(legend.position = "none",
            axis.title = element_blank(), 
            axis.text.y = element_blank(), 
            axis.ticks = element_blank(), 
            panel.grid = element_blank(),
            panel.border = element_blank(),
            strip.background = element_rect(color=NA, fill="grey95", size=1.5, linetype="solid")      )+
      scale_x_discrete(expand = c(0,0))+#scale_y_continuous(expand = c(0,0)) + 
      scale_fill_manual(name = "", values = c("Very high" = "darkred","High" = "orange4", "Moderate" = "dodgerblue4", "Low" = "grey") )
})
(biome_plt <- gridExtra::grid.arrange(grobs = plot_list))
#ggsave(plot = biome_plt, "figs/biome_sankey_4b.png", dpi = 1200, width = 8, height = 8)


# library(ggplot2);library(moonBook);library(webr);library(ggforce);library(grid)
# 
# source("Codes/others/2_donut.R")
# df1$vulB <- factor(df1$vulB, levels = c("Low", "Moderate", "High"))
# df1$vulF <- factor(df1$vulF, levels = c("Low", "Moderate", "High"))
# 
# summary(df1$vulB)/nrow(df1)
# summary(df1$vulF)/nrow(df1)
# 
# tiff("figs/Risk.tiff", width = 5000, height = 5000, res=600,pointsize=3,bg = NA)
# enrniD <- df1 %>% 
#   mutate(Risk = vulB) %>%
#   RiskDonut(
#     aes(pies = Risk, donuts = vulF), r0 = .7,
#     ratioByGroup = T, pieLabelSize = 5,
#     donutLabelSize = 5, titlesize = 13, color = "white", pieAlpha = 0.7,
#     labelpositionThreshold = 0.3)
# dev.off() 



#import bivariate codes
source("codes/colmatrix.R")

# Define the number of breaks
nBreaks <- 10

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
# bivar.data <- cbind.data.frame( ECO_ID = testing$ECO_ID, 
#                                 binY = ntile(testing$risk.rcp85, 10),
#                                 binX = ntile(testing$std_landuse, 10)) %>% inner_join(y = lgdBiv, by = c("binY", "binX"))

testing <- testing %>% 
  mutate( luse.bin = scale(std_landuse)[,1],
          vul.bin = scale(risk.rcp85)[,1],
          
          luse.bin = ifelse(luse.bin >= 1.282, 3, ifelse(luse.bin >= 0.675, 2, 0)),
          vul.bin = ifelse(vul.bin >= 1.282, 3, ifelse(vul.bin >= 0.675, 1, 0)),
          
          risk.bin = sqrt(vul.bin*luse.bin))
# summary(testing$risk.bin)/nrow(testing)

risk.grad <- merge(dat_coods,testing[,c("ECO_ID","risk.bin")],by ="ECO_ID" ,all.y = TRUE) %>%
  filter(!risk.bin %in% "0") %>% ggplot()+
  geom_sf(data=bbox,fill=NA,lwd=.05,colour="grey50")+
  geom_raster(aes(x,y,fill=risk.bin))+
  geom_sf(data=land, fill=NA, colour = "grey50", lwd=0.01)+
  scale_fill_viridis_c()+
  #scale_fill_manual(name ="", values = c("0" = "grey95", "1" = "#FFDEAD", "2" = "#CD853F", "3" = "#D2691E", "4" = "#8B4513"), breaks = c(0,1,2,3,4), label = c("Others", "high", "", "", "Critical"))+
  theme_void(base_size = 8)+labs(title = "A")+
  theme(legend.position = c(0.2,0.5),
        #legend.spacing.x = unit(0, 'cm'),
        legend.key.width = unit(.25, 'cm'),
        legend.key.height = unit(.5, 'cm'))

ggsave(plot = risk.grad,"figs/ssp_grad.png", width = 5.5, height = 3)



bivar.data <- cbind.data.frame( ECO_ID = testing$ECO_ID,
                                climate = testing$vul.bin,
                                landuse = testing$luse.bin)

baseline_shp <- merge(dat_coods, bivar.data, by = "ECO_ID", all.y = TRUE)

(risk.map.c1 <- ggplot()+
    geom_sf(data=bbox, fill=NA,lwd=.1,colour="grey50")+
    geom_sf(data=land, fill=NA, colour = "grey20", lwd=0.1)+
    geom_raster(data=baseline_shp, aes(x,y,fill=BivCol))+
    #labs(title = "(a)")+
    scale_fill_identity() + theme_void(base_size = 14)+theme(legend.position = "none"))
ggsave(plot = risk.map.c1,"figs/ssp585x-N.png", width = 5.5, height = 3)

#Generate data for legend
xR <- range(testing$std_landuse);yR <- range(testing$risk.rcp85)
(legCI <- expand.grid(
  y = seq(round(yR[1],3),round(yR[2],3), diff(yR)/1000),
  x = seq(round(xR[1],3),round(xR[2],3), diff(xR)/1000)
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
    geom_point(data=testing, aes(x=round(std_landuse,2), y=round(risk.rcp85,2),shape = "Future"), stroke = 0.5, size=1) +
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

(qbrks1 =  quantile((testing$climate), probs = seq(0,1,0.1)))
testing$exp.bin <- cut(testing$clim_rcp85_prop, breaks =qbrks1, include.lowest = TRUE)
levels(testing$exp.bin) <- c("<0.02", round(qbrks1[-c(1,2,11)], 2),"<1")

(qbrks2 =  quantile(abs(testing$integrity), probs = seq(0,1,0.1)))
testing$ac.bin <- cut(abs(testing$integrity), breaks =qbrks2, include.lowest = TRUE)
levels(testing$ac.bin) <- c(".48", round(qbrks2[-c(1,2,11)],2),"10")

(qbrks3 <- quantile(testing$landuse, probs = seq(0,1,0.1)))
testing$deg.bin <- cut(testing$landuse,breaks =qbrks3, include.lowest = TRUE)
levels(testing$deg.bin) <- c(".21", round(qbrks3[-c(1,2,11)], 2),"<1")

#PLOT EXPOSURE TO FOREST LOSS
risk.exp <- ggplot()+
  geom_sf(data=bbox,fill=NA,lwd=.05,colour="grey50")+
  geom_raster(data=merge(dat_coods,testing[,c("ECO_ID","exp.bin")],by ="ECO_ID" ,all.y = TRUE),
              aes(x,y,fill=exp.bin))+
  geom_sf(data=land, fill=NA, colour = "grey50", lwd=0.01)+
  scale_fill_manual(name ="", values = colorRampPalette(brewer.pal(9, "Reds"))(10))+
  theme_void(base_size = 8)+labs(title = "A")+
  theme(legend.position = "right",
        legend.spacing.x = unit(0, 'cm'),
        legend.key.width = unit(.1, 'cm'),
        legend.key.height = unit(.25, 'cm'))
rsk.inset1 <- ggplot(testing, aes(x = clim_rcp85_prop, after_stat(count))) + geom_density(n=10, colour = "red4", lwd = .3)+
  theme_classic(base_size = 3)+theme(panel.grid.major = element_blank())+labs(x = "Degree of exposure", y = "Number of forest types")
risk.exp <- risk.exp + annotation_custom(ggplotGrob(rsk.inset1), ymin=-.2e6, ymax = -8e6, xmin = -17e6, xmax = -8.5e6)

#PLOT ADAPTIVE CAPACITY
risk.ac <- ggplot() +
  geom_sf(data=bbox,fill=NA,lwd=.05,colour="grey50")+
  geom_raster(data=merge(dat_coods,testing[,c("ECO_ID","ac.bin")],by ="ECO_ID" ,all.y = TRUE), 
              aes(x,y,fill=ac.bin))+
  geom_sf(data=land, fill=NA, colour = "grey50", lwd=0.01)+
  scale_fill_manual(name ="", values = colorRampPalette(brewer.pal(9, "Greens"))(10))+
  theme_void(base_size = 8)+labs(title = "B")+
  theme(legend.position = "right",
        legend.spacing.x = unit(0, 'cm'),
        legend.key.width = unit(.1, 'cm'),
        legend.key.height = unit(.25, 'cm'))
rsk.inset3 <- ggplot(testing, aes(x = integrity, after_stat(count))) + geom_density(n=10, colour = "green4", lwd = .3)+
  theme_classic(base_size = 3)+theme(panel.grid.major = element_blank())+labs(x = "Degree of forest integrity", y = "Number of forest types")
risk.ac <- risk.ac + annotation_custom(ggplotGrob(rsk.inset3), ymin=-.2e6, ymax = -8e6, xmin = -17e6, xmax = -8.5e6)


#PLOT EXPOSURE TO FOREST LOSS
risk.degr <- ggplot() +
  geom_sf(data=bbox,fill=NA,lwd=.05,colour="grey50")+
  geom_tile(data=merge(dat_coods,testing[,c("ECO_ID","deg.bin")],by ="ECO_ID" ,all.y = TRUE), 
            aes(x,y,fill=deg.bin))+
  geom_sf(data=land, fill=NA, colour = "grey50", lwd=0.01)+
  #scale_fill_viridis_c(name ="", option = "A",trans="sqrt", limits = c(0,1))+
  scale_fill_manual(name ="", values = colorRampPalette(brewer.pal(9, "Blues"))(10))+
  theme_void(base_size = 8)+labs(title = "C")+
  theme(legend.position = "right",
        legend.spacing.x = unit(0, 'cm'),
        legend.key.width = unit(.1, 'cm'),
        legend.key.height = unit(.25, 'cm'))

rsk.inset3 <- ggplot(testing, aes(x = landuse, after_stat(count))) + geom_density(n=10, colour = "dodgerblue4", lwd = .3)+
  theme_classic(base_size = 3)+theme(panel.grid.major = element_blank())+labs(x = "Degree of integrity loss", y = "Number of forest types")
risk.degr <- risk.degr + annotation_custom(ggplotGrob(rsk.inset3), ymin=-.2e6, ymax = -8e6, xmin = -17e6, xmax = -8.5e6)

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
