#WIO VELCOTIY FOR MAINA JM
rm(list = ls())

land <- ne_download(scale = 50, type = "countries", category ="cultural", returnclass = 'sf')
land <- land %>% filter(!CONTINENT %in% c("Seven seas (open ocean)", "Antarctica"))

library(ggplot2)
library(ncdf4)
library(raster)
library(sf)
library(rgeos)
library(rasterVis)
library(gridExtra)
library(doParallel)
library(foreach)
library(scales)
library(data.table)
library(repmis)
library(rgdal)
library(xts)
library(geosphere)

mainDir <- "E:/PUBLISHED/CONS_FOREST_RISKS"
setwd(mainDir)

#lsfolder = "E:/LARGE FILE STORAGE (LFS)/CLIM-CCVA/"
luh_ref = "E:/PUBLISHED/NCC_PROTECTED_AREAS/LULC velocity/ProcessingHausz/RawDatasets/LUH2 v2h 850-2015 AD/states.nc"
luh_fut = "E:/PUBLISHED/NCC_PROTECTED_AREAS/LULC velocity/ProcessingHausz/RawDatasets/LUH2_v2 - 2015-2100/rcp8.5/ssp585_states_2015-2100.nc"

nc_varname <- c("primf","secdf","urban","c3ann","c4ann","c3per","c4per","c3nfx","pastr","range")

#Impact coefs
imp_coef <- c(1,0.855,0.645,0.460,0.460,0.635,0.635,0.585,0.533,0.550)
imp_coef <- scales::rescale(imp_coef, c(0.1, 1))

# Loop over the filenames, convert NetCDFs to raster stack
# reference period = 1971-2000
r_base <- lapply(1:length(nc_varname),function(x){
  raster::brick(luh_ref, varname = nc_varname[[x]])[[1121:1150]]*imp_coef[[x]]
  })
names(r_base) <- nc_varname
r_base <- stackApply(stack(r_base), indices = rep(seq(1:30), 10), sum, na.rm = TRUE)
ref_mean <- (1-calc(r_base, mean, na.rm = TRUE))
ref_mean <- mask(ref_mean, land)
plot(ref_mean)

#Projected period = 2021-2050
r_futu <- lapply(1:length(nc_varname),function(x){
  raster::brick(luh_fut, varname = nc_varname[[x]])[[6:35]]*imp_coef[[x]]
  })
names(r_futu) <- nc_varname
r_futu <- stackApply(stack(r_futu), indices = rep(seq(1:30), 10), sum, na.rm = TRUE)
r_futu <- mask(r_futu, land)
r_futu <- (1-r_futu)
plot(r_futu)

# Functions to calculate gradient-based velocity of climate change
# After Garcia Molinos et al.2020 (https://doi.org/10.1111/2041-210X.13295)
source(paste(mainDir, "/codes/velocity/v1_spatGrad.R", sep=""))
source(paste(mainDir, "/codes/velocity/v2_tempTrend.R", sep=""))
source(paste(mainDir, "/codes/velocity/v3_sumSeries.R", sep=""))
source(paste(mainDir, "/codes/velocity/v4_angulo.R", sep=""))

#The VoCC functionhttp://127.0.0.1:19395/graphics/plot_zoom_png?width=3491&height=1390
gVel <- function(r, ref_lyr){
  # temporal trend
  vt <- slpFUN(r)
  # spatial gradient
  vg <- spatGrad(ref_lyr, th = 0.001, projected = FALSE)
  # climate velocity
  VoCC <- vt[[1]]/vg[[1]]
    # velocity angles have opposite direction to the spatial climatic gradient if warming and same direction (cold to warm) if cooling
    ind <- which(getValues(VoCC) > 0)
    VoCCang <- vg[[2]]
    VoCCang[ind] <- vg[[2]][ind] + 180
    VoCCang[] <- ifelse(VoCCang[] >= 360, VoCCang[] - 360, VoCCang[])
    output <- stack(VoCC,VoCCang)
    names(output) <- c("voccMag", "voccAng")
    return(output)
}

#Run Velocity for each Scenario and both periods
lu_vel <- gVel(r_futu, ref_mean)
plot(log(abs(lu_vel[[1]])+1))
hist(lu_vel[[1]], breaks = 50)

#write to drive as NC file
writeRaster(lu_vel[[1]],
            filename=file.path("vocc_luh2_ssp585.tif"), overwrite=TRUE)





