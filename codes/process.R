

library(dplyr)
names(ecoregion_data)
rr <- ecoregion_data |> 
  group_by(PageName, ECO_ID) |> 
  summarise(Shape_Area = sum(eco_area, na.rm = TRUE)) |>
  
  ungroup()


rr <- rr |>
  group_by(PageName, ECO_ID) |> 
  arrange(-desc(PageName)) 

hist(rr$Shape_Area)

rr <- rr |>
  group_by(PageName) |> 
  mutate(dd = sum(Shape_Area, na.rm = TRUE)) |>  ungroup() |> 
  mutate(prop = Shape_Area/dd )

rr <- rr |>
  group_by(PageName) |>
  slice_max(prop, with_ties = FALSE)

rr <- ungroup(na.omit(rr))

rr <- select(rr, c(PageName, ECO_ID))
rr <- merge(rr, dplyr::select(ecoregion_data, -c(eco_area, ECO_ID)), by = "PageName", all.x = TRUE)
write.csv(rr, "ecoregion_attr.csv", row.names = FALSE)




# primf, primn, secdf, secdn, urban, c3ann, c4ann, c3per, c4per, c3nfx, pastr, range, secmb, secma
library(raster)
sfiles <- "D:/LUH/LUH2_v2 - 2015-2100/rcp8.5/ssp585_states_2015-2100.nc"
vlst <- c("primf", "primn", "secdf", "secdn", "urban", "c3ann", "c4ann", "c3per", "c4per", "c3nfx", "pastr", "range", "secmb", "secma")
rr2015 <- c()
rr2050 <- c()
for(var in vlst) {
  rr <- raster::stack(sfiles, varname = var)
  rr2015[[var]] <- rr[[1]]
  rr2050[[var]] <- rr[[45]]
}

rr2015 <- stack(rr2015)
names(rr2015) <- paste0(vlst, "_2015")
raster::writeRaster(rr2015, paste0(names(rr2015), ".tif"), bylayer = TRUE)

rr2050 <- stack(rr2050)
names(rr2050) <- paste0(vlst, "_2050")
raster::writeRaster(rr2050, paste0(names(rr2050), ".tif"), bylayer = TRUE)


plot(rr2050 - rr2015)

