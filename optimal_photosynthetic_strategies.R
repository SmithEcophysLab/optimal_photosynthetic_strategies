# global_optimal_traits.R
## script to predict global optimal leaf traits
## still a bit of work to do!

## load functions
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

## load libraries
library(raster)
library(RColorBrewer)
library(maps)
library(tidyverse)
library(factoextra)
library(sp)
library(usmap)
library(car)
library(relaimpo)

## source optimality model and related functions
source('../optimal_vcmax_r/calc_optimal_vcmax.R')
sourceDirectory('../optimal_vcmax_r/functions')

## test optimality model
calc_optimal_vcmax()

## read in global environmental data
cru_growingseason_data <- read.csv('data/cru_growingseason/cru_growingseason.csv')[,-1]
global_z_data <- read.csv('data/watch_elevation/z_globe.csv')[,-1]

## read in phenology data

## combine datasets by lat/lon
global_data_all <- left_join(cru_growingseason_data, global_z_data)

## subset out places with tmp < 10C (lowest leaf temp for Bernacchi 2003)
nrow(subset(global_data_all, tmp > 10))/nrow(global_data_all) # 67% of all data
global_data <- subset(global_data_all, tmp > 10)

## read in MODIS land cover data (to filter out non-veg sites)
modis_2001 = raster('data/landCoverMODIS/LC_hd_global_2001.tif')
modis_2002 = raster('data/landCoverMODIS/LC_hd_global_2002.tif')
modis_2003 = raster('data/landCoverMODIS/LC_hd_global_2003.tif')
modis_2004 = raster('data/landCoverMODIS/LC_hd_global_2004.tif')
modis_2005 = raster('data/landCoverMODIS/LC_hd_global_2005.tif')
modis_2006 = raster('data/landCoverMODIS/LC_hd_global_2006.tif')
modis_2007 = raster('data/landCoverMODIS/LC_hd_global_2007.tif')
modis_2008 = raster('data/landCoverMODIS/LC_hd_global_2008.tif')
modis_2009 = raster('data/landCoverMODIS/LC_hd_global_2009.tif')
modis_2010 = raster('data/landCoverMODIS/LC_hd_global_2010.tif')
modis_2011 = raster('data/landCoverMODIS/LC_hd_global_2011.tif')
modis_2012 = raster('data/landCoverMODIS/LC_hd_global_2012.tif')

modis = overlay(modis_2001, modis_2002, modis_2003, modis_2004, modis_2005, modis_2006, 
                modis_2007, modis_2008, modis_2009, modis_2010, modis_2011, modis_2012, fun = mean)
modis[modis == 16] <- 0 #barren
modis[modis > 0] <- 1 # vegetated

## create global data raster
global_data_raster <- rasterFromXYZ(cbind(global_data$lon, global_data$lat, global_data[,3:7]))
global_data_veg_raster <- global_data_raster * modis
global_data_veg <- as.data.frame(rasterToPoints(global_data_veg_raster))
colnames(global_data_veg) <- c('lon', 'lat', 'par', 'tmp', 'vpd', 'f', 'z')

## read isotope values, calculate beta distributions
### read isotope data
iso_data <- read.csv('data/iso_data/iso_data.csv') # from cheaib et al. (2025; file = data_clean_beta.csv)
head(iso_data)

### separate into c3 and c4
iso_data_c3 <- subset(iso_data, PS_pathway == 'C3')
iso_data_c4 <- subset(iso_data, PS_pathway == 'C4')
nrow(iso_data_c3)
nrow(iso_data_c4)

### calculate beta for c3
iso_data_c3$gammastar_pa <- calc_gammastar_pa(temp = iso_data_c3$tmp, z = iso_data_c3$z)
iso_data_c3$km_pa <- calc_km_pa(temp = iso_data_c3$tmp, z = iso_data_c3$z)
iso_data_c3$nstar <- calc_nstar(temp = iso_data_c3$tmp, z = iso_data_c3$z)
iso_data_c3$vpd_kpa <- calc_vpd(temp = iso_data_c3$tmp, z = iso_data_c3$z, vpdo = iso_data_c3$vpd)
iso_data_c3$ca <- iso_data_c3$CO2 * 1e-6 * calc_patm(iso_data_c3$z)
iso_data_c3$a_frac <- 4.4
iso_data_c3$b_frac <- 28
iso_data_c3$f_frac <- 12
iso_data_c3$chi <- (iso_data_c3$big_D13 - (iso_data_c3$a_frac + (iso_data_c3$f_frac * (iso_data_c3$gammastar_pa/iso_data_c3$ca))))/
  (iso_data_c3$b_frac - iso_data_c3$a_frac)
# hist(iso_data_c3$chi)
iso_data_c3$beta <- 1.6 * iso_data_c3$nstar * iso_data_c3$vpd_kpa * 1000 *
  (((iso_data_c3$chi - (iso_data_c3$gammastar_pa/iso_data_c3$ca))^2)/
                                                 (((1- iso_data_c3$chi)^2) * (iso_data_c3$km_pa + iso_data_c3$gammastar_pa)))
# hist(subset(iso_data_c3, chi < 0.95 & chi > 0.2)$beta)
# hist(log(subset(iso_data_c3, chi < 0.95 & chi > 0.2)$beta))

### calculate beta for c4
iso_data_c4$gammastar_pa <- calc_gammastar_pa(temp = iso_data_c4$tmp, z = iso_data_c4$z)
iso_data_c4$km_pa <- calc_km_pa(temp = iso_data_c4$tmp, z = iso_data_c4$z)
iso_data_c4$kp_pa <- calc_kp_temp_pa(temp = iso_data_c4$tmp, z = iso_data_c4$z)
iso_data_c4$nstar <- calc_nstar(temp = iso_data_c4$tmp, z = iso_data_c4$z)
iso_data_c4$vpd_kpa <- calc_vpd(temp = iso_data_c4$tmp, z = iso_data_c4$z, vpdo = iso_data_c4$vpd)
iso_data_c4$ca <- iso_data_c4$CO2 * 1e-6 * calc_patm(iso_data_c4$z)
iso_data_c4$a_frac <- 4.4
iso_data_c4$b_frac <- -5.7+0.2*30
iso_data_c4$f_frac <- 12
iso_data_c4$chi <- (iso_data_c4$big_D13 - (iso_data_c4$a_frac + (iso_data_c4$f_frac * (iso_data_c4$gammastar_pa/iso_data_c4$ca))))/
  (iso_data_c4$b_frac - iso_data_c4$a_frac)
# hist(iso_data_c4$chi)
# hist(subset(iso_data_c4, chi > 0)$chi)
iso_data_c4$beta <- 1.6 * iso_data_c4$nstar * iso_data_c4$vpd_kpa * 1000 *
  (((iso_data_c4$chi)^2)/
     (((1- iso_data_c4$chi)^2) * (iso_data_c4$kp_pa)))
# hist(subset(iso_data_c4, chi < 0.95 & chi > 0.1)$beta)
# hist(log(subset(iso_data_c4, chi < 0.95 & chi > 0.1)$beta))

### calculate mean and stdev beta for c3 and c4
beta_c3_mean <- mean(log(subset(iso_data_c3, chi < 0.95 & chi > 0.2)$beta))
beta_c3_sd <- sd(log(subset(iso_data_c3, chi < 0.95 & chi > 0.2)$beta))
beta_c4_mean <- mean(log(subset(iso_data_c4, chi < 0.95 & chi > 0.1)$beta))
beta_c4_sd <- sd(log(subset(iso_data_c4, chi < 0.95 & chi > 0.1)$beta))

## read in and clean neon data
### read data
neon_data <- read.csv('data/neon/neon_core_terrestrial_metadata.csv')

### calculate nearest latitude
latitude_values <- global_data_veg$lat
neon_closest_lat <- c()
for(i in 1:length(neon_data$latitude)){
  
  temp_lat <- neon_data$latitude[i]
  closest_lat_position <- which(abs(latitude_values - temp_lat) == min(abs(latitude_values - temp_lat)))[1]
  closest_lat <- latitude_values[closest_lat_position]
  neon_closest_lat <- c(neon_closest_lat, closest_lat)
  
}
neon_data$closest_latitude <- neon_closest_lat

### calculone nearest longitude
longitude_values <- global_data_veg$lon
neon_closest_lon <- c()
for(i in 1:length(neon_data$longitude)){
  
  temp_lon <- neon_data$longitude[i]
  closest_lon_position <- which(abs(longitude_values - temp_lon) == min(abs(longitude_values - temp_lon)))[1]
  closest_lon <- longitude_values[closest_lon_position]
  neon_closest_lon <- c(neon_closest_lon, closest_lon)
  
}
neon_data$closest_longitude <- neon_closest_lon

### add climate to neon dataset
neon_data_clim <- left_join(neon_data, global_data_veg, by = c('closest_latitude' = 'lat', 'closest_longitude' = 'lon'))


#############################
### Run various models ##
#############################

## run model for c3 deciduous plants
global_optimal_traits_c3_deciduous <- calc_optimal_vcmax(pathway = 'C3',
                                                         deciduous = 'yes',
                                                         tg_c = global_data_veg$tmp, 
                                                         vpdo = global_data_veg$vpd,
                                                         paro = global_data_veg$par,
                                                         z = global_data_veg$z,
                                                         f = global_data_veg$f)

## add lat/lon
global_optimal_traits_c3_deciduous$lat <- global_data_veg$lat
global_optimal_traits_c3_deciduous$lon <- global_data_veg$lon

## run model for c3 evergreen plants
global_optimal_traits_c3_evergreen <- calc_optimal_vcmax(pathway = 'C3',
                                                         deciduous = 'no',
                                                         tg_c = global_data_veg$tmp, 
                                                         vpdo = global_data_veg$vpd,
                                                         paro = global_data_veg$par,
                                                         z = global_data_veg$z,
                                                         f = global_data_veg$f)

## add lat/lon
global_optimal_traits_c3_evergreen$lat <- global_data_veg$lat
global_optimal_traits_c3_evergreen$lon <- global_data_veg$lon

## run model for c4 deciduous plants
global_optimal_traits_c4_deciduous <- calc_optimal_vcmax(pathway = 'C4',
                                                        deciduous = 'yes',
                                                        tg_c = global_data_veg$tmp, 
                                                        vpdo = global_data_veg$vpd,
                                                        paro = global_data_veg$par,
                                                        z = global_data_veg$z,
                                                        f = global_data_veg$f)

## add lat/lon
global_optimal_traits_c4_deciduous$lat <- global_data_veg$lat
global_optimal_traits_c4_deciduous$lon <- global_data_veg$lon

## run model for neon sites with varying beta values
cbind(neon_data_clim$site_id, neon_data_clim$dominant_nlcd_classes)
### representative sites (2 of each type)
# c4 sites = cper[3], konz[6]
# c3 deciduous sites = scbi[12], unde[17]
# c3 evergreen sites = puum[11], wref[19]
# c3 mixed = harv[5], tall[15]

### cper
global_optimal_traits_cper <- data.frame()
beta_cper <- exp(rnorm(10000, beta_c4_mean, beta_c4_sd))
tmp_cper <- rnorm(10000, neon_data_clim$tmp[3], neon_data_clim$tmp[3] * 0.1)
vpd_cper <- rnorm(1000, neon_data_clim$vpd[3], neon_data_clim$vpd[3] * 0.1)
par_cper <- rnorm(10000, neon_data_clim$par[3], neon_data_clim$par[3] * 0.1)
for(i in 1:length(beta_cper)){
  
  optimal_traits <- calc_optimal_vcmax(pathway = 'C4',
                                       deciduous = 'yes',
                                       tg_c = tmp_cper[i], 
                                       vpdo = vpd_cper[i],
                                       paro = par_cper[i],
                                       z = neon_data_clim$z[3],
                                       f = neon_data_clim$f[3],
                                       beta = beta_cper[i])
  
  global_optimal_traits_cper <- rbind(global_optimal_traits_cper, optimal_traits)
  
}

## add lat/lon
global_optimal_traits_cper$lat <- neon_data_clim$closest_latitude[3]
global_optimal_traits_cper$lon <- neon_data_clim$closest_longitude[3]

### konz
global_optimal_traits_konz <- data.frame()
beta_konz <- exp(rnorm(10000, beta_c4_mean, beta_c4_sd))
tmp_konz <- rnorm(10000, neon_data_clim$tmp[6], neon_data_clim$tmp[6] * 0.1)
vpd_konz <- rnorm(10000, neon_data_clim$vpd[6], neon_data_clim$vpd[6] * 0.1)
par_konz <- rnorm(10000, neon_data_clim$par[6], neon_data_clim$par[6] * 0.1)
for(i in 1:length(beta_konz)){
  
  optimal_traits <- calc_optimal_vcmax(pathway = 'C4',
                                       deciduous = 'yes',
                                       tg_c = tmp_konz[i], 
                                       vpdo = vpd_konz[i],
                                       paro = par_konz[i],
                                       z = neon_data_clim$z[6],
                                       f = neon_data_clim$f[6],
                                       beta = beta_konz[i])
  
  global_optimal_traits_konz <- rbind(global_optimal_traits_konz, optimal_traits)
  
}

## add lat/lon
global_optimal_traits_konz$lat <- neon_data_clim$closest_latitude[6]
global_optimal_traits_konz$lon <- neon_data_clim$closest_longitude[6]

### scbi
global_optimal_traits_scbi <- data.frame()
beta_scbi <- exp(rnorm(10000, beta_c3_mean, beta_c3_sd))
tmp_scbi <- rnorm(10000, neon_data_clim$tmp[12], neon_data_clim$tmp[12] * 0.1)
vpd_scbi <- rnorm(10000, neon_data_clim$vpd[12], neon_data_clim$vpd[12] * 0.1)
par_scbi <- rnorm(10000, neon_data_clim$par[12], neon_data_clim$par[12] * 0.1)
for(i in 1:length(beta_scbi)){
  
  optimal_traits <- calc_optimal_vcmax(pathway = 'C3',
                                       deciduous = 'yes',
                                       tg_c = tmp_scbi[i], 
                                       vpdo = vpd_scbi[i],
                                       paro = par_scbi[i],
                                       z = neon_data_clim$z[12],
                                       f = neon_data_clim$f[12],
                                       beta = beta_scbi[i])
  
  global_optimal_traits_scbi <- rbind(global_optimal_traits_scbi, optimal_traits)
  
}

## add lat/lon
global_optimal_traits_scbi$lat <- neon_data_clim$closest_latitude[12]
global_optimal_traits_scbi$lon <- neon_data_clim$closest_longitude[12]

### unde
global_optimal_traits_unde <- data.frame()
beta_unde <- exp(rnorm(10000, beta_c3_mean, beta_c3_sd))
tmp_unde <- rnorm(10000, neon_data_clim$tmp[17], neon_data_clim$tmp[17] * 0.1)
vpd_unde <- rnorm(10000, neon_data_clim$vpd[17], neon_data_clim$vpd[17] * 0.1)
par_unde <- rnorm(10000, neon_data_clim$par[17], neon_data_clim$par[17] * 0.1)
for(i in 1:length(beta_unde)){
  
  optimal_traits <- calc_optimal_vcmax(pathway = 'C3',
                                       deciduous = 'yes',
                                       tg_c = tmp_unde[i], 
                                       vpdo = vpd_unde[i],
                                       paro = par_unde[i],
                                       z = neon_data_clim$z[17],
                                       f = neon_data_clim$f[17],
                                       beta = beta_unde[i])
  
  global_optimal_traits_unde <- rbind(global_optimal_traits_unde, optimal_traits)
  
}

## add lat/lon
global_optimal_traits_unde$lat <- neon_data_clim$closest_latitude[17]
global_optimal_traits_unde$lon <- neon_data_clim$closest_longitude[17]

### puum
global_optimal_traits_puum <- data.frame()
beta_puum <- exp(rnorm(10000, beta_c3_mean, beta_c3_sd))
tmp_puum <- rnorm(10000, neon_data_clim$tmp[11], neon_data_clim$tmp[11] * 0.1)
vpd_puum <- rnorm(10000, neon_data_clim$vpd[11], neon_data_clim$vpd[11] * 0.1)
par_puum <- rnorm(10000, neon_data_clim$par[11], neon_data_clim$par[11] * 0.1)
for(i in 1:length(beta_puum)){
  
  optimal_traits <- calc_optimal_vcmax(pathway = 'C3',
                                       deciduous = 'no',
                                       tg_c = tmp_puum[i], 
                                       vpdo = vpd_puum[i],
                                       paro = par_puum[i],
                                       z = neon_data_clim$z[11],
                                       f = neon_data_clim$f[11],
                                       beta = beta_puum[i])
  
  global_optimal_traits_puum <- rbind(global_optimal_traits_puum, optimal_traits)
  
}

## add lat/lon
global_optimal_traits_puum$lat <- neon_data_clim$closest_latitude[11]
global_optimal_traits_puum$lon <- neon_data_clim$closest_longitude[11]

### sjer
global_optimal_traits_sjer <- data.frame()
beta_sjer <- exp(rnorm(10000, beta_c3_mean, beta_c3_sd))
tmp_sjer <- rnorm(10000, neon_data_clim$tmp[13], neon_data_clim$tmp[13] * 0.1)
vpd_sjer <- rnorm(10000, neon_data_clim$vpd[13], neon_data_clim$vpd[13] * 0.1)
par_sjer <- rnorm(10000, neon_data_clim$par[13], neon_data_clim$par[13] * 0.1)
for(i in 1:length(beta_sjer)){
  
  optimal_traits <- calc_optimal_vcmax(pathway = 'C3',
                                       deciduous = 'no',
                                       tg_c = tmp_sjer[i], 
                                       vpdo = vpd_sjer[i],
                                       paro = par_sjer[i],
                                       z = neon_data_clim$z[13],
                                       f = neon_data_clim$f[13],
                                       beta = beta_sjer[i])
  
  global_optimal_traits_sjer <- rbind(global_optimal_traits_sjer, optimal_traits)
  
}

## add lat/lon
global_optimal_traits_sjer$lat <- neon_data_clim$closest_latitude[13]
global_optimal_traits_sjer$lon <- neon_data_clim$closest_longitude[13]

### harv_deciduous
global_optimal_traits_harv_deciduous <- data.frame()
beta_harv_deciduous <- exp(rnorm(10000, beta_c3_mean, beta_c3_sd))
tmp_harv_deciduous <- rnorm(10000, neon_data_clim$tmp[5], neon_data_clim$tmp[5] * 0.1)
vpd_harv_deciduous <- rnorm(10000, neon_data_clim$vpd[5], neon_data_clim$vpd[5] * 0.1)
par_harv_deciduous <- rnorm(10000, neon_data_clim$par[5], neon_data_clim$par[5] * 0.1)
for(i in 1:length(beta_harv_deciduous)){
  
  optimal_traits <- calc_optimal_vcmax(pathway = 'C3',
                                       deciduous = 'yes',
                                       tg_c = tmp_harv_deciduous[i], 
                                       vpdo = vpd_harv_deciduous[i],
                                       paro = par_harv_deciduous[i],
                                       z = neon_data_clim$z[5],
                                       f = neon_data_clim$f[5],
                                       beta = beta_harv_deciduous[i])
  
  global_optimal_traits_harv_deciduous <- rbind(global_optimal_traits_harv_deciduous, optimal_traits)
  
}

## add lat/lon
global_optimal_traits_harv_deciduous$lat <- neon_data_clim$closest_latitude[5]
global_optimal_traits_harv_deciduous$lon <- neon_data_clim$closest_longitude[5]

### harv_evergreen
global_optimal_traits_harv_evergreen <- data.frame()
beta_harv_evergreen <- exp(rnorm(10000, beta_c3_mean, beta_c3_sd))
tmp_harv_evergreen <- rnorm(10000, neon_data_clim$tmp[5], neon_data_clim$tmp[5] * 0.1)
vpd_harv_evergreen <- rnorm(10000, neon_data_clim$vpd[5], neon_data_clim$vpd[5] * 0.1)
par_harv_evergreen <- rnorm(10000, neon_data_clim$par[5], neon_data_clim$par[5] * 0.1)
for(i in 1:length(beta_harv_evergreen)){
  
  optimal_traits <- calc_optimal_vcmax(pathway = 'C3',
                                       deciduous = 'no',
                                       tg_c = tmp_harv_evergreen[i], 
                                       vpdo = vpd_harv_evergreen[i],
                                       paro = par_harv_evergreen[i],
                                       z = neon_data_clim$z[5],
                                       f = neon_data_clim$f[5],
                                       beta = beta_harv_evergreen[i])
  
  global_optimal_traits_harv_evergreen <- rbind(global_optimal_traits_harv_evergreen, optimal_traits)
  
}

## add lat/lon
global_optimal_traits_harv_evergreen$lat <- neon_data_clim$closest_latitude[5]
global_optimal_traits_harv_evergreen$lon <- neon_data_clim$closest_longitude[5]

### tall_deciduous
global_optimal_traits_tall_deciduous <- data.frame()
beta_tall_deciduous <- exp(rnorm(10000, beta_c3_mean, beta_c3_sd))
tmp_tall_deciduous <- rnorm(10000, neon_data_clim$tmp[15], neon_data_clim$tmp[15] * 0.1)
vpd_tall_deciduous <- rnorm(10000, neon_data_clim$vpd[15], neon_data_clim$vpd[15] * 0.1)
par_tall_deciduous <- rnorm(10000, neon_data_clim$par[15], neon_data_clim$par[15] * 0.1)
for(i in 1:length(beta_tall_deciduous)){
  
  optimal_traits <- calc_optimal_vcmax(pathway = 'C3',
                                       deciduous = 'yes',
                                       tg_c = tmp_tall_deciduous[i], 
                                       vpdo = vpd_tall_deciduous[i],
                                       paro = par_tall_deciduous[i],
                                       z = neon_data_clim$z[15],
                                       f = neon_data_clim$f[15],
                                       beta = beta_tall_deciduous[i])
  
  global_optimal_traits_tall_deciduous <- rbind(global_optimal_traits_tall_deciduous, optimal_traits)
  
}

## add lat/lon
global_optimal_traits_tall_deciduous$lat <- neon_data_clim$closest_latitude[15]
global_optimal_traits_tall_deciduous$lon <- neon_data_clim$closest_longitude[15]

### tall_evergreen
global_optimal_traits_tall_evergreen <- data.frame()
beta_tall_evergreen <- exp(rnorm(10000, beta_c3_mean, beta_c3_sd))
tmp_tall_evergreen <- rnorm(10000, neon_data_clim$tmp[15], neon_data_clim$tmp[15] * 0.1)
vpd_tall_evergreen <- rnorm(10000, neon_data_clim$vpd[15], neon_data_clim$vpd[15] * 0.1)
par_tall_evergreen <- rnorm(10000, neon_data_clim$par[15], neon_data_clim$par[15] * 0.1)
for(i in 1:length(beta_tall_evergreen)){
  
  optimal_traits <- calc_optimal_vcmax(pathway = 'C3',
                                       deciduous = 'no',
                                       tg_c = tmp_tall_evergreen[i], 
                                       vpdo = vpd_tall_evergreen[i],
                                       paro = par_tall_evergreen[i],
                                       z = neon_data_clim$z[15],
                                       f = neon_data_clim$f[15],
                                       beta = beta_tall_evergreen[i])
  
  global_optimal_traits_tall_evergreen <- rbind(global_optimal_traits_tall_evergreen, optimal_traits)
  
}

## add lat/lon
global_optimal_traits_tall_evergreen$lat <- neon_data_clim$closest_latitude[15]
global_optimal_traits_tall_evergreen$lon <- neon_data_clim$closest_longitude[15]

#############################
### global sim PCA ##
#############################

##global pft plots
### combine model outputs with new category
global_optimal_traits_c3_deciduous$pft <- 'c3_deciduous'
global_optimal_traits_c3_evergreen$pft <- 'c3_evergreen'
global_optimal_traits_c4_deciduous$pft <- 'c4_deciduous'

global_optimal_traits_all <- rbind(global_optimal_traits_c3_deciduous, 
                                   global_optimal_traits_c3_evergreen, 
                                   global_optimal_traits_c4_deciduous)

global_optimal_traits_all_select <- as.data.frame(dplyr::select(subset(global_optimal_traits_all, par > 0 & vpd > 0 & tg_c > 0), 
                                                         lma, Anet, wue, gsw, chi, nue, nphoto, narea, nmass, vcmax25, jmax25, rd25,
                                                         vpd, tg_c, par, pft))
global_optimal_traits_all_select_nona <- na.omit(global_optimal_traits_all_select)

### fit pca
global_optimal_traits_all_pca <- prcomp(global_optimal_traits_all_select_nona[,c(1:7,12)], scale = T, center = T)
summary(global_optimal_traits_all_pca)
global_optimal_traits_all_pca$rotation[,1:3]

### plot results
arrow_scale_all <- 5 # Scale factor for arrows and labels to extend from origin
label_scale_all <- 5.5

loadings_pca_all <- as.data.frame(global_optimal_traits_all_pca$rotation[, 1:3]) 
loadings_pca_all$trait <- rownames(loadings_pca_all) ## ad trait column
loadings_pca_all$label_PC1 <- with(loadings_pca_all, PC1 * label_scale_all)
loadings_pca_all$label_PC2 <- with(loadings_pca_all, PC2 * label_scale_all)
loadings_pca_all$label_PC3 <- with(loadings_pca_all, PC3 * label_scale_all)
loadings_pca_all$label_name = c('LMA', 'Anet', 'iWUE', 'gsw', 'χ', 'PNUE', 'Nphoto', 'Rd25')

pca_scores_all <- as.data.frame(global_optimal_traits_all_pca$x) # get scores
pca_scores_all$pft <- global_optimal_traits_all_select_nona$pft

global_optimal_traits_all_pca_plot_PC1PC2 <- ggplot(pca_scores_all, 
                                                    aes(x = PC1, y = PC2, group = pft, color = pft)) +
  theme_minimal(base_size = 14) +
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 14),
        axis.line = element_line(color = "black", linewidth = 0.6),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  geom_point(alpha = 0.1, size = 0.5) +
  scale_color_manual(values = c("cyan", "purple", "orange"),
                     labels = c(expression('C'[3] * ' deciduous'), expression('C'[3] * ' evergreen'), expression('C'[4] * ' deciduous'))) +
  stat_density_2d(aes(fill = after_stat(level)), geom = "polygon", contour = TRUE, alpha = 0.3) +
  scale_fill_viridis_c(option = "turbo") +
  geom_segment(data = loadings_pca_all,
               aes(x = 0, y = 0, xend = PC1 * arrow_scale_all, yend = PC2 * arrow_scale_all, group = NULL, color = NULL),
               arrow = arrow(length = unit(0.25, "cm")), color = "black", linewidth = 0.6) +
  geom_text(data = loadings_pca_all,
            aes(x = label_PC1, y = label_PC2, label = label_name, group = NULL, color = NULL),
            size = 4, fontface = "bold", parse = TRUE, show.legend = FALSE) +
  labs(x = "PC1 (51%)", y = "PC2 (28%)", color = '', tag = '(a)') +
  guides(fill = guide_colorbar(title = "Density level"))

global_optimal_traits_all_pca_plot_PC2PC3 <- ggplot(pca_scores_all, aes(x = PC2, y = PC3, group = pft, color = pft)) +
  theme_minimal(base_size = 14) +
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 14),
        axis.line = element_line(color = "black", linewidth = 0.6),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  geom_point(alpha = 0.1, size = 0.5) +
  scale_color_manual(values = c("cyan", "purple", "orange"),
                     labels = c(expression('C'[3] * ' deciduous'), expression('C'[3] * ' evergreen'), expression('C'[4] * ' deciduous'))) +
  stat_density_2d(aes(fill = ..level..), geom = "polygon", contour = TRUE, alpha = 0.3) +
  scale_fill_viridis_c(option = "turbo") +
  geom_segment(data = loadings_pca_all,
               aes(x = 0, y = 0, xend = PC2 * arrow_scale_all, yend = PC3 * arrow_scale_all, group = NULL, color = NULL),
               arrow = arrow(length = unit(0.25, "cm")), color = "black", linewidth = 0.6) +
  geom_text(data = loadings_pca_all,
            aes(x = label_PC2, y = label_PC3, label = label_name, group = NULL, color = NULL),
            size = 4, fontface = "bold", parse = TRUE, show.legend = FALSE) +
  labs(x = "PC2 (28%)", y = "PC3 (16%)", color = '', tag = '(b)') +
  guides(fill = guide_colorbar(title = "Density level"))

# jpeg('results/plots/global_optimal_traits_all_pca_plot_PC1PC2.jpeg', width = 10, height = 10, units = 'in', res = 600)
# plot(global_optimal_traits_all_pca_plot_PC1PC2)
# dev.off()
# 
# jpeg('results/plots/global_optimal_traits_all_pca_plot_PC2PC3.jpeg', width = 10, height = 10, units = 'in', res = 600)
# plot(global_optimal_traits_all_pca_plot_PC2PC3)
# dev.off()

# jpeg('results/plots/global_optimal_traits_all_pca_plot_PC1PC2_PC2PC3.jpeg', width = 20, height = 10, units = 'in', res = 600)
# multiplot(global_optimal_traits_all_pca_plot_PC1PC2, global_optimal_traits_all_pca_plot_PC2PC3, cols = 2)
# dev.off()

#############################
### environmental predictors of global PCs ##
#############################
pca_scores_all_env <- cbind(global_optimal_traits_all_select_nona, pca_scores_all)
head(pca_scores_all_env)

hist(pca_scores_all_env$PC1)
hist(pca_scores_all_env$PC2)
hist(pca_scores_all_env$PC3)

global_optimal_traits_all_pca$rotation[,1:3]
# PC1: -wue, +chi, -nphoto
# PC2: +Anet, +gsw, +nue
# PC3: +lma, +rd25

PC1_env_lm <- lm(PC1 ~ tg_c + par + vpd + factor(pft), data = pca_scores_all_env)
Anova(PC1_env_lm)
summary(PC1_env_lm) # +tg_c, +par, -vpd, c3e > c3d > c4
PC1_env_lm_relaimp <- calc.relimp(PC1_env_lm, type = c("lmg"), rela = F) # pft and temperature driven
PC1_env_lm_relaimp_data <- c(PC1_env_lm_relaimp$lmg, PC1_env_lm_relaimp$R2)

PC2_env_lm <- lm(PC2 ~ tg_c + par + vpd + factor(pft), data = pca_scores_all_env)
Anova(PC2_env_lm)
summary(PC2_env_lm) # +tg_c, +par, -vpd, c4 > c3d > c3e
PC2_env_lm_relaimp <- calc.relimp(PC2_env_lm, type = c("lmg"), rela = F) # temperature and pft driven, with contributions from par and vpd
PC2_env_lm_relaimp_data <- c(PC2_env_lm_relaimp$lmg, PC2_env_lm_relaimp$R2)

PC3_env_lm <- lm(PC3 ~ tg_c + par + vpd + factor(pft), data = pca_scores_all_env)
Anova(PC3_env_lm)
summary(PC3_env_lm) # -tg_c, +par, +vpd, c3e > c3d > c4
PC3_env_lm_relaimp <- calc.relimp(PC3_env_lm, type = c("lmg"), rela = F) # par and vpd driven
PC3_env_lm_relaimp_data <- c(PC3_env_lm_relaimp$lmg, PC3_env_lm_relaimp$R2)

#### bind all data together and output
PC_env_lm_relaimp_data <- cbind(PC1_env_lm_relaimp_data, PC2_env_lm_relaimp_data, PC3_env_lm_relaimp_data)
colnames(PC_env_lm_relaimp_data) <- c("PC1", "PC2", "PC3")
rownames(PC_env_lm_relaimp_data) <- c("PFT", "Tg", "Ig", "Dg", "R2")
#write.csv(PC_env_lm_relaimp_data, 'results/tables/PC_env_lm_relaimp_data.csv')

#############################
### maps of PCs ##
#############################
# is this even possible given that there are different PC scores for each PFT??

# C3D PC raster

#############################
### global sim trait histograms ##
#############################

global_optimal_traits_all_select_nona$c3c4[global_optimal_traits_all_select_nona$pft=='c3_deciduous'] <- 'c3'
global_optimal_traits_all_select_nona$c3c4[global_optimal_traits_all_select_nona$pft=='c3_evergreen'] <- 'c3'
global_optimal_traits_all_select_nona$c3c4[global_optimal_traits_all_select_nona$pft=='c4_deciduous'] <- 'c4'

global_optimal_traits_all_hist_Anet <- ggplot(data = global_optimal_traits_all_select_nona, 
                                              aes(x = Anet, fill = c3c4)) +
  theme_minimal(base_size = 14) +
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 14),
        axis.line = element_line(color = "black", linewidth = 0.6),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  geom_density(alpha = 0.3) +
  scale_fill_manual(values = c('blue', 'orange'), labels = c(expression('C'[3]), expression('C'[4]))) +
  labs(x = expression('A'[net] * ' (µmol m'^'2'*' s'^'-1'*')'), y = 'Density', fill = '', tag = '(a)')

global_optimal_traits_all_hist_gsw <- ggplot(data = global_optimal_traits_all_select_nona, 
                                              aes(x = gsw, fill = c3c4)) +
  theme_minimal(base_size = 14) +
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 14),
        axis.line = element_line(color = "black", linewidth = 0.6),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  geom_density(alpha = 0.3) +
  scale_fill_manual(values = c('blue', 'orange'), labels = c(expression('C'[3]), expression('C'[4]))) +
  labs(x = expression('g'[sw] * ' (mol m'^'2'*' s'^'-1'*')'), y = 'Density', fill = '', tag = '(b)')

global_optimal_traits_all_hist_wue <- ggplot(data = global_optimal_traits_all_select_nona, 
                                             aes(x = wue, fill = c3c4)) +
  theme_minimal(base_size = 14) +
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 14),
        axis.line = element_line(color = "black", linewidth = 0.6),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  geom_density(alpha = 0.3) +
  scale_fill_manual(values = c('blue', 'orange'), labels = c(expression('C'[3]), expression('C'[4]))) +
  labs(x = expression('iWUE' * ' (µmol mol'^'-1'*')'), y = 'Density', fill = '', tag = '(f)')

global_optimal_traits_all_hist_rd25 <- ggplot(data = global_optimal_traits_all_select_nona, 
                                              aes(x = rd25, fill = c3c4)) +
  theme_minimal(base_size = 14) +
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 14),
        axis.line = element_line(color = "black", linewidth = 0.6),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  geom_density(alpha = 0.3) +
  scale_fill_manual(values = c('blue', 'orange'), labels = c(expression('C'[3]), expression('C'[4]))) +
  labs(x = expression('R'[d25] * ' (µmol m'^'2'*' s'^'-1'*')'), y = 'Density', fill = '', tag = '(d)')

global_optimal_traits_all_hist_chi <- ggplot(data = global_optimal_traits_all_select_nona, 
                                              aes(x = chi, fill = c3c4)) +
  theme_minimal(base_size = 14) +
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 14),
        axis.line = element_line(color = "black", linewidth = 0.6),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  geom_density(alpha = 0.3) +
  scale_fill_manual(values = c('blue', 'orange'), labels = c(expression('C'[3]), expression('C'[4]))) +
  labs(x = expression('χ' * ' (mol mol'^'-1'*')'), y = 'Density', fill = '', tag = '(e)')

global_optimal_traits_all_hist_nphoto <- ggplot(data = global_optimal_traits_all_select_nona, 
                                             aes(x = nphoto, fill = c3c4)) +
  theme_minimal(base_size = 14) +
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 14),
        axis.line = element_line(color = "black", linewidth = 0.6),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  geom_density(alpha = 0.3) +
  scale_fill_manual(values = c('blue', 'orange'), labels = c(expression('C'[3]), expression('C'[4]))) +
  labs(x = expression('N'[photo] * ' (g m'^'-2'*')'), y = 'Density', fill = '', tag = '(c)')

global_optimal_traits_all_hist_lma <- ggplot(data = global_optimal_traits_all_select_nona, 
                                                aes(x = lma, fill = pft)) +
  theme_minimal(base_size = 14) +
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 14),
        axis.line = element_line(color = "black", linewidth = 0.6),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  geom_density(alpha = 0.3) +
  scale_fill_manual(values = c("cyan", "purple", "orange"), 
                    labels = c(expression('C'[3] * ' deciduous'), expression('C'[3] * ' evergreen'), expression('C'[4] * ' deciduous'))) +
  labs(x = expression('LMA' * ' (g m'^'-2'*')'), y = 'Density', fill = '', tag = '(h)')

global_optimal_traits_all_hist_nue <- ggplot(data = global_optimal_traits_all_select_nona, 
                                             aes(x = nue, fill = pft)) +
  theme_minimal(base_size = 14) +
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 14),
        axis.line = element_line(color = "black", linewidth = 0.6),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  geom_density(alpha = 0.3) +
  scale_fill_manual(values = c("cyan", "purple", "orange"), 
                    labels = c(expression('C'[3] * ' deciduous'), expression('C'[3] * ' evergreen'), expression('C'[4] * ' deciduous'))) +
  labs(x = expression('PNUE' * ' (µmol gN'^'-1' * ' s'^'-1' * ')'), y = 'Density', fill = '', tag = '(g)')

# jpeg(filename = "results/plots/global_optimal_traits_hist_all.jpeg", 
#      width = 26, height = 14, units = 'in', res = 600)
# multiplot(global_optimal_traits_all_hist_Anet,
#           global_optimal_traits_all_hist_chi,
#           global_optimal_traits_all_hist_gsw,
#           global_optimal_traits_all_hist_wue,
#           global_optimal_traits_all_hist_nphoto,
#           global_optimal_traits_all_hist_nue,
#           global_optimal_traits_all_hist_rd25,
#           global_optimal_traits_all_hist_lma,
#           cols = 4)
# dev.off()

#############################
### site-level PCA ##
#############################

##site plots
### combine model outputs with new category
global_optimal_traits_cper$site <- 'cper'
global_optimal_traits_konz$site <- 'konz'
global_optimal_traits_scbi$site <- 'scbi'
global_optimal_traits_unde$site <- 'unde'
global_optimal_traits_puum$site <- 'puum'
global_optimal_traits_sjer$site <- 'sjer'
global_optimal_traits_harv_deciduous$site <- 'harv'
global_optimal_traits_harv_evergreen$site <- 'harv'
global_optimal_traits_tall_deciduous$site <- 'tall'
global_optimal_traits_tall_evergreen$site <- 'tall'


global_optimal_traits_sites <- rbind(global_optimal_traits_cper, 
                                     global_optimal_traits_konz,
                                   global_optimal_traits_scbi,
                                   global_optimal_traits_unde,
                                   global_optimal_traits_unde,
                                   global_optimal_traits_puum,
                                   global_optimal_traits_sjer,
                                   global_optimal_traits_harv_deciduous,
                                   global_optimal_traits_harv_evergreen,
                                   global_optimal_traits_tall_deciduous,
                                   global_optimal_traits_tall_evergreen)

global_optimal_traits_sites_select <- as.data.frame(dplyr::select(subset(global_optimal_traits_sites, par > 0 & vpd > 0 & tg_c > 0), 
                                                         lma, Anet, wue, gsw, chi, nue, nphoto, rd25, site))
global_optimal_traits_sites_select_nona <- na.omit(global_optimal_traits_sites_select)

### fit pca
global_optimal_traits_sites_pca <- prcomp(global_optimal_traits_sites_select_nona[,c(1:8)], scale = T, center = T)
summary(global_optimal_traits_sites_pca)
global_optimal_traits_sites_pca$rotation[,1:3]

### plot results
arrow_scale_sites <- 5 # Scale factor for arrows and labels to extend from origin
label_scale_sites <- 5.5

loadings_pca_sites <- as.data.frame(global_optimal_traits_sites_pca$rotation[, 1:3]) 
loadings_pca_sites$trait <- rownames(loadings_pca_sites) ## ad trait column
loadings_pca_sites$label_PC1 <- with(loadings_pca_sites, PC1 * label_scale_sites)
loadings_pca_sites$label_PC2 <- with(loadings_pca_sites, PC2 * label_scale_sites)
loadings_pca_sites$label_PC3 <- with(loadings_pca_sites, PC3 * label_scale_sites)
loadings_pca_sites$label_name = c('LMA', 'Anet', 'iWUE', 'gsw', 'χ', 'PNUE', 'Nphoto', 'Rd25')

pca_scores_sites <- as.data.frame(global_optimal_traits_sites_pca$x) # get scores
pca_scores_sites$site <- global_optimal_traits_sites_select_nona$site

#### get average PC1 and PC2 for plotting purposes
pca_scores_sites_group_by <- group_by(pca_scores_sites, site)
pca_scores_sites_summarise <- summarise(pca_scores_sites_group_by, mean_PC1 = mean(PC1), mean_PC2 = mean(PC2), mean_PC3 = mean(PC3))


global_optimal_traits_sites_pca_plot_PC1PC2 <- ggplot(pca_scores_sites, 
                                                    aes(x = PC1, y = PC2, group = site, color = site)) +
  theme_minimal(base_size = 14) +
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 14),
        axis.line = element_line(color = "black", linewidth = 0.6),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  geom_point(alpha = 0.1, size = 0.5) +
  scale_color_manual(values = c("orange1", "red1", "orange4", "purple1", "cyan1", "purple4", "red4", "cyan4"),
                     labels = c(expression("CPER (C"[4]*"D)"), expression("HARV (C"[3]*"M)"), expression("KONZ (C"[4]*"D)"), expression("PUUM (C"[3]*"E)"),
                                expression("SCBI (C"[3]*"D)"), expression("SJER (C"[3]*"E)"), expression("TALL (C"[3]*"M)"), expression("UNDE (C"[3]*"D)"))) +
  stat_density_2d(aes(fill = after_stat(level)), geom = "polygon", contour = TRUE, alpha = 0.3) +
  scale_fill_viridis_c(option = "turbo") +
  geom_segment(data = loadings_pca_sites,
               aes(x = 0, y = 0, xend = PC1 * arrow_scale_sites, yend = PC2 * arrow_scale_sites, group = NULL, color = NULL),
               arrow = arrow(length = unit(0.25, "cm")), color = "black", linewidth = 0.6) +
  geom_text(data = loadings_pca_sites,
            aes(x = label_PC1, y = label_PC2, label = label_name, group = NULL, color = NULL),
            size = 4, fontface = "bold", parse = TRUE, show.legend = FALSE) +
  labs(x = "PC1 (55%)", y = "PC2 (20%)", color = 'Site', tag = '(a)') +
  geom_point(data = pca_scores_sites_summarise, aes(x = mean_PC1, y = mean_PC2), size = 5, shape = 1, stroke = 3, show.legend = FALSE) +
  geom_point(data = pca_scores_sites_summarise, aes(x = mean_PC1, y = mean_PC2), size = 2, shape = 16, color = 'black', show.legend = FALSE) +
  guides(fill = guide_colorbar(title = "Density level"))

global_optimal_traits_sites_pca_plot_PC2PC3 <- ggplot(pca_scores_sites, 
                                                      aes(x = PC2, y = PC3, group = site, color = site)) +
  theme_minimal(base_size = 14) +
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 14),
        axis.line = element_line(color = "black", linewidth = 0.6),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  geom_point(alpha = 0.1, size = 0.5) +
  scale_color_manual(values = c("orange1", "red1", "orange4", "purple1", "cyan1", "purple4", "red4", "cyan4"),
                     labels = c(expression("CPER (C"[4]*"D)"), expression("HARV (C"[3]*"M)"), expression("KONZ (C"[4]*"D)"), expression("PUUM (C"[3]*"E)"),
                                expression("SCBI (C"[3]*"D)"), expression("SJER (C"[3]*"E)"), expression("TALL (C"[3]*"M)"), expression("UNDE (C"[3]*"D)"))) +
  stat_density_2d(aes(fill = after_stat(level)), geom = "polygon", contour = TRUE, alpha = 0.3) +
  scale_fill_viridis_c(option = "turbo") +
  geom_segment(data = loadings_pca_sites,
               aes(x = 0, y = 0, xend = PC2 * arrow_scale_sites, yend = PC3 * arrow_scale_sites, group = NULL, color = NULL),
               arrow = arrow(length = unit(0.25, "cm")), color = "black", linewidth = 0.6) +
  geom_text(data = loadings_pca_sites,
            aes(x = label_PC2, y = label_PC3, label = label_name, group = NULL, color = NULL),
            size = 4, fontface = "bold", parse = TRUE, show.legend = FALSE) +
  labs(x = "PC2 (20%)", y = "PC3 (16%)", color = 'Site', tag = '(b)') +
  geom_point(data = pca_scores_sites_summarise, aes(x = mean_PC2, y = mean_PC3), size = 5, shape = 1, stroke = 3, show.legend = FALSE) +
  geom_point(data = pca_scores_sites_summarise, aes(x = mean_PC2, y = mean_PC3), size = 2, shape = 16, color = 'black', show.legend = FALSE) +
  guides(fill = guide_colorbar(title = "Density level"))

# jpeg('results/plots/global_optimal_traits_sites_pca_plot_PC1PC2.jpeg', width = 10, height = 10, units = 'in', res = 600)
# plot(global_optimal_traits_sites_pca_plot_PC1PC2)
# dev.off()

# jpeg('results/plots/global_optimal_traits_sites_pca_plot_PC1PC2_PC2PC3.jpeg', width = 20, height = 10, units = 'in', res = 600)
# multiplot(global_optimal_traits_sites_pca_plot_PC1PC2, global_optimal_traits_sites_pca_plot_PC2PC3, cols = 2)
# dev.off()

#############################
### site-level histograms ##
#############################
global_optimal_traits_sites_hist_Anet <- ggplot(data = global_optimal_traits_sites_select_nona, 
                                              aes(x = Anet, fill = site)) +
  theme_minimal(base_size = 14) +
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 14),
        axis.line = element_line(color = "black", linewidth = 0.6),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  geom_density(alpha = 0.3, linewidth = 1) +
  scale_fill_manual(values = c("orange1", "red1", "orange4", "purple1", "cyan1", "purple4", "red4", "cyan4"),
                    labels = c(expression("CPER (C"[4]*"D)"), expression("HARV (C"[3]*"M)"), expression("KONZ (C"[4]*"D)"), expression("PUUM (C"[3]*"E)"),
                               expression("SCBI (C"[3]*"D)"), expression("SJER (C"[3]*"E)"), expression("TALL (C"[3]*"M)"), expression("UNDE (C"[3]*"D)"))) +
  labs(x = expression('A'[net] * ' (µmol m'^'2'*' s'^'-1'*')'), y = 'Density', fill = '', tag = '(a)')

global_optimal_traits_sites_hist_gsw <- ggplot(data = global_optimal_traits_sites_select_nona, 
                                             aes(x = gsw, fill = site)) +
  theme_minimal(base_size = 14) +
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 14),
        axis.line = element_line(color = "black", linewidth = 0.6),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  geom_density(alpha = 0.3, linewidth = 1) +
  scale_fill_manual(values = c("orange1", "red1", "orange4", "purple1", "cyan1", "purple4", "red4", "cyan4"),
                    labels = c(expression("CPER (C"[4]*"D)"), expression("HARV (C"[3]*"M)"), expression("KONZ (C"[4]*"D)"), expression("PUUM (C"[3]*"E)"),
                               expression("SCBI (C"[3]*"D)"), expression("SJER (C"[3]*"E)"), expression("TALL (C"[3]*"M)"), expression("UNDE (C"[3]*"D)"))) +
  labs(x = expression('g'[sw] * ' (mol m'^'2'*' s'^'-1'*')'), y = 'Density', fill = '', tag = '(b)')

global_optimal_traits_sites_hist_wue <- ggplot(data = global_optimal_traits_sites_select_nona, 
                                             aes(x = wue, fill = site)) +
  theme_minimal(base_size = 14) +
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 14),
        axis.line = element_line(color = "black", linewidth = 0.6),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  geom_density(alpha = 0.3, linewidth = 1) +
  scale_fill_manual(values = c("orange1", "red1", "orange4", "purple1", "cyan1", "purple4", "red4", "cyan4"),
                    labels = c(expression("CPER (C"[4]*"D)"), expression("HARV (C"[3]*"M)"), expression("KONZ (C"[4]*"D)"), expression("PUUM (C"[3]*"E)"),
                               expression("SCBI (C"[3]*"D)"), expression("SJER (C"[3]*"E)"), expression("TALL (C"[3]*"M)"), expression("UNDE (C"[3]*"D)"))) +
  labs(x = expression('iWUE' * ' (µmol mol'^'-1'*')'), y = 'Density', fill = '', tag = '(f)')

global_optimal_traits_sites_hist_rd25 <- ggplot(data = global_optimal_traits_sites_select_nona, 
                                              aes(x = rd25, fill = site)) +
  theme_minimal(base_size = 14) +
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 14),
        axis.line = element_line(color = "black", linewidth = 0.6),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  geom_density(alpha = 0.3, linewidth = 1) +
  scale_fill_manual(values = c("orange1", "red1", "orange4", "purple1", "cyan1", "purple4", "red4", "cyan4"),
                    labels = c(expression("CPER (C"[4]*"D)"), expression("HARV (C"[3]*"M)"), expression("KONZ (C"[4]*"D)"), expression("PUUM (C"[3]*"E)"),
                               expression("SCBI (C"[3]*"D)"), expression("SJER (C"[3]*"E)"), expression("TALL (C"[3]*"M)"), expression("UNDE (C"[3]*"D)"))) +
  labs(x = expression('R'[d25] * ' (µmol m'^'2'*' s'^'-1'*')'), y = 'Density', fill = '', tag = '(d)')

global_optimal_traits_sites_hist_chi <- ggplot(data = global_optimal_traits_sites_select_nona, 
                                             aes(x = chi, fill = site)) +
  theme_minimal(base_size = 14) +
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 14),
        axis.line = element_line(color = "black", linewidth = 0.6),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  geom_density(alpha = 0.3, linewidth = 1) +
  scale_fill_manual(values = c("orange1", "red1", "orange4", "purple1", "cyan1", "purple4", "red4", "cyan4"),
                    labels = c(expression("CPER (C"[4]*"D)"), expression("HARV (C"[3]*"M)"), expression("KONZ (C"[4]*"D)"), expression("PUUM (C"[3]*"E)"),
                               expression("SCBI (C"[3]*"D)"), expression("SJER (C"[3]*"E)"), expression("TALL (C"[3]*"M)"), expression("UNDE (C"[3]*"D)"))) +
  labs(x = expression('χ' * ' (mol mol'^'-1'*')'), y = 'Density', fill = '', tag = '(e)')

global_optimal_traits_sites_hist_nphoto <- ggplot(data = global_optimal_traits_sites_select_nona, 
                                                aes(x = nphoto, fill = site)) +
  theme_minimal(base_size = 14) +
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 14),
        axis.line = element_line(color = "black", linewidth = 0.6),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  geom_density(alpha = 0.3, linewidth = 1) +
  scale_fill_manual(values = c("orange1", "red1", "orange4", "purple1", "cyan1", "purple4", "red4", "cyan4"),
                    labels = c(expression("CPER (C"[4]*"D)"), expression("HARV (C"[3]*"M)"), expression("KONZ (C"[4]*"D)"), expression("PUUM (C"[3]*"E)"),
                               expression("SCBI (C"[3]*"D)"), expression("SJER (C"[3]*"E)"), expression("TALL (C"[3]*"M)"), expression("UNDE (C"[3]*"D)"))) +
  labs(x = expression('N'[photo] * ' (g m'^'-2'*')'), y = 'Density', fill = '', tag = '(c)')

global_optimal_traits_sites_hist_lma <- ggplot(data = global_optimal_traits_sites_select_nona, 
                                             aes(x = lma, fill = site)) +
  theme_minimal(base_size = 14) +
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 14),
        axis.line = element_line(color = "black", linewidth = 0.6),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  geom_density(alpha = 0.3, linewidth = 1) +
  scale_fill_manual(values = c("orange1", "red1", "orange4", "purple1", "cyan1", "purple4", "red4", "cyan4"),
                     labels = c(expression("CPER (C"[4]*"D)"), expression("HARV (C"[3]*"M)"), expression("KONZ (C"[4]*"D)"), expression("PUUM (C"[3]*"E)"),
                                expression("SCBI (C"[3]*"D)"), expression("SJER (C"[3]*"E)"), expression("TALL (C"[3]*"M)"), expression("UNDE (C"[3]*"D)"))) +
  labs(x = expression('LMA' * ' (g m'^'-2'*')'), y = 'Density', fill = '', tag = '(h)')

global_optimal_traits_sites_hist_nue <- ggplot(data = global_optimal_traits_sites_select_nona, 
                                             aes(x = nue, fill = site)) +
  theme_minimal(base_size = 14) +
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 14),
        axis.line = element_line(color = "black", linewidth = 0.6),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  geom_density(alpha = 0.3, linewidth = 1) +
  scale_fill_manual(values = c("orange1", "red1", "orange4", "purple1", "cyan1", "purple4", "red4", "cyan4"),
                    labels = c(expression("CPER (C"[4]*"D)"), expression("HARV (C"[3]*"M)"), expression("KONZ (C"[4]*"D)"), expression("PUUM (C"[3]*"E)"),
                               expression("SCBI (C"[3]*"D)"), expression("SJER (C"[3]*"E)"), expression("TALL (C"[3]*"M)"), expression("UNDE (C"[3]*"D)"))) +
  labs(x = expression('PNUE' * ' (µmol gN'^'-1' * ' s'^'-1' * ')'), y = 'Density', fill = '', tag = '(g)')

# jpeg(filename = "results/plots/global_optimal_traits_hist_sites.jpeg", 
#     width = 26, height = 14, units = 'in', res = 600)
# multiplot(global_optimal_traits_sites_hist_Anet,
#          global_optimal_traits_sites_hist_chi,
#          global_optimal_traits_sites_hist_gsw,
#          global_optimal_traits_sites_hist_wue,
#          global_optimal_traits_sites_hist_nphoto,
#          global_optimal_traits_sites_hist_nue,
#          global_optimal_traits_sites_hist_rd25,
#          global_optimal_traits_sites_hist_lma,
#          cols = 4)
# dev.off()


#############################
### map of sites ##
#############################

neon_data_sub <- subset(neon_data, site_id %in% c('CPER', 'HARV', 'KONZ', 'PUUM', 'SCBI', 'SJER', 'TALL', 'UNDE'))
neon_data_sub$site_id
neon_data_sub_latlon <- select(neon_data_sub, longitude, latitude)
neon_data_sub_latlon_trans <- usmap_transform(neon_data_sub_latlon,
                                              input_names = c("longitude", "latitude"))

site_map <- plot_usmap("states", color = 'grey') +
  geom_sf(data = neon_data_sub_latlon_trans, 
             #color = c("orange1", "red1", "orange4", "purple1", "cyan1", "purple4", "red4", "cyan4"),
             color = c("orange1", "red1", "orange4", "cyan1", "purple4", "red4", "cyan4", "purple1"),
             size = 3,
          shape = 1, stroke = 2) +
  geom_sf(data = neon_data_sub_latlon_trans, 
          color = 'black',
          size = 1) +
  geom_sf_text(data = neon_data_sub_latlon_trans,
                aes(label = c('CPER', 'HARV', 'KONZ', 'SCBI', 'SJER', 'TALL', 'UNDE', 'PUUM')),
               position = position_nudge(y=-100000), size = 5, fontface = 'bold')
  
# jpeg(filename = "results/plots/global_optimal_traits_site_map.jpeg", width = 10, height = 10, units = 'in', res = 600)
# plot(site_map)
# dev.off()

#############################
### map of global growing season climate ##
#############################
# map to show spatial variability in tmp, vpd, and par
global_data_veg_select <- subset(global_data_veg, tmp > 0 & vpd >0 &par >0)
hist(global_data_veg_select$tmp)
hist(global_data_veg_select$vpd)
hist(global_data_veg_select$par)

# jpeg(filename = "results/plots/cru_temperature_map.jpeg", width = 10, height = 10, units = 'in', res = 600)
# temperature_raster <- rasterFromXYZ(global_data_veg_select[,c(1, 2, 4)])
# plot(temperature_raster, col=rev(brewer.pal(9,'Spectral')), breaks=seq(10,40,3),
#      legend.args=list(text=expression('T'[g])))
# maps::map('world',col='black',fill=F, add = T) #see line 180 of global leaf evaluation code (or see below)
# dev.off()

# jpeg(filename = "results/plots/cru_vpd_map.jpeg", width = 10, height = 10, units = 'in', res = 600)
# vpd_raster <- rasterFromXYZ(global_data_veg_select[,c(1, 2, 5)])
# plot(vpd_raster, col=rev(brewer.pal(9,'Spectral')), breaks=seq(0,3,0.5),
#      legend.args=list(text=expression('D'[g])))
# maps::map('world',col='black',fill=F, add = T) #see line 180 of global leaf evaluation code (or see below)
# dev.off()

# jpeg(filename = "results/plots/cru_par_map.jpeg", width = 10, height = 10, units = 'in', res = 600)
# par_raster <- rasterFromXYZ(global_data_veg_select[,c(1, 2, 3)])
# plot(par_raster, col=rev(brewer.pal(9,'Spectral')), breaks=seq(400,1400,100),
#      legend.args=list(text=expression('D'[g])))
# maps::map('world',col='black',fill=F, add = T) #see line 180 of global leaf evaluation code (or see below)
# dev.off()

# par(mfrow=c(1,1), oma=c(4,4,1,2), mar=c(1,1,1,1))
# plot(Vcmax_globe_ras_mod, col=cols, breaks = seq(0, 140, 5), cex.axis=1.5, yaxt = 'n', xaxt = 'n', lab.breaks = seq(0, 140, 5), ylim = c(-90, 90), legend.args=list(text=expression(italic('V'*"'")[cmax]*' (µmol m'^'-2'*' s'^'-1'*')'), line = 4, side = 4, las = 3, cex = 1.5), legend = T, xlim = c(-180, 180), axis.args = arg)
# map('world',col='black',fill=F, add = T, ylim = c(-180, 180))
# axis(2, at = seq(-90, 90, 30), labels = T, cex.axis = 1.5)
# axis(1, at = seq(-180, 180, 90), labels = T, cex.axis = 1.5)