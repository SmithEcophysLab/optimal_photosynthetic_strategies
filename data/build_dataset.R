# build_dataset

## combine and build cru datasets

cru_latlon <- read.csv('cru_growingseason/cru_par_climExtract_growingseason_globe.csv')[,2:3]
cru_par <- read.csv('cru_growingseason/cru_par_clim_growingseason.csv')[,2]
cru_tmp <- read.csv('cru_growingseason/cru_tmp_clim_growingseason.csv')[,2]
cru_vpd <- read.csv('cru_growingseason/cru_vpd_clim_growingseason.csv')[,2]
cru_f <- read.csv('cru_growingseason/cru_f_clim_growingseason.csv')[,2]/12

cru_all <- cbind(cru_latlon, cru_par, cru_tmp, cru_vpd, cru_f)

colnames(cru_all) <- c('lon', 'lat', 'par', 'tmp', 'vpd', 'f')

write.csv(cru_all, 'cru_growingseason/cru_growingseason.csv')
