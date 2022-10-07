# Data management
library(xlsx)
library(readxl)
library(tidyr)
library(dplyr)

# Analysis
library(sdmTMB)
library(ape)

# Mapping
library(sf)

# Visualization
library(visreg)
library(ggplot2)

setwd("C://KDale/Projects/Phenology/Data/")

all_tows_roms <- read.csv("AllTows_200nm.csv") %>% subset(., year >= 1995 & year <= 2019)

northAmerica <- read_sf("C://KDale/GIS/North_South_America.shp")
northAmerica <- sf::st_transform(northAmerica, crs = 5070)

# GET SPECIES DATA -------------------------------------------------------------
getspeciesData <- function(species) {
  speciesData <- read.csv("AllCruises_Combined_200nm.csv") %>%
    subset(scientific_name == species & year >= 1995 & year <= 2019) %>%
    mutate(., presence = 1)
  
  # Get maximum extents of positive tows
  maxLat = max(speciesData$latitude) 
  minLat = min(speciesData$latitude)
  maxLon = max(speciesData$longitude)
  minLon = min(speciesData$longitude)
  
  # Merge with tow data
  speciesData <- merge(all_tows_roms, speciesData[c("towID", "scientific_name", "density_anomalies", "presence")], all.x = TRUE)
  
  # Add the minimum to all density anomaly values plus a very small value to avoid zeroes
  speciesData$density_anomaly_positive = speciesData$density_anomalies + abs(min(speciesData$density_anomalies, na.rm = T)) + 0.01
  
  # Replace all NA density anomaly values with zero
  speciesData <-  mutate(speciesData, presence = replace_na(presence, 0),
                         density_anomaly_positive = replace_na(density_anomaly_positive, 0), 
                        density_anomalies = replace_na(density_anomalies, 0),
                         scientific_name = species) %>% 
    subset(., !is.na(sst_roms) & !is.na(ssh_roms) & !is.na(salinity_roms) & # Remove any tows without ROMS data
             latitude <= maxLat & latitude >= minLat & # Subset to within area bounded by positive tows
             longitude >= minLon & longitude <= maxLon) %>%
    mutate(., sst_scaled = scale(sst_roms)[,1], ssh_scaled = scale(ssh_roms)[,1], salinity_scaled = scale(salinity_roms)[,1]) # Center and scale enviro data
  
  # Add UTM in km coordinates
  #speciesData <- sdmTMB::add_utm_columns(dat = speciesData, ll_names = c("longitude", "latitude"))
  
  # Add Albert Equal Area coordinates
  coords= cbind(speciesData$longitude, speciesData$latitude)
  scale = 1000
  albert_equal_area <- sf::sf_project(from = "EPSG:4326", to = 'EPSG:5070', pts = coords)/scale
  speciesData$X = albert_equal_area[,1]
  speciesData$Y = albert_equal_area[,2]
  
  # Add time block
  speciesData$timeblock = 0
  for (i in 1:nrow(speciesData)) {
    year = speciesData$year[i] 
    
    if(year < 2000) {
      speciesData$timeblock[i] = "1995-1999"
    } else if (year < 2005 ) {
      speciesData$timeblock[i] = "2000-2004"
    } else if (year < 2010) {
      speciesData$timeblock[i] = "2005-2009"
    } else if (year < 2015) {
      speciesData$timeblock[i] = "2010-2014"
    } else if (year < 2020) {
      speciesData$timeblock[i] = "2015-2019"
    } 
  }
  
  # Add latitudinal region
  speciesData$region = 0
  for (i in 1:nrow(speciesData)) {
    latitude = speciesData$latitude[i] 
    
    if(latitude < 34.5) {
      speciesData$region[i] = "Southern CCE"
    } else if (latitude < 42 ) {
      speciesData$region[i] = "Central CCE"
    } else if (latitude < 48.3) {
      speciesData$region[i] = "OR/WA"
    } else if (latitude < 54.4) {
      speciesData$region[i] = "British Columbia"
    } else {
      speciesData$region[i] = "Gulf of Alaska"
    } 
  }
  
  speciesData$region <- factor(speciesData$region, levels = c("Southern CCE", "Central CCE", "OR/WA", "British Columbia", "Gulf of Alaska"))
  
  return(speciesData)
}

# RUN MODELS --------------------------------------------------------------------------------
runModel <- function(data, type) {
  
  if (type == "presence") {
    # Presence
    # Step across months
    fit <- sdmTMB(
      presence ~ 0 + as.factor(timeblock) + s(ssh_scaled) + s(sst_scaled) + s(salinity_scaled),
      data = data,
      mesh = mesh,
      family = binomial(link = "logit"),
      spatial = "on", spatiotemporal = "RW",
      silent = FALSE,
      time = "month")
    
    save(fit, data, mesh, file=paste0("C://KDale/Projects/Phenology/Results/",species,"_month_presence.rdata"))
    
  } else if (type == "density") {
    # Density anomalies
    # Step across months
    fit <- sdmTMB(
      density_anomaly_positive ~ 0 + as.factor(timeblock) + s(ssh_scaled) + s(sst_scaled) + s(salinity_scaled),
      data = data,
      mesh = mesh,
      family = tweedie(link = "log"),
      spatial = "on", spatiotemporal = "RW",
      silent = FALSE,
      time = "month")
    
    save(fit, data, mesh, file="C://KDale/Projects/Phenology/Results/",species,"_month_density.rdata")
  }
  
  return(fit)
}

# MAKE GRID -------------------------------------------------------------------
makeGrid <- function(data) {
  
  # Load land shapefiles
  northSouthAmerica <- st_read("C://KDale/GIS/North_South_America.shp") %>%
    st_union() %>% st_make_valid()
  
  # Convert species data to multipoint, and get positive stations for study area creation
  # Casting to multipoint is required for creating convex hull
  data <- st_as_sf(data, coords = c("longitude", "latitude"), crs = 4326) %>% 
    st_cast(., to = "MULTIPOINT") %>%
    st_combine(.)
  
  # Create study area polygon
  studyArea <- st_convex_hull(data) %>%
    st_buffer(., dist = 1000, singleSide = FALSE)
  
  # Create grid
  
  hex <- st_make_grid(studyArea, crs = 4326, cellsize = 1.5, what = "polygons", square = FALSE) %>%
    st_make_valid() %>%
    st_difference(., northSouthAmerica) %>% # remove hexes that are on land
    st_intersection(., studyArea) %>% 
    st_sf(.) %>%
    mutate(., GridID = 1:length(.))
  
  # Add hex ID to dataset
  data.hex <- st_join(data, hex)
  
  # Print hex with catch info
  print(ggplot() +
          geom_sf(data = hex, fill = "gray95") +
          geom_sf(
            data = data,
            mapping = aes(col = density_anomalies),
            size = 3,
            alpha = 0.4) +
          scale_color_gradient("Standardized # of larvae", low = "yellow3", high = "red3") +
          theme_classic())
  
}


#INLA MESH TO SF ---------------------------------------------------------------
# Convert INLA mesh to sf (experimental function from fmesher)
inlaMesh_to_sf <- function(x, ...) {
  stopifnot(inherits(x, "inla.mesh"))
  geom <- sf::st_geometry(
    sf::st_multipolygon(
      lapply(
        seq_len(nrow(x$graph$tv)),
        function(k) {
          list(x$loc[x$graph$tv[k, c(1, 2, 3, 1)], , drop = FALSE])
        }
      ),
      dim = "XYZ"
    )
  )
  #sf::st_crs(geom) <- fm_as_sf_crs(x$crs)
  geom
}
# FREQUENCY OF SAMPLING -------------------------------------------------------
frequencyOfSampling <- function() {
  
  # Frequency of sampling across months/years by program
  freq_months <- speciesData %>% group_by_at(., c("timeblock", "month", "program")) %>%
    summarize(n = n())
  
  print(ggplot(freq_months) +
    geom_tile(aes(month, timeblock, fill = n)) +
    facet_wrap(vars(program)) +
    scale_fill_gradient2("Number \n of tows", low = "white", mid = "cornflowerblue", high = "darkblue") +
    scale_x_continuous(n.breaks = 12) +
    labs(x = "Month", y = "Program") +
    coord_cartesian(expand = FALSE) +
    theme_classic())
  
  # Frequency of sampling across months/years by region
  freq_months <- speciesData %>% group_by_at(., c("timeblock", "month", "region")) %>%
    summarize(n = n())

  print(ggplot(freq_months) +
    geom_tile(aes(month, timeblock, fill = n)) +
    facet_wrap(vars(region)) +
    scale_fill_gradient2("Number \n of tows", low = "white", mid = "cornflowerblue", high = "darkblue") +
    scale_x_continuous(n.breaks = 12) +
    labs(x = "Month", y = "Region") +
    coord_cartesian(expand = FALSE) +
    theme_classic())
}

# RUN FUNCTIONS -----------------------------------------------------------------

## Get species data -----
species = "Tarletonbeania crenularis"
speciesData <- getspeciesData(species)
speciesData <- subset(speciesData, program != "RREAS" & program != "PRS_larvae" & program != "PRS_juveniles")

## Make mesh -------
mesh <- make_mesh(speciesData, xy_cols = c("X",  "Y"), n_knots = 200, type= "cutoff_search")
speciesData$meshNum <- mesh$mesh$idx$loc

# Extract mesh to do things with it -- still working on this
# mesh.sf <- inlaMesh_to_sf(mesh$mesh)
# mesh.df <- st_make_valid(mesh.sf)
# %>% st_sf(., crs = "EPSG:5070")  %>% mutate(., GridID = 1:length(.))

## Plot mesh w/ empirical data -----
ggplot() +
  inlabru::gg(mesh$mesh) +
  geom_point(data = subset(speciesData, density_anomaly_positive = 0), aes(x = X, y = Y), col = "gray50") +
  #geom_point(data = subset(speciesData, density_anomaly_positive > 0), aes(x = X, y = Y, col = density_anomaly_positive)) +
  scale_color_gradient(low = "darkblue", high = "cornflowerblue") +
  theme_classic() +
  coord_equal()

## Run model----
fit <- runModel(data = speciesData, type = "presence")

## Reload model run -----
load(file = paste0("C://KDale/Projects/Phenology/Results/",species,"_month_presence.rdata"))

# extract parameters as a dataframe
# range: A derived parameter that defines the distance at which 2 points are effectively independent
# (actually about 13% correlated). If the share_range argument is changed to FALSE then the spatial
# and spatiotemporal ranges will be unique, otherwise the default is for both to share the same range.
# phi: Observation error scale parameter (e.g., SD in Gaussian).
# sigma_O: SD of the spatial process ("Omega").
# sigma_E: SD of the spatiotemporal process ("Epsilon").
# tweedie_p: Tweedie p (power) parameter; between 1 and 2.

## Run basic checks -----
sanity(fit)

fit

# View confidence intervals and extract parameters as a dataframe
tidy(fit, effects = "ran_pars", conf.int = TRUE)

# Randomized quantile residuals - quick
speciesData$resids <- residuals(fit)
qqnorm(speciesData$resids)
qqline(speciesData$resids)

# RQRs - MCMC
mcmc_res <- residuals(fit, type = "mle-mcmc", mcmc_iter = 201, mcmc_warmup = 200)
qqnorm(mcmc_res);qqline(mcmc_res)

# Plot smooother on variables in link space with randomized quantile partial residuals
# takes a few min to run
visreg(fit, xvar = "sst_scaled")
visreg(fit, xvar = "ssh_scaled")
visreg(fit, xvar = "salinity_scaled")

# Plot the response scale
visreg::visreg(fit, xvar = "sst_scaled", scale = "response")
visreg::visreg(fit, xvar = "ssh_scaled", scale = "response")
visreg::visreg(fit, xvar = "salinity_scaled", scale = "response")

## Predict on original data -----
p.obj <- predict(fit, return_tmb_object = T)
p <- predict(fit)

plot(x = speciesData$sst_scaled, y = speciesData$p)

# Quantifying predictive performance for presence/absence model
speciesData$p <- predict(fit)$est
rocr <- ROCR::prediction(exp(speciesData$p), speciesData$presence)
ROCR::performance(rocr, measure = "auc")@y.values[[1]] # Values near 0.5 ~ random, want close to 1

## Center of gravity --------
# works on object returned from predict()
cog = get_cog(p, format = "wide") # Does not work on prediction without an sdmTMB object returned
ggplot(cog, aes(est_x, est_y, colour = year)) +
  geom_pointrange(aes(xmin = lwr_x, xmax = upr_x)) +
  geom_pointrange(aes(ymin = lwr_y, ymax = upr_y)) +
  scale_colour_viridis_c()

## Plot in natural space ---------
# Hint: inverse link function also available in `fit$family$linkinv()`
jpeg(filename = "C://KDale/Projects/Phenology/Figures/Empirical.jpg", width = 7, height = 4, units = "in", res = 500)
ggplot(northAmerica) +
  geom_sf() +
  facet_wrap(~month, nrow = 2) +
  geom_point(data = p, mapping = aes(x = (X*1000), y = (Y*1000), col = plogis(est))) +
  scale_fill_viridis_c() +
  #geom_point(data = speciesData, aes(X*1000, Y*1000, col = presence), pch = 20, inherit.aes = FALSE) +
  scale_color_steps("Predicted presence", low = "gray70", high =  "firebrick") +
  xlim(min(speciesData$X)*1000-1000, max(speciesData$X)*1000+1000) +
  ylim(min(speciesData$Y)*1000-1000, max(speciesData$Y)*1000+1000) +
  labs(x = "Longitude", y = "Latitude") +
  theme_classic()
dev.off()

ggplot() +
  scale_fill_viridis_c() +
  geom_point(data = p, mapping = aes(x = (X*1000), y = (Y*1000), color = plogis(est))) 
  # trim extreme high values to make spatial variation more visible
  #facet_wrap(~timeblock)

# CENTRAL TENDENCY ----------------------------------------------------------

# Average densities
speciesData.timeblocks = group_by_at(speciesData, c("timeblock", "month", "region")) %>%
  summarize(avg_density = mean(density_anomaly_positive), avg_presence = mean(presence)) %>% 
  mutate(density_x_month = avg_density * month, presence_x_month = avg_presence * month)

# Calculate CT for each block
# sum of (month # * mean abundance in month)/(sum of mean abundances across all months)
speciesData.ct <- group_by_at(speciesData.timeblocks, c("timeblock", "region")) %>% 
  summarize(ct_density = sum(density_x_month)/sum(avg_density),
            ct_presence = sum(presence_x_month)/sum(avg_presence))

ggplot(speciesData.timeblocks) +
  geom_line(mapping = aes(x = month, y = avg_density, fill = timeblock, color = timeblock), lwd = 1) +
  facet_wrap(vars(region)) +
  scale_color_brewer("Time period", palette = "Greens" ) +
  scale_fill_brewer("Time period", palette = "Greens" ) +
  scale_x_continuous(n.breaks = 12) +
  labs(x = "Month", y = "Average density anomaly \n (positive transformation)") +
  theme_classic() 

# Frequency of sampling plots
frequencyOfSampling()

# SANDBOX --------------------------------------------------------------------

# Test if a spatial model should be used (Moran's I)
inv_dists <- as.matrix(dist(speciesData[,c("X","Y")]))
diag(inv_dists) <- 0
Moran.I(speciesData$density_anomaly_positive, inv_dists)



