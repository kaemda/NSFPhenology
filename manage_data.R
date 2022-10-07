# data management
library(xlsx)
library(readxl) 
library(tidyr)
library(dplyr)
library(lubridate)
library(stringr)
library(RANN)

# data access
library(ncdf4)
library(rerddap)
library(rerddapXtracto)

# analysis
library(qdap)

# mapping and plotting
library(sf)
library(ggplot2)

setwd("C://KDale/Projects/Phenology/Data/")

speciesInfo <- read_xlsx(path = "species_of_interest.xlsx", sheet = 1)
speciesNames <- as.vector(speciesInfo$scientific_name)

colors = c("firebrick3",  "coral3", "darkgoldenrod2", "darkseagreen3","cornflowerblue", "deepskyblue3",  "darkslateblue","darkslategray")

# ERDDAP function ---------------------------------------------------------------
getERDDAP <- function(name) {
  
  data <- tabledap(x = name, url = "https://coastwatch.pfeg.noaa.gov/erddap") %>% 
    separate(., col = time, into = c("date", "time"), sep = c("T")) %>% 
    mutate(., .after = "date", year = year(date), month = month(date), day = day(date))
  
  data$cruise <- as.character(data$cruise)
  data$date <- as.Date(data$date)
  
  return(data)
}

# IMECOCAL ---------------------------------------------------------------------
getIMECOCAL <- function() {
  species_codes <- read.csv("Species_codes.csv", as.is = TRUE)
  
  sheet_names = excel_sheets("OriginalDatasets/Imecocal all species1.xlsx")
  nsheets = length(sheet_names)
  
  for(i in 1:nsheets) {
    
    sheet <- read_xlsx("OriginalDatasets/Imecocal all species1.xlsx", sheet = i)
    
    # Go from wide to long
    sheet_long <- pivot_longer(sheet, cols = -colnames(sheet[1:10]), names_to = "species_code", values_to = "larvae_10m2") %>% 
      subset(., !is.na(larvae_10m2) & larvae_10m2 > 0) %>% # Remove any blank columns (artifact of Excel)
      mutate(., .before = "STATION", CRUISE = sheet_names[i]) # Add in cruise code
    
    # Use character date/time formats -- I noticed these are different
    sheet_long$HOUR = as.character(sheet_long$HOUR)
    sheet_long$DATE <- as.Date(sheet_long$DATE, tryFormats = c("%d/%m/%Y", "%m/%d/%Y"))
    
    if (i == 1) {
      all_data = sheet_long #If this is the first time around the loop, make this all_data
    } else {
      all_data = bind_rows(all_data, sheet_long) # Otherwise, add the new sheet to all_data
    }
  }
  
  # Cast the integer value in species names to a character (necessary to merge and keep non-numeric character labels)
  species_codes$no = as.character(species_codes$no)
  
  # R creates new column names for duplicate columns (i.e., duplicate CalCOFI species codes) in the format "870...1".
  # So to return to the original CalCOFI number in order to link with species names, we need to remove the extra suffix
  all_data$species_code = vapply(strsplit(all_data$species_code, "...", fixed = TRUE), "[", "", 1)
  
  colnames(all_data) = tolower(colnames(all_data))
  
  # Link dataset to species names, convert to scientific names, add program, separate into line/station, standardize columns, extract month/day/year 
  all_data = merge(all_data, species_codes, by.x = "species_code", by.y = "no", all.x = TRUE) %>% 
    mutate(., .before = 1, program = "IMECOCAL") %>% 
    rename(., 
      latitude = lat,
      longitude = long,
      day_night = "hour (0=night; 1= day)",
      time = hour,
      surface_temp_oC = 'sup temp',
      surface_sal_psu = 'sup sal',
      tow_depth_m = 'tow depth (m)') %>%
    separate(., col = station, into = c("line", "station"), sep = c("[.]")) %>%   # the brackets are necessary because "." is a special character 
    mutate(., .after = "date", year = year(date), month = month(date), day = day(date), across(species_code, as.numeric), across(c("line", "station"), as.character)) %>% 
    mutate(., day_night = replace(day_night, day_night == 1, "D"), gear = "CB") %>% 
    mutate(., day_night = replace(day_night, day_night == 0, "N"))
  
  # Write file
  write.csv(x = all_data, file = "Imecocal.csv", row.names = FALSE)
  
  return(all_data)
}

# CalCOFI ----------------------------------------------------------------------
getcalcofi <- function() {
  
  calcofi <- getERDDAP(name = "erdCalCOFIlrvcntpos") %>% 
    as.data.frame(.) %>%
    mutate(., .before = 1, program = "CalCOFI") %>% 
    rename(volume_sampled_m3 = volume_sampled, gear = net_type) %>% 
    mutate(., across(c(volume_sampled_m3, latitude, longitude, larvae_10m2, larvae_1000m3), as.numeric) , across(c(cruise, line, station), as.character))
  
  hydros <- getERDDAP(name = "erdNOAAhydros") %>%
    subset(., standard_depth <= 210) %>%
    mutate(., across(c(cruise, line, station), as.character), across(c(standard_depth, temperature, salinity, density, oxygen, dynamic_height, percent_saturation, latitude, longitude), as.numeric)) %>%
    group_by_at(., c("cruise","line", "station", "date")) %>%
    summarise(., surface_sal_psu = mean(salinity),
              surface_temp_oC = mean(temperature),
              density_kg_m3 = mean(density),
              dissolved_oxygen_mL_L = mean(oxygen),
              dynamic_height_m = mean(dynamic_height))
  
  calcofi <- merge(
      x = calcofi,
      y = hydros,
      by.x = c("date", "line", "station", "cruise"),
      by.y = c("date", "line", "station", "cruise"),
      all.x = TRUE) 
    # select(., -cruise.y) %>%
    #rename(., species_code = "calcofi_species_code", cruise = cruise.x) %>% 
  
  write.csv(x = calcofi, file = "CalCOFI.csv", row.names = FALSE)
  
  #write.csv(calcofi_hydros, file = "CalCOFI_env.csv", row.names = FALSE)
  
  return(calcofi)
}

# RREAS ------------------------------------------------------------------------
getrreas <- function() {
  
  # Catch data
  rreas <- getERDDAP(name = "FED_Rockfish_Catch") %>% 
    mutate(., .before = 1, program = "RREAS") %>% 
    rename(., tow_number = haul_no, scientific_name = sci_name, larvae_count = catch) %>% 
    mutate(., across(c(latitude, longitude, bottom_depth), as.numeric) , across(c(station), as.character), day_night = "N") %>%   # all tows in RREAS occur at night
    mutate(., gear = "Cobb MWT") # all sampling is via a cobb midwater trawl
  
  # Environmental data
  rreas.env <- getERDDAP(name = "erdFedRockfishCtd") %>% 
    subset(.,depth <= 30) %>%
    mutate(., across(c(cruise, station), as.character),
           across(c(depth,
                    temperature,
                    salinity,
                    density,
                    transmissivity,
                    fluor_volt,
                    oxygen,
                    dyn_hgt,
                    chlorophyll,
                    latitude,
                    longitude),as.numeric)) %>%
    group_by_at(., c("cruise","station", "date")) %>%
    summarise(.,
              surface_sal_psu = mean(salinity),
              surface_temp_oC = mean(temperature),
              density_kg_m3 = mean(density),
              dissolved_oxygen_mL_L = mean(oxygen),
              dynamic_height_m = mean(dyn_hgt),
              irradiance_umol_m2_s = mean(irradiance),
              fluor_volt = mean(fluor_volt),
              oxygen_volt = mean(oxygen_volt),
              trans_percent = mean(transmissivity))
  
  rreas <- merge(rreas, rreas.env, all.x = TRUE)
  
  write.csv(x = rreas, file = "RREAS.csv", row.names = FALSE)
  
  return (rreas)
}

# PRERECRUIT -------------------------------------------------------------------
getPrerecruit <- function() {
  
  # LARVAE --------------
  # bongo data
  # Pivot original dataset to long format -> reformat dates -> add program column
  prs.bongo <- read_xlsx("OriginalDatasets/PreRecruit Bongo.xlsx", sheet = "Number Master") %>% 
    pivot_longer(data = ., cols = c(20:76), values_to = "larvae_count", names_to = "scientific_name") %>%
    mutate(., .before = 1, program = "PRS_larvae")
  # Standardize columns
  colnames(prs.bongo) = tolower(colnames(prs.bongo))
  prs.bongo <-
    rename(prs.bongo,
           station = 'station (new)',
           line = 'transect (new)',
           volume_sampled_m3 = 'volume filtered (m3)',
           tow_depth_m = 'haul depth (m)',
           bottom_depth = 'sta depth (m)',
           original_station = 'original station',
           original_transect = 'n/a transect') %>%
    mutate(., across(c("original_station", "line", "station"), as.character)) %>% 
    mutate(., .after = larvae_count, larvae_m3 = larvae_count/volume_sampled_m3)
  
  # JUVENILES --------------
  # Midwater trawl data
  prs.mwt = read_xlsx(path = "OriginalDatasets/PreRecruit MWT.xlsx", sheet = "Catch")
  prs.mwt.haul = read_xlsx(path = "OriginalDatasets/PreRecruit MWT.xlsx", sheet = "Haul")
  
  colnames(prs.mwt) = tolower(colnames(prs.mwt))
  colnames(prs.mwt.haul) = tolower(colnames(prs.mwt.haul))
  
  prs.mwt <- merge(prs.mwt, prs.mwt.haul[c("distance from shore (km)", "start latitude", "start longitude", "start time", "year", "month", "day", "station (new)", "transect (new)", "total revs")]) %>% 
    rename( .,
            latitude = 'start latitude', longitude = 'start longitude',
            line = 'transect (new)',
            time = 'start time',
            scientific_name = taxon,
            larvae_count = number,
            tow_depth_m = 'start depth (m)',
            station_notes = comments,
            original_station = 'original station',
            station = 'station (new)',
            tow_number = 'new haul #',
            surface_temp_oC = 'surface temp (oc)') %>%
    mutate(., .before = 1, program = "PRS_juveniles", across(c("line", "station"), as.character), gear = "Cobb MWT")
  
  # Combine both datasets
  prerecruit <- bind_rows(prs.bongo, prs.mwt) %>% mutate(., day_night = "N")
  
  # Add date column
  prerecruit <- mutate(prerecruit, .before = year, date = as.Date(paste(prerecruit$year, prerecruit$month, prerecruit$day, sep = "-"), "%Y-%m-%d"))
  
  # Import environmental data
  prerecruit.env <- read_xlsx("Environment/Prerecruit_CTD.xlsx", sheet =  "CTD Master")
  colnames(prerecruit.env) = tolower(colnames(prerecruit.env)) 
  
  # Add date column, rename columns
  prerecruit.env <- mutate(prerecruit.env, .before = year, date = as.Date(paste(prerecruit.env$year, prerecruit.env$month, prerecruit.env$day, sep = "-"), "%Y-%m-%d")) %>% 
    rename(
      line = 'transect (new)',
      station = 'station (new)',
      depth_m = 'depth (m)',
      surface_temp_oC = 'temperature (oc)',
      surface_sal_psu = 'salinity (psu)',
      density_kg_m3 = 'density (sigma-theta, kg/m3)',
      trans_percent = 'beam transmission (%)',
      fluor_volt = 'fluorescence (v)',
      fluor_mg_m3 = 'fluorescence (mg/m3)',
      dissolved_oxygen_mL_L = 'dissolved oxygen (ml/l)'
    ) %>% 
    subset(., depth_m < 100) %>% 
    group_by_at(., c("date", "line", "station")) %>% 
    summarize(surface_temp_oC = mean(surface_temp_oC))
  
  prerecruit <- merge(prerecruit, prerecruit.env, by = c("date", "line", "station"), all.x = TRUE)
  
  write.csv(prerecruit, file = "C://KDale/Projects/Phenology/Data/Prerecruit_full.csv")
  
  return(prerecruit)
}

# NH LINE ----------------------------------------------------------------------
getNHLine <- function() {
  
  nhline <- read_xlsx("OriginalDatasets/NHLine.xlsx", sheet = 1)

  colnames(nhline) = tolower(colnames(nhline))

  nhline <- rename(nhline, line = transect, volume_sampled_m3 = volume_m3_flowmeter, scientific_name = species, bottom_depth = 'station depth (m)', larvae_m3 = number_per_m3) %>% 
    mutate(., program = "NH-Line", .before = 1) %>% 
    mutate(., .after = "date", year = year(date), month = month(date), day = day(date)) %>% 
    mutate(., across(c("line", "station"), as.character)) %>% 
    mutate(larvae_count = larvae_m3 * volume_sampled_m3)
  
  # Environmental data
  nhline.env <- read.csv(file = "Environment/NHLine_CTD_lessthan100m.csv")
  colnames(nhline.env) <- tolower(colnames(nhline.env))
  nhline.env <- 
    rename(nhline.env,
           'station code' = station.code,
           cruise = cruiseid,
           line = transect,
           date = sample.date,
           latitude = lat,
           longitude = long,
           time = time_local) %>%
    mutate(., date = as.POSIXct(date, tryFormats = "%m/%d/%y"), across('station', as.character)) %>% 
    group_by_at(., c("station code")) %>% 
    summarize(., latitude = max(latitude),
              longitude = max(longitude),
              surface_sal_psu = mean(salinity),
              surface_temp_oC = mean(temperature),
              density_kg_m3 = mean(density),
              dissolved_oxygen_mL_L = mean(oxygen),
              fluor_volt = mean(fluor_volt),
              trans_percent = mean(tran_percent)) %>% 
    ungroup(.)
  
  nhline <- merge(nhline, nhline.env, by = c("station code"), all.x = TRUE)
  
  write.csv(nhline, file = "NHLine.csv", row.names = FALSE)
  
  return(nhline)
  
}

# CANADA -----------------------------------------------------------------------
getCanada <- function() {
  
  canada <- read_xlsx("OriginalDatasets/Canada.xlsx", sheet = "Fish")
  colnames(canada) = tolower(colnames(canada))
  canada <- rename(canada, time = stn_time, volume_sampled_m3 = 'volume filtered(m3)', latitude = lat, longitude = lon, scientific_name = name, bottom_depth = 'bottom depth(m)', day_night = twilight, gear = net_type, depth_strata = depth_strt1, larvae_m3 = 'abundance(#/m3)', phylum = 'phylum:', class = 'class:', order = 'order:', family = 'family:') %>% 
    mutate(., .after = "date", year = year(date), month = month(date), day = day(date)) %>% 
    mutate(., day_night = replace(day_night, day_night == "Daylight", "D") , day_night = replace(day_night, day_night == "Night", "N")) %>% 
    separate(., col = scientific_name, into = c("genus", "species", "maturity"), sep = c(" ", " ", " "), extra = "merge") %>% 
    mutate(., species = replace(species, species == "*sp.", NA)) %>% 
    unite(., col = scientific_name, genus, species, na.rm = TRUE, sep  = " ") %>% 
    mutate(., across(c(latitude, longitude, volume_sampled_m3, depth_strata, bottom_depth, depth_end1), as.numeric) , across(c(time), as.character)) %>% 
    mutate(., .before = larvae_m3, larvae_count = larvae_m3 * volume_sampled_m3) %>%
    mutate(., .before = 1, program = "Canada") 
  
  write.csv(x = canada, file = "Canada.csv", row.names = FALSE)
  
  return(canada)
}

# EcoFCOI ----------------------------------------------------------------------
getEcoFOCI <- function() {
  
  ecofoci <- read_xlsx("OriginalDatasets/EcoFOCI.xlsx", sheet = "DuplicatesRemoved")
  
  colnames(ecofoci) <- tolower(colnames(ecofoci))
  
  ecofoci <- separate(ecofoci, col = haul_id, sep = c(" "), into = c(NA, "station", NA, "gear", "net_num"), remove = TRUE) %>% 
    rename(., tow_number = haul_name, latitude = lat, longitude = lon, date = gmt_date_time, scientific_name = species_name, region_name = geographic_area, larvae_10m2 = larvalcatchper10m2, larvae_1000m3 = larvalcatchper1000m3, larvae_count = number_measured_counted) %>% 
    mutate(., date = as.Date(date), .before = 1, program = "EcoFOCI") %>%
    mutate(., .after = date, month = month(date), day = day(date))
  
  write.csv(ecofoci, "EcoFOCI.csv", row.names = FALSE)
  
  return(ecofoci)
}
# GET TOWS ---------------------------------------------------------------------
getTows <- function() {
  
  # calcofi_tows <- getERDDAP(name = "erdCalCOFItows") %>%
  #   mutate(across(c(latitude, longitude), as.numeric)) %>%
  #   mutate(., .before = 1, program = "CalCOFI", across(c(percent_sorted, latitude, longitude), as.numeric)) %>%
  #   rename(gear = net_type) %>%
  #   group_by_at(., c("program", "date", "line", "station", "gear")) %>%
  #   summarize(latitude = max(latitude), longitude = max(longitude), total_larvae = sum(total_larvae), larvae_10m2 = total_larvae * standard_haul_factor * percent_sorted) %>%
  #   ungroup(.) %>%
  #   mutate(.,  across(c("line", "station"), as.character))
 
  # rreas_tows <- group_by_at(rreas, c("program", "date", "station", "gear")) %>%
  #   subset(., !is.na(larvae_count)) %>% summarize(., latitude = max(latitude), longitude = max(longitude), total_larvae = sum(larvae_count))

  all_tows <- group_by_at(all_data_eez, c("towID", "program", "date", "latitude", "longitude", "year", "month", "day")) %>% 
    summarize(
      surface_temp_oC = max(surface_temp_oC),
      surface_sal_psu = max(surface_sal_psu),
      dissolved_oxygen_mL_L = max(dissolved_oxygen_mL_L),
      dynamic_height_m = max(dynamic_height_m),
      fluor_volt = max(fluor_volt),
      density_kg_m3 = max(density_kg_m3),
      trans_percent = max(trans_percent)
    ) %>% 
    mutate(day_of_year = yday(date))
  
  write.csv(x = all_tows, file = "AllTows.csv", row.names = FALSE)
  
  return(all_tows)
}

# SELECT WITHIN EEZ --------------------------------------------------------------
selectWithinEEZ <- function(data) {
  
  # Load EEZ shapefile
  eez <- read_sf("C://KDale/GIS/NorthAmerica_EEZ.shp")
  
  # Subset tows to within 200 km of land (EEZ)
  data_sf <- subset(data, !is.na(latitude)) %>% st_as_sf(., coords = c("longitude", "latitude"), remove = FALSE)
  data_sf$program = factor(data_sf$program, levels = c("CalCOFI","IMECOCAL","RREAS","PRS_juveniles","PRS_larvae", "NH-Line","Canada", "EcoFOCI"))
  st_crs(data_sf) <- st_crs(eez)
  eez_data <- st_intersection(data_sf, eez)
  
  return(eez_data)
  
}

# MAKE MAP ----------------------------------------------------------------------
createMap <- function(data) { 
  
  # Load North America shapefile
  northAmerica <- read_sf("C://KDale/GIS/NorthAmerica.shp")
  
  data = st_as_sf(data, coords = c("longitude", "latitude"))
  st_crs(data) = st_crs(northAmerica)
  
  jpeg("C://KDale/Projects/Phenology/Figures/All_samples.jpg",units = "in", width = 8, height = 8, res = 200)
  print(ggplot() +
          geom_sf(data = data, mapping = aes(color = program)) +
          geom_sf(data = northAmerica) +
          coord_sf(expand = FALSE) +
          scale_color_manual("Program", values = colors) +
          theme_bw(base_size = 14))
  dev.off()
}

# CREATE SPECIES TABLE ---------------------------------------------------------
createSpeciesTable <- function(data) {
  
  data = data %>% as.data.frame()
  
  speciesOfInterest <- filter(data, scientific_name %in% speciesNames) %>% subset(., larvae_count > 0)
  
  speciesAbundances <-
    group_by_at(speciesOfInterest, c("program","scientific_name")) %>%
    summarize(., n = sum(larvae_count)) %>% 
    pivot_wider(., id_cols = "scientific_name", names_from = "program", values_from = "n", names_prefix = "n_")
  speciesAbundances$nMin = min(speciesAbundances[,which(colnames(speciesAbundances)=="n_IMECOCAL"):which(colnames(speciesAbundances)=="n_EcoFOCI")])
  speciesAbundances$nMax = max(speciesAbundances[,which(colnames(speciesAbundances)=="n_IMECOCAL"):which(colnames(speciesAbundances)=="n_EcoFOCI")])
  
  years <- group_by_at(data, c("program", "year")) %>%
    summarise(n = n()) %>%
    group_by(., program) %>%
    summarise(program_years = n())
  
  speciesFrequency <-
    group_by_at(speciesOfInterest, c("program", "year", "scientific_name")) %>%
    summarise(n = n()) %>%
    group_by_at(., c("program","scientific_name")) %>%
    summarise(., seen_years = n()) %>% 
    merge(., years, by = "program") %>% 
    mutate(., freq_across_years = seen_years / program_years) %>% 
    pivot_wider(., id_cols = "scientific_name", names_from = "program", values_from = "freq_across_years", names_prefix = "freq_")
  speciesFrequency$freqMin = min(speciesFrequency[,which(colnames(speciesFrequency)=="freq_IMECOCAL"):which(colnames(speciesFrequency)=="freq_EcoFOCI")])
  speciesFrequency$freqMax = max(speciesFrequency[,which(colnames(speciesFrequency)=="freq_IMECOCAL"):which(colnames(speciesFrequency)=="freq_EcoFOCI")])
  
  speciesTable <- merge(speciesInfo,speciesAbundances,  by = "scientific_name") %>% 
    merge(., speciesFrequency, by = "scientific_name") %>%
    replace(is.na(.), 0)
  
  write.xlsx(speciesTable, file = "Species_table.xlsx")
  
  return(speciesTable)
}

# ROMS -------------------------------------------------------------------------
linkroms <- function(tows) {
  roms <- ncdf4::nc_open(filename = "ROMS/nep_srf_1995-2019.nc")
  
  # lon_rho, lat_rho: longitude and latitude of grid points
  # mask_rho: land-sea mask at grid points (0 = land, 1 = ocean)
  # h: total water depth at grid points
  # zeta: monthly sea surface height at grid points
  # temp: monthly sea surface temperature at grid points
  
  longitude=ncvar_get(roms, varid = "lon_rho")
  latitude=ncvar_get(roms, varid = "lat_rho")
  sst=ncvar_get(roms, varid = "temp")
  ssh=ncvar_get(roms, varid = "zeta")
  salinity = ncvar_get(roms, varid = "salt")
  
  time=data.frame(date = as.Date(ncvar_get(roms, varid = "ocean_time")/86400, origin = '1900-01-01')) %>% 
    mutate(., year = year(date), month = month(date)) %>% 
    mutate(., year_month = paste(year,month, sep = "-"))
  
  nc_close(roms) # close the netCDF file
  
  # create a "year-month" column (ROMS output indexed by month)
  tows$year_month = paste(tows$year, tows$month, sep = "-") 
  
  # progress bar
  pb <- txtProgressBar(min = 0, max = nrow(tows), char = "=", style = 3)
  tows$sst_roms = 0
  tows$ssh_roms = 0
  tows$salinity_roms = 0
  
  ### Match satellite data to sampling locations
  for (i in 1:nrow(tows)) {
    
    # Cancel current iteration if year is out of range
    if (tows$year[i] < min(time$year) | tows$year[i] > max(time$year)) 
      next
    
    targetLat = tows$latitude[i]
    
    # Get the best date match (ROMS output is monthly)
    dateIndex <-
      which(tows$year_month[i] == time$year_month)
    
    # Find all ROMS lat/lons within 0.1 deg of the sampling point
    nearbyLatitudeIndices <-
      which(abs(latitude - tows$latitude[i]) <= min(abs(latitude - tows$latitude[i])+0.1))
    nearbyLongitudeIndices <-
      which(abs(longitude - tows$longitude[i]) <= min(abs(longitude - tows$longitude[i])+0.1))
    
    # Find potential nearby points with valid latitude/longitudes
    matchIndices <-
      which(nearbyLatitudeIndices %in% nearbyLongitudeIndices == TRUE)
    
    # Calculate the best match (assumed to be the one with the closest latitude)
    differences = latitude[nearbyLatitudeIndices[matchIndices]] - targetLat
    bestMatch = matchIndices[which(abs(differences) == min(abs(differences)))]
    
    # Get the matching latitude
    matchLat = latitude[nearbyLatitudeIndices[bestMatch]]
    
    # Get row/column info (necessary for accessing the sst/ssh matrices)
    row = which(latitude == matchLat, arr.ind = TRUE)[1]
    column = which(latitude == matchLat, arr.ind = TRUE)[2]
    
    # Extract sst and ssh data
    tows$sst_roms[i] = sst[row, column, dateIndex]
    tows$ssh_roms[i] = ssh[row, column, dateIndex]
    tows$salinity_roms[i] = salinity[row, column, dateIndex]
    
    setTxtProgressBar(pb, i)
    
  }
  
  close(pb)
  
  tows[tows == 0] <- NA

  # Return
  return(tows)
  
}
# RUN FUNCTIONS ----------------------------------------------------------------
## Standardize datasets  -----
imecocal <- getIMECOCAL()
calcofi <- getcalcofi()
rreas <- getrreas()
prerecruit <- getPrerecruit()
nhline <- getNHLine()
canada <- getCanada()
ecofoci <- getEcoFOCI()

# Alternatively, import already-cleaned CSV versions (use mutate/across to assign correct data types)
imecocal <- read.csv("Imecocal.csv") %>% mutate(., across(date, as.Date)) %>% mutate(., across(c(cruise, line, station), as.character))
calcofi <- read.csv("CalCOFI.csv") %>% mutate(., across("date", as.Date)) %>% mutate(., across(c("cruise", "line", "station"), as.character))
rreas <- read.csv("RREAS.csv") %>% mutate(., across("date", as.Date)) %>% mutate(., across(c("cruise", "station"), as.character))
prerecruit <- read.csv("Prerecruit_full.csv") %>% mutate(., across("date", as.Date)) %>% mutate(., across(c("line", "station"), as.character))
nhline <- read.csv("NHLine.csv") %>% mutate(., across("date", as.Date)) %>% mutate(., across("station", as.character))
ecofoci <- read.csv("ecofoci.csv") %>% mutate(., across("date", as.Date)) %>% mutate(., across("gear", as.character))
canada <- read.csv("canada.csv") 
canada <- mutate(canada, across("date", as.Date))

## Combine datasets, select relevant columns -----
all_data <- bind_rows(list(imecocal, calcofi, rreas, prerecruit, nhline, canada, ecofoci)) %>%
  select(., program, cruise, line, station, date, year, month, day, time,
                   day_night, gear, latitude, longitude, scientific_name, larvae_count,
                   volume_sampled_m3, larvae_10m2, larvae_m3, larvae_1000m3, maturity,
                   tow_depth_m, bottom_depth,
                   surface_temp_oC, surface_sal_psu, dissolved_oxygen_mL_L, dynamic_height_m, fluor_volt, density_kg_m3, trans_percent)

## Select within 200nm EEZ -----
all_data_eez <- selectWithinEEZ(all_data) %>% st_drop_geometry(.) %>% 
  mutate(., larvae_count_scaled = scale(larvae_count)[,1], # scale catch data columns
         larvae_10m2_scaled = scale(larvae_10m2)[,1],
         larvae_m3_scaled = scale(larvae_m3)[,1],
         larvae_1000m3_scaled= scale(larvae_1000m3)[,1]) %>% 
  mutate(.,density_anomalies = coalesce( # Create one column of density anomalies
      larvae_10m2_scaled,
      larvae_m3_scaled,
      larvae_1000m3_scaled, 
      larvae_count_scaled))

## Add a unique tow ID (for cruises without "line" column, remove NAs)
all_data_eez <- mutate(all_data_eez, .after = program, towID = paste(all_data_eez$program, all_data_eez$date, all_data_eez$line, all_data_eez$station, all_data_eez$gear, sep = "_")) %>% 
  mutate(., .after = day, day_of_year = yday(date))
all_data_eez$towID = gsub(pattern = "_NA_", replacement = "_", x = all_data_eez$towID)

## Get tows and select within 200nm EEZ ----
## Add ROMS output to tows and add unique towID ----
all_tows <- getTows()
all_tows_eez <- 
  selectWithinEEZ(all_tows) %>%
  st_drop_geometry(.) %>%
  linkroms(.)

## Read/write ----
# Catch data
write.csv(all_data, file = "AllCruises_Combined.csv")
write.csv(all_data_eez, file = "AllCruises_Combined_200nm.csv", row.names = F)
all_data_eez <- read.csv(file = "AllCruises_Combined_200nm.csv")

# Tows
write.csv(all_tows, file = "AllTows.csv")
write.csv(all_tows_eez, file = "AllTows_200nm.csv", row.names=F) 
all_tows_eez <- read.csv(file = "AllTows_200nm.csv")

## Species table & map ----
speciesTable = createSpeciesTable(all_data_eez)
createMap(data = all_tows_eez)

# SPATIAL OPERATIONS ----------------------------------------------------------

# Create EEZ buffer for North America
# eez <- read_sf("C://KDale/GIS/World_EEZ_v11_20191118/World_EEZ_v11_20191118/eez_v11.shp") %>% 
#   st_make_valid(.) %>%  
#   st_crop(., xmin = -170, xmax = -109,
#           ymin = 20, ymax = 62,
#           expand=FALSE) %>% 
#   subset(., TERRITORY1 != "Hawaii" & TERRITORY1 != "Johnston Atoll") %>% 
#   st_union(.) %>% 
#   write_sf(., "C://KDale/GIS/NorthAmerica_EEZ.shp")

# northAmerica <- read_sf("C://KDale/GIS/North_South_America.shp") %>% 
#   st_crop(., xmin = -100, xmax = -180,
#           ymin = 15, ymax = 65,
#           expand=TRUE) %>% 
#   st_union(northAmerica)
# write_sf(northAmerica, "C://KDale/GIS/NorthAmerica.shp")

# SANDBOX --------------------------

