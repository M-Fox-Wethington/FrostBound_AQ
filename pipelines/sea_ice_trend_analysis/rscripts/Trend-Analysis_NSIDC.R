library(terra)

nsidc <- rast("D:/Manuscripts_localData/FrostBound_AQ/Datasets/NSIDC-Sea-Ice-Index/stack/NSIDC_25km_Study-Area.nc")
study_areas <- vect("D:/Manuscripts_localData/FrostBound_AQ/Datasets/gis-layers/study-area/Frostbound_AQ_Subregions_EPSG_3976.shp")

# Retrieve time attribute as dates
dates <- time(nsidc)
str(dates)

# Ensure the dates are in Date format (this is redundant but here for safety)
dates <- as.Date(dates)
str(dates)

# Filter for Antarctic winter months: June, July, August, September
winter_months <- which(format(dates, "%m") %in% c("06", "07", "08", "09")) #indexes
nsidc_winter <- nsidc[[winter_months]] #extract spatrasters layers by date indexs


for(i in unique(study_areas$Region)){
  
  region_maskl <- study_areas[[study_areas$Region == i,]]
  
  
}
