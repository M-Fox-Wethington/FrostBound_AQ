# Load required libraries
library(terra)
library(rhdf5)

# Set working directory and specify HDF5 file path
wd <- "D:/Manuscripts_localData/FrostBound_AQ/Datasets" 
h5_file <- paste0(wd,"/178568954/AMSR_U2_L3_SeaIce12km_B04_20120704.he5")

# View subdatasets within the HDF5 file
he_file <- "D:/Manuscripts_localData/FrostBound_AQ/Datasets/178568954/AMSR_U2_L3_SeaIce12km_B04_20120704.he5"
View(h5ls(he_file, all = TRUE))

# Read latitude and longitude data
h5grid <- h5read(h5_file, "/HDFEOS/GRIDS/SpPolarGrid12km")
lat <- h5grid$lat
str(lat)
lon <- h5grid$lon
str(lon)

# Read SeaIce data
SeaIce <- h5read(h5_file, "/HDFEOS/GRIDS/SpPolarGrid12km/Data Fields/SI_12km_SH_ICECON_DAY")
str(SeaIce)

# Read XDim and YDim data
XDim <- h5read(h5_file, "/HDFEOS/GRIDS/SpPolarGrid12km/XDim")
XDim
YDim <- h5read(h5_file, "/HDFEOS/GRIDS/SpPolarGrid12km/YDim")
YDim

# Read detailed metadata about the dataset
Core_Metadata <- h5read(h5_file, "/HDFEOS INFORMATION/CoreMetadata.0")
Core_Metadata


Struct_Metadata <- h5read(h5_file, "/HDFEOS INFORMATION/StructMetadata.0")
Struct_Metadata
str(Struct_Metadata)
