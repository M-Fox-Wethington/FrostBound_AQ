# Sea Ice Concentration Data Processing

This project automates the processing of Sea Ice Concentration (SIC) data from HDF5 files, converting it into georeferenced GeoTIFF format for spatial analysis. The primary focus is on batch processing of HDF5 files to extract SIC data, georeference it using the Antarctic Polar Stereographic South projection (EPSG:3412), and save the output as GeoTIFF files organized in date-based subdirectories.

## Features

- **Automated HDF5 Processing**: Batch process multiple HDF5 files containing SIC data.
- **Data Transformation**: Convert SIC data into a spatial format (sf) and rasterize it to GeoTIFF.
- **Structured Data Organization**: Automatically create date-based subdirectories for organized data storage.
- **Spatial Analysis Ready**: Outputs are georeferenced GeoTIFFs ready for spatial analysis.

## Getting Started

### Prerequisites

Ensure you have R installed on your system, along with the following R packages:

- `terra` for spatial data analysis and raster manipulation.
- `rhdf5` for reading HDF5 files.
- `sf` for handling spatial vector data.

You can install these packages using R commands like:

install.packages(c("terra", "Rhdf5lib", "sf"))

### Installation

Clone this repository to your local machine using:

bash

git clone https://github.com/M-Fox-Wethington/FrostBound_AQ.git

### Usage:

Place your HDF5 files containing Sea Ice Concentration data in the staged directory within D:/Manuscripts_localData/FrostBound_AQ/Datasets/AMSR-E_2/.

Run the ProcessSeaIceDataToGeoTIFF.R script in R:

source("Path/To/ProcessSeaIceDataToGeoTIFF.R")

Processed GeoTIFF files will be saved in the Processed directory, organized into subdirectories by year and month.


```r
