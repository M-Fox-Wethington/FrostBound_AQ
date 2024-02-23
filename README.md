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
Downloading Data

    Visit the NSIDC data page for the AMSR-E/AMSR2 dataset.
    Download the HDF5 files for the dates and measurements you are interested in analyzing.
    Save the downloaded HDF5 files to a directory on your local machine.

Running the Script

    Clone this repository or download the ProcessSeaIceDataToGeoTIFF.R script to your local machine.

    bash

git clone https://github.com/M-Fox-Wethington/FrostBound_AQ.git

Modify the hdf5_dir variable in the script to point to the directory where you saved the downloaded HDF5 files.

Open an R session and set the working directory to where the script is located, or use setwd() to change the working directory.

Run the script in your R session:

r

    source("Path/To/ProcessSeaIceDataToGeoTIFF.R")

The script will automatically process each HDF5 file, creating GeoTIFF outputs organized into year and month subdirectories within a specified Processed directory.
For Novice R Users

If you are new to R, here are a few tips to get started:

    R Installation: Download and install R from CRAN, the Comprehensive R Archive Network.
    RStudio: Consider using RStudio, an integrated development environment (IDE) for R, to make script editing and execution easier.
    Learning R: Numerous free resources are available online to learn R, including R for Data Science and the Introduction to R manual.


```r
