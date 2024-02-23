markdown

# Sea Ice Concentration Data Processing Project

This project provides an automated R script to process AMSR-E/AMSR2 Unified L3 Daily 12.5 km Brightness Temperatures, Sea Ice Concentration, Motion & Snow Depth Polar Grids, Version 1 data files. It converts the data from HDF5 format to georeferenced GeoTIFF format, suitable for spatial analysis. The processed data is then organized into date-based subdirectories for easy access and use.

## Data Source

The data is sourced from the National Snow and Ice Data Center (NSIDC) and pertains specifically to the AMSR-E/AMSR2 Unified L3 Daily 12.5 km Polar Grids, Version 1 dataset.

- **Dataset Information**: [AMSR-E/AMSR2 Unified L3 Daily 12.5 km Brightness Temperatures, Sea Ice Concentration, Motion & Snow Depth Polar Grids, Version 1](https://nsidc.org/data/au_si12/versions/1)

## Getting Started

### Prerequisites

Before running the script, ensure you have R installed on your system. You will also need the following R packages:

- `terra`: For spatial data analysis and raster manipulation.
- `rhdf5`: For reading HDF5 files.
- `sf`: For handling spatial vector data.

Install these packages using the following R command:

```r
install.packages(c("terra", "Rhdf5lib", "sf"))

Downloading Data

    Visit the NSIDC data page for the AMSR-E/AMSR2 dataset at NSIDC AU_SI12.
    Download the HDF5 (.he5) files for the dates and measurements you are interested in analyzing.
    Save the downloaded HDF5 files to a directory on your local machine, such as D:/Manuscripts_localData/FrostBound_AQ/Datasets/AMSR-E_2/staged/tmp.

Running the Script

    Clone this repository or download the ProcessSeaIceDataToGeoTIFF.R script to your local machine.

    bash

git clone https://github.com/M-Fox-Wethington/FrostBound_AQ.git

Modify the hdf5_dir variable in the script to point to the directory where you saved the downloaded HDF5 files.

Open an R session and set the working directory to where the script is located, or use setwd() to change the working directory to the script's location.

Run the script in your R session:

r

    source("Path/To/ProcessSeaIceDataToGeoTIFF.R")

The script will automatically process each HDF5 file, creating GeoTIFF outputs organized into year and month subdirectories within a specified Processed directory.
For Novice R Users

If you are new to R, here are a few tips to get started:

    R Installation: Download and install R from CRAN, the Comprehensive R Archive Network.
    RStudio: Consider using RStudio, an integrated development environment (IDE) for R, to make script editing and execution easier.
    Learning R: Numerous free resources are available online to learn R, including R for Data Science and the Introduction to R manual.

Contributing

We welcome contributions to improve this project. Please follow the standard GitHub flow: fork the repo, make your changes, and submit a pull request.
License

This project is licensed under the MIT License - see the LICENSE file for details.
Acknowledgments

    Thanks to the National Snow and Ice Data Center (NSIDC) for providing the valuable sea ice concentration data.
    Appreciation goes to the contributors and maintainers of the terra, rhdf5, and sf R packages.

vbnet


Please copy the entire block above into your GitHub project's README.md file. This mark
