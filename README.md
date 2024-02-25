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
```

### Downloading Data

    Visit the NSIDC data page for the AMSR-E/AMSR2 dataset at NSIDC AU_SI12.
    Download the HDF5 (.he5) files from (https://nsidc.org/data/au_si12/versions/1)
    Save the downloaded HDF5 files to a directory on your local machine.


### Running the Script

After downloading the `ProcessSeaIceDataToGeoTIFF.R` script and placing your HDF5 files in the desired directory, follow these steps to process your data:

1. **Open RStudio or another R environment** of your choice. If you don't have RStudio, you can download it from [RStudio's official site](https://www.rstudio.com/products/rstudio/download/).

2. **Set your working directory** to the location where you saved the script. You can do this by using the `setwd()` function in R. For example:

    ```r
    setwd("C:/path/to/script")
    ```

    Replace `"C:/path/to/script"` with the actual path where you have saved the `ProcessSeaIceDataToGeoTIFF.R` script.

3. **Load the script** into your R session. You can do this by using the `source()` function, which will execute the script and load the `process_hdf5_file` function into your environment:

    ```r
    source("ProcessSeaIceDataToGeoTIFF.R")
    ```

4. **Prepare the list of HDF5 files** you want to process. The script expects a directory containing HDF5 (.he5) files. Ensure you have modified the `hdf5_dir` variable in the script to point to your directory of HDF5 files, or adjust the script to use a variable path.

5. **Call the `process_hdf5_file` function** for each HDF5 file you wish to process. If you have multiple files and want to automate this process, you can use the `list.files()` function to generate a list of file paths and then use `lapply()` to apply the processing function to each file, as shown in the script:

    ```r
    hdf5_files <- list.files(hdf5_dir, pattern = "\\.he5$", full.names = TRUE)
    lapply(hdf5_files, process_hdf5_file)
    ```

    This command lists all `.he5` files in the specified directory and then processes each file, converting it to a GeoTIFF and saving it in the designated output directory.

### Tips for Novice R Users

- **Ensure all required libraries are installed**: Before running the script, make sure you have installed all required R packages (`terra`, `rhdf5`, `sf`) using the `install.packages()` function.
- **Check your file paths**: File paths in R use forward slashes (`/`) or double backslashes (`\\`). Ensure your paths are correctly formatted to avoid errors.
- **Use RStudio for a better experience**: RStudio provides a user-friendly interface for R, making script execution and debugging more manageable.

By following these steps, even those new to R should be able to successfully run the script and process their sea ice concentration data.



Acknowledgments

    Thanks to the National Snow and Ice Data Center (NSIDC) for providing the valuable sea ice concentration data.
    Appreciation goes to the contributors and maintainers of the terra, rhdf5, and sf R packages.



