Contact: muhammad.farhan.ul.moazzam@uqtr.ca 
By: Muhammad Farhan Ul Moazzam
Date: 02/26/2026

Surface Urban Heat Island Dynamics Across Australian Cities (2003–2025)

Overview

This repository contains the complete workflow used to analyze Surface Urban Heat Island (SUHI) intensity and urban thermal stress across Adelaide, Brisbane, Melbourne, Perth, and Sydney from 2003 to 2025. The project integrates MODIS land surface temperature, land cover, and vegetation products with ERA5 reanalysis data within Google Earth Engine, followed by statistical analysis and visualization in Python. The structure is designed to ensure transparency, reproducibility, and straightforward extension to other regions.

Monthly SUHI and Seasonal Comparisons

Monthly SUHI is derived from the contrast between urban and rural land surface temperatures within each city boundary, alongside associated urban–rural vegetation differences. Seasonal summaries are generated for Annual (January–December), Summer (DJF), and Winter (JJA) conditions. To evaluate temporal evolution, the study period is divided into an Early phase (2003–2013) and a Late phase (2014–2025), and their differences are computed for both daytime and nighttime SUHI. The workflow exports city-level time series, seasonal summary tables, and multi-band GeoTIFF stacks representing all season–period–diurnal combinations.

Spatial SUHI and UTFVI Mapping

Beyond city-mean metrics, the repository produces pixel-level SUHI maps to capture intra-urban spatial variability. Each urban pixel is evaluated relative to the city’s rural baseline temperature, and seasonal averages are exported as structured raster stacks. In parallel, Urban Thermal Field Variance Index (UTFVI) maps are generated to quantify normalized thermal stress within the urban footprint. These products allow assessment of both relative urban–rural contrast and internal thermal heterogeneity.

Meteorological Drivers (ERA5)

Monthly ERA5 meteorological variables are extracted and spatially averaged over each city to provide climatic context. The dataset includes 2 m air temperature, relative humidity, wind speed, total precipitation, cloud cover fraction, and boundary layer height. These variables are temporally aligned with the SUHI series and support regression and attribution analyses linking urban heat intensity to background atmospheric conditions.

Repository Structure and Reproducibility

The repository includes Google Earth Engine scripts for data extraction and export, along with Jupyter notebooks for downstream analysis, persistence assessment, robustness testing, and figure generation. All masking strategies, thresholds, temporal definitions, and aggregation steps are implemented directly in code to ensure clarity and reproducibility. The framework is modular and scalable, enabling application to additional cities or comparative regional studies with minimal adjustment.

