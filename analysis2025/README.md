# CINAR State of the Ecosystem - OA Data Synthesis

## Glider- and Vessel-based datasets and processing
Full-resolution delayed-mode glider datasets can be downloaded from the [IOOS Glider DAC](https://gliders.ioos.us/erddap/index.html). A table of all pH glider datasets used in this analysis can be found [here](https://github.com/lgarzio/cinar-soe/blob/master/analysis2025/glider_deployment_summary_2019_2024.csv).

Vessel-based data were mined from the [Coastal Ocean Data Analysis Product in North America](https://essd.copernicus.org/articles/13/2777/2021/), version v2021 (Jiang et al. 2021). More recent vessel-based datasets that were not included in CODAP-NA v2021 were ollected during more recent NOAA NEFSC Ecosystem Monitoring (EcoMon) surveys: June 2019 (Cruise ID HB1902), August 2019 (Cruise ID GU1902), October 2019 (Cruise ID GU1905), May 2021 (Cruise ID GU2102), August 2021 (Cruise ID PC2104), October 2021 (Cruise ID PC2106), June 2022 (Cruise ID HB2204), November 2022 (Cruise ID PC2205), June 2023 (Cruise ID HB2302), and the East Coast Ocean Acidification ECOA-3 Cruise (Cruise ID ECOA3). These additional datasets were downloaded via the [NCEI Ocean Carbon and Acidification Data Portal](https://www.ncei.noaa.gov/access/ocean-carbon-acidification-data-system-portal/).

For MAB glider datasets, total alkalinity was calculated from salinity using a linear relationship determined from in situ water sampling data taken during glider deployment and recovery in addition to ship-based water samples (Wright-Fairbanks et al. 2020). For the Gulf of Maine glider datasets, total alkalinity (provided in the files) was calculated from salinity using the linear relationship (TA = 47.6 * salinity + 643.0) taken from Hunt et al 2021, figure 6, Historical regression for Gulf of Maine: https://doi.org/10.1016/j.marchem.2021.103960. Calculations for aragonite were then conducted using PyCO2SYS (Humphreys et al. 2020) with inputs of pressure, temperature, salinity, total alkalinity, and pH.

For vessel-based datasets, when aragonite was unavailable it was calculated using PyCO2SYS (Humphreys et al. 2020) with inputs of pressure, temperature, salinity, total alkalinity, and pH.

Bottom water pH and aragonite values were defined as the median of the measurements (or calculated aragonite values) within the deepest 1m of a glider profile or, for vessel-based measurements, the deepest measurement of a vertical CTD/Rosette cast where water samples were collected, for profiles deeper than 10m. In order to validate whether the deepest depth was at or near the bottom, the sampling depth was compared to water column depth (when provided) or water depths extracted from a GEBCO bathymetry grid based on the sample collection coordinates. Any glider profiles/vessel-based casts with the deepest measurement shallower than the bottom 20% of total water column depth were removed. This allowed for a sliding scale instead of providing a strict cut off (e.g., 1 m above the bottom). 

Surface water pH and aragonite values were defined as the median of the measurements (or calculated ΩArag values) recorded at the top of a glider profile (between 2-4m depth) for profiles deeper than 10m or, for vessel-based measurements the shallowest measurement of a vertical CTD/Rosette cast where water samples were collected provided the measurement depth was <10m, or all of the surface measurements provided for underway datasets.

Synthesized datasets that contain bottom- and surface-water pH and aragonite saturation state from glider- and vessel-based data sources can be found [here](https://marine.rutgers.edu/~lgarzio/cinar_soe/2025_submission/data_files/). Code used to generated the synthesized vessel-based [bottom](https://github.com/lgarzio/cinar-soe/blob/master/analysis2025/data_wrangler_bottom_vessel.py) and [surface](https://github.com/lgarzio/cinar-soe/blob/master/analysis2025/data_wrangler_surface_vessel.py) datasets and [glider-based datasets](https://github.com/lgarzio/cinar-soe/blob/master/analysis2025/data_wrangler_bottom_surface_glider.py).

## Satellite-derived data
Satellite-derived chlorophyll and sea surface temperature data files were downloaded from the NASA Ocean Biology Distributed Active Archive Center (OB.DAAC) [Level 3 & 4 Browser](https://oceancolor.gsfc.nasa.gov/l3/) with the following selections:
- Product Status: Standard
- Instrument: SNPP-VIIRS
- Product: Chlorophyll concentration, Triple window sea surface temperature
- Period: Monthly
- Resolution: 4km
- Start Date: 2012-03-01
- End Date: 2024-08-31

## Bottom temperature data
Monthly data files of bottom temperature from the E.U. Copernicus Marine Service Information (CMEMS): [Global Ocean Physics Reanalysis (GLORYS12V1)](https://data.marine.copernicus.eu/product/GLOBAL_MULTIYEAR_PHY_001_030/description) dataset were downloaded using the Python Copernicus Marine Toolbox API for the following time ranges:
- Dataset ID: cmems\_mod\_glo\_phy\_my\_0.083deg_P1M-m: March 2012 to June 2021
- Dataset ID: cmems\_mod\_glo\_phy\_myint\_0.083deg_P1M-m: July 2021 to August 2024

### Citations
[CODAP-NA](https://essd.copernicus.org/articles/13/2777/2021/): Jiang, L.-Q., Feely, R. A., Wanninkhof, R., Greeley, D., Barbero, L., Alin, S., Carter, B. R., Pierrot, D., Featherstone, C., Hooper, J., Melrose, C., Monacci, N., Sharp, J. D., Shellito, S., Xu, Y.-Y., Kozyr, A., Byrne, R. H., Cai, W.-J., Cross, J., Johnson, G. C., Hales, B., Langdon, C., Mathis, J., Salisbury, J., and Townsend, D. W.: Coastal Ocean Data Analysis Product in North America (CODAP-NA) – an internally consistent data product for discrete inorganic carbon, oxygen, and nutrients on the North American ocean margins, Earth Syst. Sci. Data, 13, 2777–2799, https://doi.org/10.5194/essd-13-2777-2021, 2021.

Humphreys, M. P., Gregor, L., Pierrot, D., van Heuven, S. M. A. C., Lewis, E. R., and Wallace, D. W. R. (2020). [PyCO2SYS](https://pypi.org/project/PyCO2SYS/): marine carbonate system calculations in Python. Zenodo. doi:10.5281/zenodo.3744275.

Lewis, E. and Wallace, D. W. R. (1998) Program Developed for CO2 System Calculations, ORNL/CDIAC-105, Carbon Dioxide Inf. Anal. Cent., Oak Ridge Natl. Lab., Oak Ridge, Tenn., 38 pp., [https://salish-sea.pnnl.gov/media/ORNL-CDIAC-105.pdf](https://salish-sea.pnnl.gov/media/ORNL-CDIAC-105.pdf)
