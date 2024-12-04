# CINAR State of the Ecosystem - OA Data Synthesis

## Data Files

### 2025 Analysis

Full-resolution delayed-mode glider datasets can be downloaded from the [IOOS Glider DAC](https://gliders.ioos.us/erddap/index.html). A table of all pH glider datasets used in this analysis can be found [here](https://github.com/lgarzio/cinar-soe/blob/master/analysis2025/glider_deployment_summary_2019_2024.csv).

Vessel-based data were mined from the [Coastal Ocean Data Analysis Product in North America](https://essd.copernicus.org/articles/13/2777/2021/), version v2021 (Jiang et al. 2021). More recent vessel-based datasets that were not included in CODAP-NA v2021 were downloaded via the [NCEI Ocean Carbon and Acidification Data Portal](https://www.ncei.noaa.gov/access/ocean-carbon-acidification-data-system-portal/).

Synthesized datasets that contain bottom- and surface-water pH and aragonite saturation state from glider- and vessel-based data sources can be found [here](https://marine.rutgers.edu/~lgarzio/cinar_soe/2025_submission/data_files/). Code used to generated the synthesized vessel-based [bottom](https://github.com/lgarzio/cinar-soe/blob/master/analysis2025/data_wrangler_bottom_vessel.py) and [surface](https://github.com/lgarzio/cinar-soe/blob/master/analysis2025/data_wrangler_surface_vessel.py) datasets and [glider-based datasets](https://github.com/lgarzio/cinar-soe/blob/master/analysis2025/data_wrangler_bottom_surface_glider.py).


Satellite-derived chlorophyll and sea surface temperature data files were downloaded from the NASA Ocean Biology Distributed Active Archive Center (OB.DAAC) [Level 3 & 4 Browser](https://oceancolor.gsfc.nasa.gov/l3/) with the following selections:
- Product Status: Standard
- Instrument: SNPP-VIIRS
- Product: Chlorophyll concentration, Triple window sea surface temperature
- Period: Monthly
- Resolution: 4km
- Start Date: 2012-03-01
- End Date: 2024-08-31

Monthly data files of bottom temperature from the E.U. Copernicus Marine Service Information (CMEMS): [Global Ocean Physics Reanalysis (GLORYS12V1)](https://data.marine.copernicus.eu/product/GLOBAL_MULTIYEAR_PHY_001_030/description) dataset were downloaded using the Python Copernicus Marine Toolbox API for the following time ranges:
- Dataset ID: cmems\_mod\_glo\_phy\_my\_0.083deg_P1M-m: March 2012 to June 2021
- Dataset ID: cmems\_mod\_glo\_phy\_myint\_0.083deg_P1M-m: July 2021 to August 2024

## Citations
[CODAP-NA](https://essd.copernicus.org/articles/13/2777/2021/): Jiang, L.-Q., Feely, R. A., Wanninkhof, R., Greeley, D., Barbero, L., Alin, S., Carter, B. R., Pierrot, D., Featherstone, C., Hooper, J., Melrose, C., Monacci, N., Sharp, J. D., Shellito, S., Xu, Y.-Y., Kozyr, A., Byrne, R. H., Cai, W.-J., Cross, J., Johnson, G. C., Hales, B., Langdon, C., Mathis, J., Salisbury, J., and Townsend, D. W.: Coastal Ocean Data Analysis Product in North America (CODAP-NA) – an internally consistent data product for discrete inorganic carbon, oxygen, and nutrients on the North American ocean margins, Earth Syst. Sci. Data, 13, 2777–2799, https://doi.org/10.5194/essd-13-2777-2021, 2021.

Humphreys, M. P., Gregor, L., Pierrot, D., van Heuven, S. M. A. C., Lewis, E. R., and Wallace, D. W. R. (2020). [PyCO2SYS](https://pypi.org/project/PyCO2SYS/): marine carbonate system calculations in Python. Zenodo. doi:10.5281/zenodo.3744275.

Lewis, E. and Wallace, D. W. R. (1998) Program Developed for CO2 System Calculations, ORNL/CDIAC-105, Carbon Dioxide Inf. Anal. Cent., Oak Ridge Natl. Lab., Oak Ridge, Tenn., 38 pp., [https://salish-sea.pnnl.gov/media/ORNL-CDIAC-105.pdf](https://salish-sea.pnnl.gov/media/ORNL-CDIAC-105.pdf)


### 2024 Analysis

Full-resolution delayed-mode glider datasets containing raw pH voltages can be found on [RUCOOL's Glider ERDDAP Server](http://slocum-data.marine.rutgers.edu/erddap/index.html).

SBU01 delayed mode files (2022 and 2023) were downloaded from the [IOOS Glider DAC](https://gliders.ioos.us/erddap/index.html).

Fully-processed pH glider datasets that are adjusted for sensor time-lag can be found [here](https://marine.rutgers.edu/~lgarzio/cinar_soe/2024_submission/glider_data/files_for_bottom_water_synthesis/).

Synthesized datasets that contain summer bottom-water pH and aragonite saturation state from all data sources can be found [here](https://marine.rutgers.edu/~lgarzio/cinar_soe/2024_submission/bottom_water_data/). Code used to generated the synthesized [vessel-based dataset](https://github.com/lgarzio/cinar-soe/blob/master/analysis2024/data_wrangler_summer_bottom_vessel.py) and [glider-based dataset](https://github.com/lgarzio/cinar-soe/blob/master/analysis2024/data_wrangler_summer_bottom_glider.py). Original vessel-based datasets were downloaded via the [CODAP documentation](https://essd.copernicus.org/articles/13/2777/2021/), and additional EcoMon datasets not included in the original CODAP dataset were downloaded via the [NCEI Ocean Carbon and Acidification Data Portal](https://www.ncei.noaa.gov/access/ocean-carbon-acidification-data-system-portal/).