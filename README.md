# CINAR State of the Ecosystem - OA Data Synthesis

## Data Files
Full-resolution delayed-mode glider datasets containing raw pH voltages can be found on [RUCOOL's Glider ERDDAP Server](http://slocum-data.marine.rutgers.edu/erddap/index.html).

SBU01 delayed mode files (2022 and 2023) were downloaded from the [IOOS Glider DAC](https://gliders.ioos.us/erddap/index.html).

Fully-processed pH glider datasets that are adjusted for sensor time-lag can be found [here](https://marine.rutgers.edu/~lgarzio/cinar_soe/glider_data/).

Synthesized datasets that contain summer bottom-water pH and aragonite saturation state from all data sources can be found [here](https://marine.rutgers.edu/~lgarzio/cinar_soe/bottom_water_data/). Code used to generated the synthesized [vessel-based dataset](https://github.com/lgarzio/cinar-soe/blob/master/data_wrangler_summer_bottom_vessel.py) and [glider-based dataset](https://github.com/lgarzio/cinar-soe/blob/master/data_wrangler_summer_bottom_glider.py). Original vessel-based datasets were downloaded via the [CODAP documentation](https://essd.copernicus.org/articles/13/2777/2021/), and additional ECOMON datasets not included in the original CODAP dataset were downloaded via the [NCEI Ocean Carbon and Acidification Data Portal](https://www.ncei.noaa.gov/access/ocean-carbon-acidification-data-system-portal/).


## Citations
[CODAP-NA](https://essd.copernicus.org/articles/13/2777/2021/): Jiang, L.-Q., Feely, R. A., Wanninkhof, R., Greeley, D., Barbero, L., Alin, S., Carter, B. R., Pierrot, D., Featherstone, C., Hooper, J., Melrose, C., Monacci, N., Sharp, J. D., Shellito, S., Xu, Y.-Y., Kozyr, A., Byrne, R. H., Cai, W.-J., Cross, J., Johnson, G. C., Hales, B., Langdon, C., Mathis, J., Salisbury, J., and Townsend, D. W.: Coastal Ocean Data Analysis Product in North America (CODAP-NA) – an internally consistent data product for discrete inorganic carbon, oxygen, and nutrients on the North American ocean margins, Earth Syst. Sci. Data, 13, 2777–2799, https://doi.org/10.5194/essd-13-2777-2021, 2021.

Humphreys, M. P., Gregor, L., Pierrot, D., van Heuven, S. M. A. C., Lewis, E. R., and Wallace, D. W. R. (2020). [PyCO2SYS](https://pypi.org/project/PyCO2SYS/): marine carbonate system calculations in Python. Zenodo. doi:10.5281/zenodo.3744275.

Lewis, E. and Wallace, D. W. R. (1998) Program Developed for CO2 System Calculations, ORNL/CDIAC-105, Carbon Dioxide Inf. Anal. Cent., Oak Ridge Natl. Lab., Oak Ridge, Tenn., 38 pp., [https://salish-sea.pnnl.gov/media/ORNL-CDIAC-105.pdf](https://salish-sea.pnnl.gov/media/ORNL-CDIAC-105.pdf)