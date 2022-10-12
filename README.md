# arcade
_Automated Rapid Climatological Azimuth Deviation Estimation_

## Purpose
This project aims to understand the feasibility of a simplified rapid azimuth deviation prediction scheme for improving the source location when implementing an association and location cross-bearings, grid-search method like IMS\_vASC [Matoza et al., 2017].
This implies using "offline", portable, and available for research purposes tools, so the calculations could be performed by any scientist in their personal workstation machine, but could also scale to a large research computation center.

## Project code
The project itself is a collection of Python scripts plus some Fortran 90 layers that communicate with 
available open source climatologies and a 3D ray-tracing algorithm. 
`ARCADE_main.py` contains the functions that serve to find and estimate the azimuth deviations, representing the core of the used methodology. 

## Current capabilities

ARCADE is cabable of different types of ray propagation modelling. The most 
basic model is range independent (spherically layered atmosphere), meaning that 
uses average atmospheric descriptions in height for a source to receiver (station)
direction slice ('profile'). This case uses only climatologies to determine the
properties of the atmosphere for a given datetime and profile. This case can be
tuned to use "perturbed" climatological profiles.

The second range intependent modeling type is a hybrid calculation that merges
ECMWF ERA 5 descriptions for the lower ~80 km, with climatologies for the 
upper atmosphere. This case needs internet connection, and could fail if
the download process hangs. 

The third type of modeling is using range dependent climatological atmospheric
descriptions. This case is more realistic than the range independent, as
it estimates the 3D atmosphere by interpolating nodal columns in height that
describe the area where the profile is. Note that this case is not 'rapid' anymore,
so it should be used as a better, more accurate estimate of the atmospheric
effect on the propagating signal. In theory it will take at least 10 times 
more time to calculate each one of this cases. Pertubed profileas have not 
been implemented in this case.

The fourth type of modeling is a range dependent hybrid set of descriptions as
in the range independent hybrid case, but for nodal columns that merge ECMWF
ERA 5 descriptions from 0 to ~80 km, with climatologies for >80 km heights. 
This case also depends on an internet connection, so it's the most realistic
case, but could take much more time (or even hang). Use it when you want to 
narrow down some effects, or you want more accurate values after identiying
interesting cases.


## Requisites
- All external repos (see below) in place and hopefully individually tested (thorugh makefiles inside each one). 
- Fortran and c++ compilers (check `makefile` to define Fortran compiler in first line)
  - For NRMLSIS2.0 you'll need gfortran version 4.8.5 or 6.3.0 or 8.1.0
  - For HWM14 you'll need fftw3 (libfftw3-dev in Ubuntu; fftw-3 with port in MacOS).
- Python 3.9 Conda: follow instructions here https://docs.conda.io/en/latest/miniconda.html
  - Create the `arcade` envinroment with `$ conda env create -f environment.yml`    
  - When running the calculations, activate the environment with `$ conda activate <env>`. **NOTE:** If using hybrid models, you need to install the dependencies here too (see below).

### External repos

The climatologies used are the Horizontal Wind Model (HWM14) [Drob et al., 2014] 
to estimate horizonal winds in height, and NRLMSIS2.0 [Emmert et al., 2021] to 
estimate other properties of the atmosphere in height, like density, temperature 
and pressure. The 3D ray-tracing algorithm used is infraGA [Blom and Waxler, 2012].

The method has two main external components: the atmospheric model (HWM14+NRLMSIS2.0) and the 3D ray tracing algorithm (infraGA/geoAC). Each one of these components has an associated repository that needs to be inside `repos` folder, with the names as defined in the `makefile`. That is, `./repos/HWM14`, `./repos/NRLMSIS2.0` and `./repos/infraGA-master`. Below are the instructions to access each one of the components.

#### NRLMSIS2.0
NRLMSIS 2.0 Code and all data samples used in this work are available at https://map.nrl.navy.mil/map/pub/nrl/NRLMSIS/NRLMSIS2.0
The associated publication is __Emmert, J. T., Drob, D. P., Picone, J. M., Siskind, D. E., Jones, M., Mlynczak, M. G., et al. (2021). NRLMSIS 2.0: A whole-atmosphere empirical model of temperature and neutral species densities. Earth and Space Science, 8, e2020EA001321. https://doi.org/10.1029/2020EA001321__


#### HWM14
Available to download from the supporting information of __Drob, D. P., Emmert, J. T., Meriwether, J. W., Makela, J. J., Doornbos, E., Conde, M., Hernandez, G., Noto, J., Zawdie, K. A., McDonald, S. E., et al. (2015), An update to the Horizontal Wind Model (HWM): The quiet time thermosphere, Earth and Space Science, 2, 301– 319, doi:10.1002/2014EA000089.__

#### infraGA/geoAC
Available from the public seosmoacoustic repository of Los Alamos National Laboratory https://github.com/LANL-Seismoacoustics/infraGA

### Hybrid models
In addition to the empirical climatologies (NRMSLSIS2.0/HWM14), you can also use hybrid atmospheric descriptions that are a combination of ERA 5 ECMWF below ~80 km and empirical climatologies at higher altitudes. These models were used to compare and validate the climatological estimations, but can still be used for realistic aproaches.

#### CDSAPI Python module
The hybrid models need the Climate Data Store (CDS) infrastructure and API module, `cdsapi`, provided by the European Centre for Medium-Range Weather Forecasts (ECMWF, https://www.ecmwf.int/). Please visit https://cds.climate.copernicus.eu/api-how-to to ensure the module is installed in your system, and properly configured API user and key to be able to download the atmospheric descriptions.

#### ecCodes to decode data
In order to decode the downloaded data and use it automatically for ray tracing with infraGA, it is necessary to install `ecCodes` provided by ECMWF. Please follow the instructions here https://confluence.ecmwf.int/display/ECC/ecCodes+installation. I recommend to use the Python binding installation (https://confluence.ecmwf.int/display/UDOC/How+to+install+ecCodes+with+Python+bindings+in+conda+-+ecCodes+FAQ).

## How to run

1. Inside repo, run `$ make all-prep`. 
2. Modify the input files `arcade_config.toml` and `discretize_parameters.toml`
to suit your model. (Currently it is setup with an example run for Puyehue-Cordon Caulle to IS02 on 2011-06-04 at 19:00:00 UTC)
3. Run `$ make run-arcade` to calculate the azimuth deviations.
4. After the modeling ends, the input configuration files and output figures
and results will be saved in the directory `path`/`name`, which is set in
`arcade_config.toml`. Inside `./output`, the text file `azimuth_deviation_table.txt` has the summary of results.

## Parameters and options

### `./input/arcade_config.toml`
Here you set up the ARCADE_main.py run parameters.
- `run_type`:
  `norm` - will use default non-perturbed atmospheric descriptions
  `pert` - will use perturbed atmospheric descriptions
- `use_thermo` : `true` or `false`, determines if using thermospheric arrivals or not
- `max_run` : maximum number of tries before giving up in searching for arrivals near each target station
- `atten_th` : attenuation threshold filter. Negative number representing the minimum attenuaton an arrival can have before not being considered in the estimations.
- `min_dist_arrv` : distance in km that sets the maximum distance at which an arrival can be to be associated with a station in the first search process (later it gets reduces to 5 km)
- `daz` : maximum distance in degrees to consider an arrival associated with a station. Measures azimuth difference.
- `dphi` : determines the aperture to launch rays around an azimuth value
- `scale` : this parameter scales the azimuth steps (the bigger the smaller the step). Check `ARCADE_main.py` to find relationship.
- `thresh`: this is the final distance in degrees to declare the search as completed.
- `k` : parameter involved in iterative search to look for arrivals around station before reaching `thresh`.

Below `[launch_parameters]` some other parameters related with the specific launch are set. These are related with the 3D ray tracing process:
- `incl_min` : minimum vertical launch angle in degrees
- `incl_max` : maximum vertical launch angle in degrees
- `incl_step` : step in degrees for launch angles, determine also the number of rays that are launched
- `drng` : maximum extra range given a source-station geographical distance to model
- `bounces` : maximum number of 'bounces' to model. A bounce is then the ray reaches a turning height and comes back to the ground.
- `src_alt` : altitude of source in kilometers
- `write_atmo` : `false` or `true`, depending if you want to write the atmospheric profiles (it slows down calculations so starts as `false`).
- `calc_amp` : `false` or `true`, if calculating the amplitudes for the rays. As it slows down calculations so it's should be `false` for testing.
- `perc_cpu` : from all available cpu/threads, how many will you use? Calculated as number of cpus/perc_cpu, so 1=100% of total, 2=50% of total, 3=33.333% of total and so on.

Below `[project_info]` a couple of parameters about the project name and location are set. These are important to save the results and running parameters after the modeling has been completed (check `wrap_proj.sh`).
- `name` : name of project as string between quotes. Example "Test_0".
- `path` : path of project as string. Example: `/home/yourname/Desktop/`. The project will be save inside this path after running with the name set in `name`. That is: `/home/yourname/Desktop/Test_0".


### `./input/discretize_parameters.toml`
Here you set up the parameters used in the discretization process. That is, when the area is defined thorugh the station locations, souce locations, grid steps, etc.
- `year` : year as YYYY. Example `2011`.
- `doys` : list of days of year to be calculated. Example `[123, 124, 125, 145]`. It can also be just one day like `[155]`.
- `sec` : time of day in seconds. Example: `68400` (19:00).
- `hmin` : minimum altitude considered in kilometers for atmospheric descriptions
- `hmax` : maximum altitude in kilometers
- `dh` : altitude step to calculate atmospheric descriptions in kilometers
- `f107a` : determines the solar activity for the climatologies, set as 'moderate' by default (value `150`).
- `f107` : related with `f107a`, set as same value.
- `apd` : determines the geomagnetic activity for the climatologies, set as 'quiet' by default (value `4`).
- `aph` : related with `apd`, set as same number.
- `sou_pos` : list of source locations in latitude-longitude pairs. Example: `[[-41.33, -72.62], [-40.59, -72.117]]`. It can be just one source as `[[-40.59, -72.117]]`.
- `sta_pos` : list of station locations. Analogous to `sou_pos`. 
- `sta_nam` : list of names of stations. Example: `["IS02", "IS08", "IS09"]`. It should correspond to `sta_pos`.
- `sta_skp` : Not used in this version. List of station numbers from list you want to skip in the plotting process.
- `ds` : step in kilometers to discretize along each path for calculating the average per height

Below `[range_dependent]` there are parameters to determine is using a range-dependent calculation or not plus the grid discretization.
- `use_rng_dep` : `false` or `true`. By default is false as the range-repdendent calculations where not used in the paper associated with this repo. It is also much slower than the range-independent cases.
- `dlat` and `dlon`: the step size in degrees for the grids.

Below `[ecmwf]` there are parameters related with the use of ERA 5 ECMWF data to construct the hybrid atmospheric descriptions.
- `use_ecmwf`: `false` or `true`, depending if using hybrid model. Of coure, `true` means using it, while `false` not using it.
- `auto_area`: `true` or `false`. Set as `false` if setting the area to be downloaded manually or `true` if determined by the source and station locations.
- `min_lon`, `max_lon`, `min_lat`, `max_lat`: sets the area if `auto_area` is `false`.
- `dlon`, `dlat`: sets the grid step of the data. It depends on ERA 5 so I recomment 1 degree.
- `h1` : sets minimum altitude to merge with climatologies. 
- `h2`: sets maximum altitude to merge with climatologies. The merging process starts at `h1` with puse ERA 5 values, and ends at `h2` with puse climatology values.

## Results

### `azimuth_deviation_table.txt`
This is the main summary of results. It contains a header and each row represents the results of a specific source-station pair. Each column of this file is described below:
- Column 1, `Year`: the year as YYYY. 
- Column 2, `DOY` : Day of Year as three digit numbed (1 to 356).
- Column 3, `SouNum` : source number starting from 1. Follows order from `sou_pos` from `discretize_parameters.toml`.
- Column 4, `StaNum` : station number starting from 1. Follows order from `sta_pos` from `discretize_parameters.toml`.
- Column 5, `TrueBaz` : geographical (true) backazimuth from station to particular source in degrees.
- Column 6, `BazDevS` : stratospheric-only calculated backazimuth deviation from true.
- Column 7, `#BazDS` : number of stratospheric ground intercepts found for this calculation.
- Column 8, `BazDevT` : thermospheric-only calculated backazimuth deviation from true.
- Column 9, `#BazDT` : number of thermospheric ground intercepts found for this calculation.
- Column 10, `BazDevA` : average of `BazDevS` and `BazDevT` to calculate the azimuth deviation estimation. 
- Column 11, `StdBDA` : standard deviation of `BazDevA`.
- Column 12, `Ill`: `True` if particular source-station pair could not be found, of `False` in the contrary case. A successful calculation should have `False` in this column.

## Thanks to...
- Robin Matoza: concept, analisis, and many ideas
- Jeremy Francoeur: proofreading, testing 
