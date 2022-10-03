# arcade
_Automated Rapid Climatological Azimuth Deviation Estimation_

## Purpose
This project aims to understand the feasibility of a simplified rapid azimuth deviation prediction scheme for improving the source location when implementing an association and location cross-bearings, grid-search method like IMS\_vASC [Matoza et al., 2017].
As a design principle, this implied using "offline", portable, and available for research purposes tools, so the calculations could be performed by any scientist in their personal workstation machine, but could also scale to a large research computation center.

## Project code
The project itself is a collection of Python scripts plus some Fortran 90 layers that communicate with 
available open source climatologies and a 3D ray-tracing algorithm. 
`ARCADE_main.py` contains the functions that serve to find and estimate the azimuth deviations, representing the core of the used methodology. 

## External repos

The climatologies used are the Horizontal Wind Model (HWM14) [Drob et al., 2014] 
to estimate horizonal winds in height, and NRLMSIS2.0 [Emmert et al., 2021] to 
estimate other properties of the atmosphere in height, like density, temperature 
and pressure. The 3D ray-tracing algorithm used is infraGA [Blom and Waxler, 2012].

The method has two main external components: the atmospheric model (HWM14+NRLMSIS2.0) and the 3D ray tracing algorithm (infraGA/geoAC). Each one of this components has an associated repository that needs to be inside `repos` folder, with the names as defined in the `makefile`. That is, `./repos/HWM14`, `./repos/NRLMSIS2.0` and `./repos/infraGA-master`. Below the instructions to access each one of the components.

### NRLMSIS2.0
NRLMSIS 2.0 Code and all data samples used in this work are available at https://map.nrl.navy.mil/map/pub/nrl/NRLMSIS/NRLMSIS2.0
The associated publicacion is __Emmert, J. T., Drob, D. P., Picone, J. M., Siskind, D. E., Jones, M., Mlynczak, M. G., et al. (2021). NRLMSIS 2.0: A whole-atmosphere empirical model of temperature and neutral species densities. Earth and Space Science, 8, e2020EA001321. https://doi.org/10.1029/2020EA001321__


### HWM14
Available to download from the supporting information of __Drob, D. P., Emmert, J. T., Meriwether, J. W., Makela, J. J., Doornbos, E., Conde, M., Hernandez, G., Noto, J., Zawdie, K. A., McDonald, S. E., et al. (2015), An update to the Horizontal Wind Model (HWM): The quiet time thermosphere, Earth and Space Science, 2, 301– 319, doi:10.1002/2014EA000089.__

### infraGA/geoAC
Available from the public seosmoacoustic repository of Los Alamos National Laboratory https://github.com/LANL-Seismoacoustics/infraGA

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
- Fortran and c++ compilers (check `makefile` to define Fortran compiler in first line)
- Python 3.9 Conda: follow instructions here https://docs.conda.io/en/latest/miniconda.html
  - Create an envinroment (replace `<env>` with the name you want) and install the requires packages listed in `arcade_conda_env.txt`
  `$ conda create --name <env> --file arcade_conda_env.txt`
  
## How to run

1. Inside repo, run `$ make all-prep`. 
2. Modify the input files `arcade_config.toml` and `discretize_parameters.toml`
to suit your model. 
3. Run `$ make run-arcade` to calculate the azimuth deviations.
4. After the modeling ends, the input configuration files and output figures
and results will be saved in the directory `path`/`name`, which is set in
`arcade_config.toml`. Inside `./output`, the text file `azimuth_deviation_table.txt` has the summary of results.

## Parameters and options

FILL ME

### Hybrid models and CDSAPI

These models were used to compare and validate the climatological estimations, but can still be used for realistic aproaches.

#### CDSAPI Python module
The hybrid models need the Climate Data Store (CDS) infrastructure and API module, `cdsapi`, provided by the European Centre for Medium-Range Weather Forecasts (ECMWF, https://www.ecmwf.int/). Please visit https://cds.climate.copernicus.eu/api-how-to to ensure the module is intalled in your system, and properly configured API user and key to be able to download the atmospheric descriptions.

#### ecCodes to decode data
In order to decode the downloaded data and use it automatically for ray tracig with infraGA, it is necessary to install `ecCodes` provided by ECMWF. Please follow the instructions here https://confluence.ecmwf.int/display/ECC/ecCodes+installation. I recommend to use the Python binding installation (https://confluence.ecmwf.int/display/UDOC/How+to+install+ecCodes+with+Python+bindings+in+conda+-+ecCodes+FAQ).

## Other
Tested on macOS Big Sur and Ubuntu Linux only!
