# arcade
_Automated Rapid Climatological Azimuth Deviation Estimation_

![Diagram climatological range-independent model](clim-range-ind-diag.jpg?raw=true "Diagram Climatological range-independent model")

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

## Related publications
- De Negri, R., & Matoza, R. S. (2023). Rapid location of remote volcanic infrasound using 3D ray tracing and empirical climatologies: Application to the 2011 Cordón Caulle and 2015 Calbuco eruptions, Chile. Journal of Geophysical Research: Solid Earth, 128, e2022JB025735. https://doi.org/10.1029/2022JB025735
- De Negri, R. S., Rose, K. M., Matoza, R. S., Hupe, P., & Ceranna, L. (2022). Long-range multi-year infrasonic detection of eruptive activity at Mount Michael volcano, South Sandwich Islands. Geophysical Research Letters, 49, e2021GL096061. https://doi.org/10.1029/2021GL096061


## Requisites
- All external repos (see below) in place and hopefully individually tested (through makefiles inside each one). 
- Fortran and c++ compilers (check `makefile` to define Fortran compiler in first line)
  - For NRMLSIS2.0 you'll need gfortran version 4.8.5 or 6.3.0 or 8.1.0. 
    
    **Note for Ubuntu**: this means using Ubuntu 20.04 to be able to install the gfortran-8 libraries. I have not figured out a way to install older libraries in the latest stable version (22.04) without risking your workstation.
  - For HWM14 you'll need fftw3 (libfftw3-dev in Ubuntu; fftw-3 with port in MacOS).
- Python 3.9 Conda: follow instructions here https://docs.conda.io/en/latest/miniconda.html
  - Create the `arcade` envinroment with `$ conda env create -f environment_cross_platform.yml`    
  - When running the calculations, activate the environment with `$ conda activate arcade`. **NOTE:** If using hybrid models, you need to install the dependencies here too (see below).
 
- A Dockerfile and Docker image are available to ease the installation process. Please refer to the section [Using Docker](#using-docker) if you want to use this solution.

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

## Example: calculating the expected azimuth deviation at IS41 for infrasound signals coming from 2011 Puyehue-Cordón Caulle, Chile eruption
1. For this example you should keep the configuration file `./input/arcade_config.toml` untouched, except for the last line (l. 52). Set the parameter `path =  "FILLME"` with the path in your workstation where you would like to save your results (Example: `path = /Users/rodrum/Desktop/`). 
  - *Note 1:* the project name is set in line 51 as `name = "PCCVC-Clim-Norm-StratoThermo"`. This will be the folder name that will contain your results.
  - *Note 2:* `perc_cpu` (line 47) is set as `1`, which means that in theory the calculations will use all the cores in your workstation. You don't have to worry about this as this example calculates just one source to one station, which will use only one core (each profile takes one core). 
2. In the parameters file `./input/discretize_parameters.toml` you need to uncomment the coordinates of the station IS41 (or I41PY), and comment all the other predefined station locations. That is, comment line 26 and uncomment line 33. Everything else is already set-up for the example. 
3. Once the input files have been modified and saved in place, just run `make run-arcade` and wait... With a 3.5 Ghz Intel Core i7 processor, this calculation should take about 19 minutes.
4. When the calculations are completed, you should have a folder named `name` in the path defined with `path` in the config. file `./input/arcade_config.toml`.

### Results
The results are organized as follows:
- `input/` will contain the same configuration files used for the model to allow for reproducibility. It should also contain the files `doys.txt`, `sources.txt`, and `stations.txt` that are automatically generated. These files are merely used for plotting during the process and not really necessary as the values come from the files `arcade_config.toml` and `discretize_parameters.toml`.
- `output/` 
  - `azimuth_deviation_table.txt` - file containing the summary of the azimuth deviations. The colums are explained in detail in the previous section. The file should contain these two lines:  
 
    
    ```
    Year     DOY  SouNum  StaNum TrueBaz BazDevS  #BazDS BazDevT  #BazDT BazDevA   StdBDA     Ill
    2011     155       1       1  217.16    6.65      15    6.67      17    6.66 6.37e-02   False
    ```
    Meaning that the average azimuth deviation is 6.66 degrees, calculated with 15 stratospheric arrivals and 17 thermospheric arrivals.
    
  - `figures/`: this folder contains useful figures that serve to have a visual idea of the locations of the source and station as well as the atmospheric winds for each profile.
    - `glob_XXX_YYYYY_ZZZZ.png` - a global map of the location of the station number `ZZZZ` and source number `YYYYY`. In this example the numbers should be `0001` and `00001` for the station and source, respectively. `XXX` describes the day of the year number (DOY) and in this case is set as `155` (June 4 of 2011).
    ![Map of source and station location](examples/PCCVC_to_IS41-clim-norm-stratothermo/glob_155_00001_0001.png?raw=true "Global Map")
    - `prof_XXX_YYYYY_ZZZZ.png` - Top: a map that contains the last figure in the top. Middle: Temperature (K), zonal winds(E-W direction, m/s), meridional winds (N-S direction, m/s), density (g/cm^3), and pressure (mbar) in height. Bottom: along-profile (source-station) winds, and across-winds (perpendicular to source-station).
    ![Atmospheric profiles](examples/PCCVC_to_IS41-clim-norm-stratothermo/prof_155_00001_0001.png?raw=true "Atmospheric profiles")
    - `arrv_XXX_YYYYY_ZZZZ.png` - Snapshot of last iteration of ground intecept (arrivals) for the profile. Top-left: colored by Transmission Loss (dB). Top-right: colored by travel time (hours). Bottom-left: colored by celerity (km/s). Bottom-right: colored by turnin height (km).
     ![Arrivals](examples/PCCVC_to_IS41-clim-norm-stratothermo/arrv_155_00001_0001.png?raw=true "Arrivals")
  - `nodes/` contains files related with the calculations, basically middle-ground files to help obtain atmospheric descriptions with HWM14/NRMLSIS2.0
  - `proc/` contains important output of the calculations that can be used lated to a deeper analysis
    - `run_arcade_out.txt` - is the output of the main calculations. Check this file if anything goes wrong and the calculations fail. It should also contain the time it took to calculate the model in the last three lines (is the output of `time` in Unix).
    - `prof_XXX_YYYYY_ZZZZ_SOMENAME.txt` - contains the output of each iteration search done with ARCADE. The name `SOMENAME` will change depending in the kind of model you choose to calculate. In this case it will be just `_out`, as we are using raw empirical climatologies. Use this file to check how the iterations turn out.
    - `arrv/` - contains the arrivals for each iteration and profile (output of 3d rat tracing infraGA). As the number of files could be huge, this feature should be disabled for large combinations of sources and stations. The format of each file name is `NUM-A_prof_XXX_YYYYY_ZZZZ.arrivals.dat`, where `NUM` is the run number (iteration number), and `A` could be `1` or `2` depending on the launch azimuth used to enclose the target station (see PAPER TBD). Use these files to understand the iteration process, or even make a video as each files is a snapshot of the resulting 3d ray tracing calculations. Note that the plot `arrv_XXX_YYYYY_ZZZZ.png` is actually the last snapshot of the profile `XXX_YYYYY_ZZZZ`.
    - `rays/` - contains the 3D ray paths for each iteration and profile (also output of infraGA). The format is similar to the files in `arrv/`. You can use this files to plot the 3D rays.
  - `profiles/` contains the atmospheric descriptions for the source-station direction. The files `.met` are the descriptions that infraGA uses to run 3D ray tracing.

## Using Docker

We published a Docker image to simplify the installation process of ARCADE. By using Docker, you can quickly set up the required environment without worrying about the dependency installations.

### Prerequisites

Before you begin, ensure that you have Docker installed on your system. You can download and install Docker by following the instructions provided in the [Docker documentation](https://docs.docker.com/get-docker/).

### Installation Steps

1. Pull the Docker image from the GitHub Container Registry (ghcr.io) by running the following command in your terminal:
```bash
docker pull ghcr.io/rodrum/arcade:main
```
2. Once the Docker image is successfully downloaded, you can create a new container and start the software using the following command:

```bash
docker run -it --name arcade ghcr.io/rodrum/arcade:main
```
	This will start interactively a new docker container based on the previously downloaded image and named "arcade". All the dependencies are already installed and you can start using the software. Nano is installed by defaults and you can use it to modify the input files. **Note that you don't need to run `conda activate arcade` before `make run-arcade`, the environment is activated by defaults.**
3. After exiting the container, you can restart it using:

```bash
docker start -i arcade
```

### Retrieving results

To retrieve the results of ARCADE, you have two solutions. 
1. Using `docker cp`: This solution should be used if you want to retrieve results from small runs or only the result table.
2. Using bind mounts: This solution should be used if you want to retrieve all the data from the output folder in case of long runs.

#### Using `docker cp`:

In a terminal type:
```bash
docker cp arcade:/app/results/path/to/folder /path/to/destination
```

This will copy the file or folder specified to the given destination.

#### Using bind mounts:

`docker cp` can be quite  especially when copying large folders/files. We propose another solution: [bind mounts](https://docs.docker.com/storage/bind-mounts/).

**Note that bind mounts may slow down file access on MacOS and Windows with WSL2 platforms.**

Bind mounts allow the docker container to directly write on the host instead of inside the docker container. To use a bind mount, you need to run the container with the following command:

```bash
docker run -it --name arcade --mount type=bind,source=/path/on/the/host,target=/app/results ghcr.io/username/docker-image:tag
```

## Thanks to...
- Robin Matoza: concept, analisis, and many ideas
- Jeremy Francoeur: proofreading, testing 
