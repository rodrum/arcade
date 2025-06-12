# arcade
_Automated Rapid Climatological Azimuth Deviation Estimation_

![How does ARCADE work](diagram.png?raw=true "How does ARCADE work")

## Purpose
This project aims to understand the feasibility of a simplified rapid azimuth
deviation prediction scheme for improving the source location when implementing
an association and location cross-bearings, grid-search method like IMS\_vASC
[Matoza et al., 2017].
This implies using "offline", portable, and available for research purposes tools,
so the calculations could be performed by any scientist in their personal
workstation machine, but could also scale to a large research computation center.

## Project code
The project itself is a collection of Python scripts plus some Fortran 90 layers
that communicate with available open source climatologies and a 3D ray-tracing
algorithm.
`ARCADE_main.py` contains the functions that serve to find and estimate the
azimuth deviations, representing the core of the used methodology.

## Current capabilities

ARCADE is cabable of different types of ray propagation modelling. The most 
basic model is range independent (spherically layered atmosphere), meaning that 
uses average atmospheric descriptions in height for a source to receiver (station)
direction slice ('profile'). This case uses only climatologies to determine the
properties of the atmosphere for a given datetime and profile. This case can be
tuned to use "perturbed" climatological profiles (see De Negri and Matoza, 2023).

The second range independent modeling type is a hybrid calculation that merges
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

Additionally, the same range-independent/range-dependent models can now use
NCPA-G2S atmospheric specifications, thanks to the work of C. Hetzer at
LANL (see https://github.com/chetzer-ncpa/ncpag2s-clc).

## Related publications to this project
- De Negri, R. S., Matoza, R. S., Hupe, P., Le Pichon, A., Rose, K. R.,
Cevuard, S., Niroa, J. J. (2025). Evaluating the temporal capability of empirical
climatologies for rapid long-range volcanic infrasound propagation estimates
using a multidecadal data set of persistent Vanuatu volcanic eruptions.
Geophysical Journal International, 241 (1), pp. 268-290.
https://doi.org/10.1093/gji/ggaf027
- De Negri, R. and Matoza, R. S. (2023). Rapid location of remote volcanic
infrasound using 3D ray tracing and empirical climatologies: Application to the
2011 Cordón Caulle and 2015 Calbuco eruptions, Chile.
Journal of Geophysical Research: Solid Earth, 128, e2022JB025735.
https://doi.org/10.1029/2022JB025735
- De Negri, R. S., Rose, K. M., Matoza, R. S., Hupe, P., & Ceranna, L. (2022).
Long-range multi-year infrasonic detection of eruptive activity at Mount Michael
volcano, South Sandwich Islands. Geophysical Research Letters, 49, e2021GL096061.
https://doi.org/10.1029/2021GL096061


## Installation on your computer

If choosing to install on your computer, you will need to manually deal with installing the necessary dependencies.
Ensure that you have:

- Fortran and c++ compilers (check `makefile` to define Fortran compiler in first line)
  - For NRMLSIS2.0 you'll need gfortran version 4.8.5 or 6.3.0 or 8.1.0. 
    Ubuntu 24.04, the gfortran version 8.x are not available, but it is possible
    to install them via `conda`, through the `conda-forge` repos (see below).
  - For HWM14 you'll need fftw3 (libfftw3-dev in Ubuntu; fftw-3 with port in MacOS).
- Python 3.9 Conda: follow instructions here https://docs.conda.io/en/latest/miniconda.html
  - Create the `arcade` envinroment with `conda env create -f environment_cross_platform.yml`    
  - When running the calculations, activate the environment with `conda activate arcade`.
  **NOTE:** If using hybrid models, you need to install the dependencies here too (see below).
- All external repos (see below) in place (inside `repos`) and hopefully individually
  tested (through makefiles inside each one).
  
### External repositories

The climatologies used are the Horizontal Wind Model (HWM14) [Drob et al., 2014] 
to estimate horizonal winds in height, and NRLMSIS2.0 [Emmert et al., 2021] to 
estimate other properties of the atmosphere in height, like density, temperature 
and pressure. The 3D ray-tracing algorithm used is infraGA [Blom and Waxler, 2012].
Additionally, a very convenient interface to obtain the NCPA ground-to-space (NCPA-G2S)
atmospheric descriptions has been recently made available by C. Hetzer
(see https://github.com/chetzer-ncpa/ncpag2s-clc).
These descriptions are similar to the hybrid models, coupling reanalysis
datasets and empirical climatologies automatically.
They are now included as part of the possible tools to model the back-azimuth
deviation adding a more realistic set of descriptions than empirical climatologies.
Note that these need an internet conection to work.

The climatological atmospheric specifications (HWM14+NRLMSIS2.0) and the 3D ray tracing algorithm
(infraGA/geoAC) require to be downloaded and placed inside `repos` folder,
with the names as defined in the `makefile`.
That is, `./repos/HWM14`, `./repos/NRLMSIS2.0` and `./repos/infraGA-master`.
For the NCPA-G2S atmospheric specifications,
the folder name should follow the `path` variable within
`input/config.toml` file, under `[discretization.ncpag2s]`
(i.e., `ncpag2s-clc-main`)

Below are the instructions to access each one of the components.

#### NRLMSIS2.0
NRLMSIS 2.0 Code and all data samples used in this work are available at https://map.nrl.navy.mil/map/pub/nrl/NRLMSIS/NRLMSIS2.0
The associated publication is __Emmert, J. T., Drob, D. P., Picone, J. M., Siskind, D. E., Jones, M., Mlynczak, M. G., et al. (2021). NRLMSIS 2.0: A whole-atmosphere empirical model of temperature and neutral species densities. Earth and Space Science, 8, e2020EA001321. https://doi.org/10.1029/2020EA001321__


#### HWM14
Available to download from the supporting information of __Drob, D. P., Emmert, J. T., Meriwether, J. W., Makela, J. J., Doornbos, E., Conde, M., Hernandez, G., Noto, J., Zawdie, K. A., McDonald, S. E., et al. (2015), An update to the Horizontal Wind Model (HWM): The quiet time thermosphere, Earth and Space Science, 2, 301– 319, doi:10.1002/2014EA000089.__

#### infraGA/geoAC
Available from the public seosmoacoustic repository of Los Alamos National Laboratory https://github.com/LANL-Seismoacoustics/infraGA

#### NCPA-G2S
Hetzer, C.H. (2024). "The NCPAG2S command line client". https://github.com/chetzer-ncpa/ncpag2s-clc. doi: 10.5281/zenodo.13345069.

### Hybrid models
In addition to the empirical climatologies (NRMSLSIS2.0/HWM14), you can also use hybrid atmospheric descriptions that are a combination of ERA 5 ECMWF below ~80 km and empirical climatologies at higher altitudes. These models were used to compare and validate the climatological estimations, but can still be used for realistic aproaches.

#### CDSAPI Python module
The hybrid models need the Climate Data Store (CDS) infrastructure and API module, `cdsapi`, provided by the European Centre for Medium-Range Weather Forecasts (ECMWF, https://www.ecmwf.int/). Please visit https://cds.climate.copernicus.eu/api-how-to to ensure the module is installed in your system, and properly configured API user and key to be able to download the atmospheric descriptions.

#### ecCodes to decode data
In order to decode the downloaded data and use it automatically for ray tracing with infraGA, it is necessary to install `ecCodes` provided by ECMWF. Please follow the instructions here https://confluence.ecmwf.int/display/ECC/ecCodes+installation. I recommend to use the Python binding installation (https://confluence.ecmwf.int/display/UDOC/How+to+install+ecCodes+with+Python+bindings+in+conda+-+ecCodes+FAQ).

## Installation using Docker 
We published a Docker image to simplify the installation process of ARCADE.
By using Docker, you can quickly set up the required environment without worrying
about the dependency installations.
**Note:** At the moment, this image is a couple of commits behind.
However, once your docker image is up and running, you can manually download the
latest version and use it inside the docker image if you want.

### Prerequisites

Before you begin, ensure that you have Docker engine installed on your system. You can download and install Docker engine by following the instructions provided in the [Docker documentation](https://docs.docker.com/engine/install/).

### Installation Steps

1. Pull the Docker image from the GitHub Container Registry (ghcr.io) by running the following command in your terminal:
```bash
docker pull ghcr.io/rodrum/arcade:main
```
2. Once the Docker image is successfully downloaded, you can create a new container and start the software using the following command:

```bash
docker run -it --name arcade ghcr.io/rodrum/arcade:main
```
This will start interactively a new docker container based on the previously
downloaded image and named "arcade".
All the dependencies are already installed and you can start using the software.
Nano is installed by defaults and you can use it to modify the input files.
**Note that you don't need to run `conda activate arcade` before `make run-arcade`, the environment is activated by default.**

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

`docker cp` can be quite slow especially when copying large folders/files. We propose another solution: [bind mounts](https://docs.docker.com/storage/bind-mounts/).

**Note that bind mounts may slow down file access on MacOS and Windows with WSL2 platforms.**

Bind mounts allow the docker container to directly write on the host instead of inside the docker container. To use a bind mount, you need to run the container with the following command:

```bash
docker run -it --name arcade --mount type=bind,source=/path/on/the/host,target=/app/results ghcr.io/username/docker-image:tag
```

## How to run in general

1. Inside repo, run `make all-prep` (not needed if using the docker image). 
2. Modify the input configuration file `config.toml` to suit your model.
(Currently it is setup with an example run for Puyehue-Cordon Caulle to IS02 on 2011-06-04 at 19:00:00 UTC)
3. Run `make run-arcade` to calculate the azimuth deviations.
4. After the modeling ends, the input configuration files and output figures
and results will be saved in the directory `path`/`name`, which is set in
`config.toml` by the variable `path` in the section `[project_info]`.
Inside `./output`, the text file `azimuth_deviation_table.txt` has the summary of results.

## Where to look if things fail

It the run stops due to an undescribed error
(e.g., `make: *** [makefile:81: run-arcade] Error 1`)
check the output processing file in `output/proc/run_arcade_out.txt`.

## Parameters and options of `./input/config.toml`

Please check the `config.toml` file for clarifications about the parameters.

## Results

### `azimuth_deviation_table.txt`
This is the main summary of results. It contains a header and each row represents the results of a specific source-station pair. Each column of this file is described below:
- Column 1, `Year`: the year as YYYY. 
- Column 2, `DOY` : Day of Year as three digit numbed (1 to 356).
- Column 3, `Seconds`: seconds of the day (e.g., 6h=6*3600=21600)
- Column 4, `SouNum` : source number starting from 1. Follows order from `sou_pos` from `discretize_parameters.toml`.
- Column 5, `StaNum` : station number starting from 1. Follows order from `sta_pos` from `discretize_parameters.toml`.
- Column 6, `StaNam`: nominal name of station
- Column 7, `TrueBaz` : geographical (true) backazimuth from station to particular source in degrees.
- Column 8, `BazDevS` : stratospheric-only calculated backazimuth deviation from true.
- Column 9, `#BazDS` : number of stratospheric ground intercepts found for this calculation.
- Column 10, `BazDevT` : thermospheric-only calculated backazimuth deviation from true.
- Column 11, `#BazDT` : number of thermospheric ground intercepts found for this calculation.
- Column 12, `BazDevA` : average of `BazDevS` and `BazDevT` to calculate the azimuth deviation estimation.
- Column 13, `StdBDA` : standard deviation of `BazDevA`.
- Column 14, `Ill`: `True` if particular source-station pair could not be found, of `False` in the contrary case. A successful calculation should have `False` in this column.

## Example: calculating the expected azimuth deviation at IS41 for infrasound signals coming from 2011 Puyehue-Cordón Caulle, Chile eruption
As shown in the section `[discretization]` of the `config.toml` file, this is
a case that simulates infrasound waves propagatin from Yasur volcano (Vanuatu archipelago)
to the station IS22 in New Caledonia at ~400 km.
The dates are for 2011, days 55 (Feb 24) and 56 (Feb 25), at 00:00:00 UTC and
00:06:00 UTC.

For this example you should keep the configuration file `./input/config.toml`
untouched, except for the line 6. Set the parameter `path = "FILL_THIS"` with the
path in your workstation where you would like to save your results
(Example: `path = /Users/rodrum/Desktop/`).
  - *Note 1:* the project name is set in line 51 as `name = "auto"`.
    This will be pick a name automatically for your project based on the
    volcano name, year, etc (see `config.toml`). If you change `auto` for a
    different name, it will pick that instead.

Once the input files have been modified and saved in place, just run
`make run-arcade` and wait...

When the calculations are completed, you should have a folder named `name` in
the path defined with `path` in the config. file `./input/config.toml`.

### Results
The results are organized as follows:
- `input/` will contain the same configuration files used for the model to allow
for reproducibility.
- `output/` 
  - `azimuth_deviation_table.txt` - file containing the summary of the azimuth
  deviations. The colums are explained in detail in the previous section.
  The file should contain these two lines:
 
  ```
   Year     DOY Seconds  SouNum  StaNum  StaNam TrueBaz BazDevS  #BazDS BazDevT  #BazDT  BazDevA  StdBDA     Ill
   2011      55       0       1       1      IS22  43.14     nan       0    1.62      37    1.62 7.00e-02   False
   2011      56       0       1       1      IS22  43.14     nan       0    1.48      37    1.48 6.59e-02   False
   2011      56   21600       1       1      IS22  43.14    2.15      28    2.28      36    2.22 5.22e-02   False
   2011      55   21600       1       1      IS22  43.14    2.56      33    2.68      35    2.62 4.72e-02   False
  ```
    Meaning that the average azimuth deviation is 1.62 +/- 7.00e-2 and
    1.44 +/- 6.59e-2 degrees (`BazDevA`) for days 56 and 56, respectively,
    at 00:00:00 UTC (`Seconds` is 0). At 00:06:00 UTC (`Seconds` is 21600),
    these values are 2.62 +/- 4.72e-2 and 2.22 +/- 5.22e-2 for days 55 and 56,
    respectively. At 00:06:00 UTC, both back-azimuth deviations are calculated
    considering both stratospheric and thermospheric arrivals
    (`BazDevS` and `BazDevT` are not `nan`).
    
  - `figures/`: this folder contains useful figures that serve to have a visual
  idea of the locations of the source and station as well as the atmospheric
  winds for each profile. The only automatically generated figures are:
    - `arrv_XXX_YYYYY_ZZZZ.png` - Snapshot of last iteration of ground intecept
    (arrivals) for the profile. Top-left: colored by Transmission Loss (dB).
    Top-right: colored by travel time (hours).
    Bottom-left: colored by celerity (km/s).
    Bottom-right: colored by turnin height (km).
     ![Arrivals](examples/Yasur-IS22-2011-55_56_0_21600-clim-range_ind/output/figures/arrv_00000_055_00001_0001.png?raw=true "Arrivals")
    These figures will be created if all `plot_arrivals` and `save_arrivals`
    are true.
  - `nodes/` contains files related with the calculations, basically
    middle-ground files to help obtain atmospheric descriptions with HWM14/NRMLSIS2.0
  - `proc/` contains important output of the calculations that can be used later
    for debugging or analysis.
    - `run_arcade_out.txt` - is the output of the main calculations.
    Check this file if anything goes wrong and the calculations fail.
    It should also contain the time it took to calculate the model in the last
    three lines (is the output of `time` in Unix).
    - `prof_XXX_YYYYY_ZZZZ_SOMENAME.txt` - contains the output of each iteration
    search done with ARCADE. The name `SOMENAME` will change depending in the kind
    of model you choose to calculate. In this case it will be just `_out`,
    as we are using raw empirical climatologies. Use this file to check how the
    iterations turn out.
    - `arrv/` - contains the arrivals for each iteration and profile
    (output of 3d rat tracing infraGA).
    As the number of files could be huge, this feature should be disabled for
    large combinations of sources and stations.
    The format of each file name is `NUM-A_prof_XXX_YYYYY_ZZZZ.arrivals.dat`,
    where `NUM` is the run number (iteration number), and `A` could be `1` or `2`
    depending on the launch azimuth used to enclose the target station.
    Use these files to understand the iteration process, or even make a video
    as each files is a snapshot of the resulting 3d ray tracing calculations.
    Note that the plot `arrv_XXX_YYYYY_ZZZZ.png` is actually the last snapshot
    of the profile `XXX_YYYYY_ZZZZ`.
    - `rays/` - contains the 3D ray paths for each iteration and profile
    (also output of infraGA).
    The format is similar to the files in `arrv/`.
    You can use this files to plot the 3D rays.
  - `profiles/` contains the atmospheric descriptions for the source-station
    direction. The files `.met` are the descriptions that infraGA uses to run 3D
    ray tracing.


## Thanks to...
- Robin Matoza: concept, analysis, and many ideas.
- Jeremy Francoeur: proofreading, testing.
- Vincent Boulenger: dealing with docker image and making output better.
