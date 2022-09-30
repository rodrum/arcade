# arcade
_Automated Rapid Climatological Azimuth Deviation Estimation_

## Purpose
This project aims to understand the feasibility of a simplified rapid azimuth deviation prediction scheme for improving the source location when implementing an association and location cross-bearings, grid-search method like IMS\_vASC.
As a design principle, this implied using "offline", portable, and available for research purposes tools, so the calculations could be performed by any scientist in their personal workstation machine, but could also scale to a large research computation center.
The empirical climatologies allow us to skip the waiting time of accessing and downloading realistic atmospheric descriptions and ensure reproducible results for reliable scientific use.

## ARCADE with NRLMSIS-2.0 and latest infraGA

Collection of Python (and some bash) scripts, that use available open source
climatologies and 3D ray-tracing to estimate arrival parameters in a rapid
manner, especially the expected azimuth deviation due to long-range propagation
of acoustic signals through a complex moving atmosphere.

The climatologies used are the Horizontal Wind Model (HWM14) [Drob et al., 2014] 
to estimate horizonal winds in height, and NRLMSIS2.0 [Emmert et al., 2021] to 
estimate other properties of the atmosphere in height, like density, temperature 
and pressure. The 3D ray-tracing algorithm used is infraGA [Blom and Waxled, 2012].

The goal of this collection of scripts is to be feed IMS_vASC [Matoza, 2017] 
with statistical estimates that are first-order robust values of expected
atmospheric effects on the azimuth deviation. 

### Current capabilities

ARCADE is cabable of different types of ray propagation modelling. The most 
basic model is range independent (spherically layered atmosphere), meaning that 
uses average atmospheric descriptions in height for a source to receiver (station)
direction slice ('profile'). This case uses only climatologies to determine the
properties of the atmosphere for a given datetime and profile. This case can be
tuned to use "perturbed" climatological profiles.

The second range intependent modeling type is a hybrid calculation that merges
ECMWF ERA 5 descriptions for the lower ~80 km, with climatologies for the 
upper atmosphere. This case needs internet connection, and could fail if
the download process hangs. For now, pertubed profiles have not been implemented.

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

- Note on `input_clim_parameters.txt`: SEC variable for NRLMSIS2.0 gets overwritten, so I had to hardcode it into `calculate_profiles.f90` for now...  
- Note on `makefile`: for `calculate_sph_nodes` it needs to be run twice (for now)

## How to run

1. Modify the input files `arcade_config.toml` and `discretize_parameters.toml`
to suit your model.
2. Use `run_all.sh` to create the atmospheric descriptions and do ray-tracing.
3. After the modeling ends, the input configuration files and output figures
and results will be saved in the directory `path`/`name`, which is set in
`arcade_config.toml`.

### IMS_vASC

#### Inside `./src/imsvasc_calbuco/IN`
1. **Cleaned** bulletin files per station in `./BULL`.
2. List of station code (ISXX), relative path to bulletin file (`./BULL/Name_of_bulletin.txt`).
3. List stations as longitude, latitude, and station code (as "ISXX") in `./IN` (called `glob_stalist.dat`).
4. List of all active volcanoes (Smithsonian) as latitude, longitude, height (m), volcano index, and volcano name (called `./IN/volcano_list_smith.txt`).

#### Run with azimuth corrections
1. `do.ims_vassc_during_azdev`: uses `ims_vascc_azdevcor`
2. `do.ims_vassc_prior`: uses `ims_vassc_cln`
3. `do.grid_proc`: uses `ums_vassc_gridproc_ev`
4. `do.readwrite_basic_during`: uses `readwrite_imsvasscbin_basic`
5. `do.readwrite_gridproc`: uses `readwrite_imsvasscbin_gridproc`

## Hybrid models and CDSAPI

These models were used to compare and validate the climatological estimations, but can still be used for realistic aproaches.

The hybrid models need the Climate Data Store (CDS) infrastructure and API module, `cdsapi`, provided by the European Centre for Medium-Range Weather Forecasts (ECMWF, https://www.ecmwf.int/). Please visit https://cds.climate.copernicus.eu/api-how-to to ensure the module is intalled in your system, and properly configured API user and key to be able to download the atmospheric descriptions.
