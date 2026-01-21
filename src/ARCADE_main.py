"""
A-utomated
R-apid
C-limatological
A-zimuth
D-eviation
E-stimation

This program looks for the right azimuth range to be launched by
infraga-sph
Based on "Search rays prototype.ipynb"
Notes:
- infraga-sph should be in PATH
By Rodrigo De Negri
@UCSB
2019/05/08
-> Modified adaptative version from 'search_azimuth.py' [2019/05/14]

Last modification on Thu Jul 11 14:06:19 PDT 2019:
    - parallelized the run in as default
    - use backazimuth instead of azimuths
"""

import numpy as np
import obspy as ob
import sys
import os
from os import mkdir, stat
from os.path import join, isdir
import subprocess
from subprocess import Popen, PIPE
from contextlib import redirect_stdout
import toml
import multiprocessing as mp
from itertools import product

gps2dist_azimuth = ob.geodetics.base.gps2dist_azimuth


def midpoint(phi1, phi2):
    """
    Function to calculate the distance and middle
    point of two angles.
    2019/08/15

    Parameters
    ----------
    :phi1: float
           angle in degrees
    :phi2: float
           angle in degrees
    """
    print("\n[midpoint]")
    # if any of the angles is np.nan, skip
    if np.isnan(phi1) or np.isnan(phi2):
        print("** W: Input angles were not numbers")
        print("      >> phi1={0:7.3f}, phi2={1:7.3f}".format(phi1, phi2))
        return np.nan
    else:
        mid_p = phi1 + np.abs(phi1-phi2)/2
        if mid_p < 0:
            mid_p += 360
        print(f"-- I: phi1={phi1:7.3f}, phi2={phi2:7.3f}")
        print(f"-- I: mid angle={mid_p:.1f}")
        return mid_p


def baz_dev(true_baz, calc_baz):
    """
    Parameters
    ----------
    :true_baz: float
           angle in degrees
    :calc_baz: float
           angle in degrees
    """
    print("\n[baz_dev]")
    # if any of the angles is np.nan, skip
    if np.isnan(true_baz) or np.isnan(calc_baz):
        print("** W: Input angles were not numbers")
        print("      >> true_baz={0:7.3f}, calc_baz={1:7.3f}".format(true_baz, calc_baz))
        return np.nan
    else:
        print(f"-- I: baz dev={true_baz-calc_baz:.1f}")

        # e.g., true_baz=359
        case_1a = np.cos(true_baz*np.pi/180)>0 and np.sin(true_baz*np.pi/180)<0
        # e.g., calc_baz=1
        case_1b = np.cos(calc_baz*np.pi/180)>0 and np.sin(calc_baz*np.pi/180)>0
        case_1 = case_1a and case_1b

        # e.g., true_baz=1
        case_2a = np.cos(true_baz*np.pi/180)>0 and np.sin(true_baz*np.pi/180)>0
        # e.g., calc_baz=359
        case_2b = np.cos(calc_baz*np.pi/180)>0 and np.sin(calc_baz*np.pi/180)<0
        case_2 = case_2a and case_2b

        if case_1 is np.True_:
            return true_baz - calc_baz - 360  # e.g., return -2
        elif case_2 is np.True_:
            return true_baz - calc_baz + 360  # e.g., return +2
        else:
            return true_baz - calc_baz


def get_profiles(doys, pert_flag=False, mix_flag=False):
    """
    Returns: (sta_lat, sta_lon, sou_lat, sou_lon, prof_name)
    where
    sta_lat      : station latitude
    sta_lon      : station longitude
    sou_lat      : source latitude
    sou_lon      : source longitude
    prof_name : name of file that has the ROOT name of the profile
    """
    print("\n[get_profiles] Adding profiles to list...")
    profiles = []
    secs        = config['discretization']['sec']
    #doys        = config['discretization']['doys']
    sources     = config['discretization']['sou_pos']
    stations    = config['discretization']['sta_pos']
    stations_name = config['discretization']['sta_name']

    end_str = '.met'
    if pert_flag is True and mix_flag is False:
        end_str = '_pert.met'
    elif pert_flag is False and mix_flag is True:
        end_str = '_mix.met'
    elif pert_flag is True and mix_flag is True:
        end_str = '_pert.met'
        # NOTE: this case is new, this is a test

    for sec, doy, isou, ista in product(secs, doys, range(len(sources)), range(len(stations))):
        sou_lat, sou_lon = sources[isou]
        sta_lat, sta_lon = stations[ista]
        sta_name = stations_name[ista]
        prof_name = f"prof_{sec:05d}_{doy:03d}_{isou+1:05d}_{ista+1:04d}" \
                    +f"{end_str}"
#        prof_name_1 = f"../output/profiles/1_{prof_name}"
#        if stat(prof_name_1).st_size > 0:  # if file is not empty
#            profiles.append((
#                sta_lat, sta_lon,
#                sou_lat, sou_lon,
#                sta_name, prof_name))
#            print(f"   > {prof_name} added.")
#        else:  # file is empty, default to normal case
#            prof_name = f"prof_{sec:05d}_{doy:03d}_{isou+1:05d}_{ista+1:04d}.met"
#            if pert_flag is True: # overwrite as mix
#                prof_name = f"prof_{sec:05d}_{doy:03d}_{isou+1:05d}_{ista+1:04d}_mix.met"
        profiles.append((
            sta_lat, sta_lon,
            sou_lat, sou_lon,
            sta_name, prof_name))
        # print("Using profile without perturbation")
        print(f"   > {prof_name} added.")
    return profiles

def get_profiles_rngdep(doys):
    """
    Returns: (sta_lat, sta_lon, sou_lat, sou_lon, prof_name)
    where
    sta_lat      : station latitude
    sta_lon      : station longitude
    sou_lat      : source latitude
    sou_lon      : source longitude
    prof_name : name of file that has the ROOT name of the profile

    Note: the profile number is not necessary, so here I set it up as zero in
            prof_num = 0. Later this part is removed from the call.
    """
    print("\n[get_profiles_rngdep] Adding profiles to list...")
    profiles = []
    secs        = config['discretization']['sec']
    #doys        = config['discretization']['doys']
    sources     = config['discretization']['sou_pos']
    stations    = config['discretization']['sta_pos']
    stations_name = config['discretization']['sta_name']

    end_str = '.met'
    # Below some flags I don't use right not, but will be useful later when
    # including perturbations to range dependent profiles
    # if pert_flag == True:
    #     end_str = '_pert.met'
    # elif pert_flag == False and mix_flag == True:
    #     end_str = '_mix.met'

    for sec, doy, isou, ista in product(
            secs, doys, range(len(sources)), range(len(stations))
            ):
        prof_name = f"prof_{sec:05d}_{doy:03d}_{isou+1:05d}_{ista+1:04d}_0" \
                    +f"{end_str}"
        profiles.append((
            stations[ista][0], stations[ista][1],
            sources[isou][0], sources[isou][1],
            stations_name[ista], prof_name))
        print(f"   > {prof_name} added.")
    return profiles

def create_launch(launch_parameters):
    """
    This function is called to create the launch configurations
    files for infraga-sph
    Each line contains one particular launch. In this case it will
    have just two lines.
    """
    print("[create_launch] Creating launch string...")
    launch = []
    prof = launch_parameters['prof_name']
    for azimuth, prof_ind in zip(launch_parameters['azimuth'],
                                 launch_parameters['prof_inds']):
        launch_txt =  "{0} -prop {1}"\
                    " incl_min={2}"\
                    " incl_max={3}"\
                    " incl_step={4}"\
                    " azimuth={5}"\
                    " rng_max={6}"\
                    " bounces={7}"\
                    " src_lat={8}"\
                    " src_lon={9}"\
                    " src_alt={10}"\
                    " write_atmo={11}"\
                    " calc_amp={12}"\
                    " freq={13}".format(
                        './infraga-sph',
                        '../output/profiles/'+str(prof_ind)+'_'+prof,
                        launch_parameters['incl_min'],
                        launch_parameters['incl_max'],
                        launch_parameters['incl_step'],
                        azimuth,
                        launch_parameters['rng_max'],
                        launch_parameters['bounces'],
                        launch_parameters['src_lat'],
                        launch_parameters['src_lon'],
                        launch_parameters['src_alt'],
                        launch_parameters['write_atmo'],
                        launch_parameters['calc_amp'],
                        launch_parameters['freq']
                        )
        launch.append(launch_txt)
    return launch

def create_launch_rngdep(atmo_type, launch_parameters):
    """
    This function is called to create the launch configurations
    files for infraga-sph
    Each line contains one particular launch. In this case it will
    have just two lines.
    infraga-sph-rngdep  [option]  profile_id  nodes-lat.loc  nodes-lon.loc  [parameters]

    Notes:
        - numbered profiles should be in ../output/profiles/
          named as "A_prof_BBB_CCCCC_DDDD_EEE.met", where
          A    : either 1 or 2, for the two azimuths that should contain aimed stat.
          BBB  : DOY with zeroes
          CCCCC: source number with zeros
          DDDD : station number with zeros
          EEE  : number of node without zeros (i.e., '1', not '001')
    """
    print("[create_launch_rngdep] Creating launch string...")
    launch = []
    prof = launch_parameters['prof_name'].split('_')
    print(f"-> {prof}")
    nsec = prof[1]
    ndoy = prof[2]
    isou = prof[3]
    ista = prof[4]
    prof = f"prof_{nsec}_{ndoy}_{isou}_{ista}_" # everything but prof num
    nodes_lon = 'nodes-lon.loc'
    nodes_lat = 'nodes-lat.loc'
    if atmo_type == 'ncpag2s':
        nodes_lon = 'ncpag2s/lons.dat'
        nodes_lat = 'ncpag2s/lats.dat'
    for azimuth, prof_ind in zip(launch_parameters['azimuth'],
                                 launch_parameters['prof_inds']):
        launch_txt = f"./infraga-sph-rngdep -prop ../output/profiles/{str(prof_ind)}_{prof} "\
            +f"../output/profiles/{nodes_lat} ../output/profiles/{nodes_lon} "\
            +f"incl_min={launch_parameters['incl_min']} "\
            +f"incl_max={launch_parameters['incl_max']} "\
            +f"incl_step={launch_parameters['incl_step']} "\
            +f"azimuth={azimuth} "\
            +f"rng_max={launch_parameters['rng_max']} "\
            +f"bounces={launch_parameters['bounces']} "\
            +f"src_lat={launch_parameters['src_lat']} "\
            +f"src_lon={launch_parameters['src_lon']} "\
            +f"src_alt={launch_parameters['src_alt']} "\
            +f"write_atmo={launch_parameters['write_atmo']} "\
            +f"calc_amp={launch_parameters['calc_amp']} "
        launch.append(launch_txt)
    return launch

def filter_atten(geo_atten, atmo_atten, atten_thresh_db=-120):
    """
    Function to filter results at 'atten_thresh_db'
    """
    print("[filter_atten] Called with attenuation threshold at {0:7.3f} db".format(atten_thresh_db))

    filt_grdInt = []
    index = 0
    for atten1, atten2 in zip(geo_atten, atmo_atten):
        atten = atten1+atten2
        if atten > atten_thresh_db:
            filt_grdInt.append(index)
        index += 1
    print("  Number of rays after filtering: {0:d} of {1:d}"
        .format(len(filt_grdInt), len(geo_atten)))
    return filt_grdInt

def filter_grdInt(grdInt_lats, grdInt_lons, gridLat, gridLon, staLat, staLon, dw, prof):
    """
    This function filter ground intercepts in a "ring" or
    "band" of width 'dw' with a radius 'r-d0', where 'r'
    is the distance from the grid node to the particular
    ground intercept, and 'd0' is the distance from the
    grid node to the station.
    """
    print(f"[filter_grdInt] called for profile {prof}.")
    print(f"  {len(grdInt_lats):d} ground intercepts.")
    print(f"  Threshold distance of {dw:2.4f} km.")


    filt_grdInt = []
    index = 0
    d0, _, _ = gps2dist_azimuth(gridLat, gridLon, staLat, staLon)
    for lat, lon in zip(grdInt_lats, grdInt_lons):
        r, _, _ = gps2dist_azimuth(gridLat, gridLon, lat, lon)
        if np.abs(r-d0)/1000.0 < dw:
            filt_grdInt.append((index, lat, lon))
        index += 1
    print("  Number of rays after filtering: {0:d} of {1:d}"
        .format(len(filt_grdInt), len(grdInt_lats)))
    lats = [tupi[1] for tupi in filt_grdInt]
    lons = [tupi[2] for tupi in filt_grdInt]
    print("  Means: lat={0:9.4f}, lon={1:9.4f}".format(np.mean(lats), np.mean(lons)))
    return filt_grdInt

def filter_stratoThermo(results_tab, filtered_grdInt, thresh_height, prof):
    """
    Function to filter Stratospheric and Thermospheric arrivals
    and return estimates of backazimuth
    [Thu Sep 19 12:12:34 PDT 2019]
    """
    print(f"[filter_stratoThermo] called for for profile {prof}")
    print(f"  Threshold height is {thresh_height:11.5f} km")
    # Calculate the mean of the filtered ground intercepts

    ind_height = [(tup[0], results_tab[tup[0], 7]) for tup in filtered_grdInt]
    ind_s = []
    ind_t = []
    for ind, h in ind_height:
        if h<thresh_height:
            ind_s.append(ind)
        else:
            ind_t.append(ind)

    lat_av_s, lon_av_s, baz_av_s, std_baz_av_s = 0, 0, 0, 0
    lat_av_t, lon_av_t, baz_av_t, std_baz_av_t = 0, 0, 0, 0
    # Stratospheric
    if len(ind_s) == 0:
        print("  WARNING: No stratospheric arrivals for this run.")
        print("           Returning an empty tuple.")
    else:
        lat_av_s = np.mean([results_tab[ind, 3] for ind in ind_s])
        lon_av_s = np.mean([results_tab[ind, 4] for ind in ind_s])
        baz_av_s = np.mean([results_tab[ind, 9] for ind in ind_s])
        std_baz_av_s = np.std(results_tab[ind_s, 9])
    # Thermospheric
    if len(ind_t) == 0:
        print("  WARNING: No thermospheric arrivals for this run.")
        print("           Returning an empty tuple.")
    else:
        lat_av_t = np.mean([results_tab[ind, 3] for ind in ind_t])
        lon_av_t = np.mean([results_tab[ind, 4] for ind in ind_t])
        baz_av_t = np.mean([results_tab[ind, 9] for ind in ind_t])
        std_baz_av_t = np.std(results_tab[ind_t, 9])
    strato_tup = (len(ind_s), lat_av_s, lon_av_s, baz_av_s, std_baz_av_s)
    thermo_tup = (len(ind_t), lat_av_t, lon_av_t, baz_av_t, std_baz_av_t)
    if len(ind_s) != 0:
        print("  Stratospheric arrivals:")
        print(f"    n   = {strato_tup[0]:d}")
        print(f"    lat = {strato_tup[1]:9.4f}")
        print(f"    lon = {strato_tup[2]:9.4f}")
        print(f"    baz = {strato_tup[3]:9.4f}")
        print(f"    std baz = {strato_tup[4]:9.4f}")
    if len(ind_t) != 0:
        print("  Thermospheric arrivals:")
        print(f"    n   = {thermo_tup[0]:d}")
        print(f"    lat = {thermo_tup[1]:9.4f}")
        print(f"    lon = {thermo_tup[2]:9.4f}")
        print(f"    baz = {thermo_tup[3]:9.4f}")
        print(f"    std baz = {thermo_tup[4]:9.4f}")
    return strato_tup, thermo_tup

def run_launch(launch_parameters, grid_params, run_num, prof_ind, atmo_type, rngdep=False, atten_th=-120, save_arrivals=False, save_raypaths=False):
    # Create launch
    launch = ''
    if rngdep is False:
        launch = create_launch(launch_parameters)
    else:
        launch = create_launch_rngdep(atmo_type, launch_parameters)
    # Run
    print("\n[run_launch] This is run number {0} for profile {1}".format(run_num, prof_ind))
    for launchi in launch:
        print("  Calculating: \n    {}".format(launchi))
        outs, errs = Popen(['bash', '-c', launchi],
                             stdout=PIPE, stderr=PIPE).communicate()
    run_num += 1
    # Load results generate with first run
    # NOTE: profiles are fixed, they should change
    prof_1, prof_2 = '', ''
    rays_1_out, rays_2_out = '', ''
    out_results = '../output/profiles/'
    if rngdep is False:
        prof_1 = '1_'+launch_parameters['prof_name'].split('.')[0]+'.arrivals.dat'
        prof_2 = '2_'+launch_parameters['prof_name'].split('.')[0]+'.arrivals.dat'
        # To save input for processing values
        rays_1_out = '1_'+launch_parameters['prof_name'].split('.')[0]+'.raypaths.dat'
        rays_2_out = '2_'+launch_parameters['prof_name'].split('.')[0]+'.raypaths.dat'
    else:
        prof = launch_parameters['prof_name'].split('_')
        # name of profile: prof_00000_001_00001_00001
        prof = prof[0]+"_"+prof[1]+"_"+prof[2]+"_"+prof[3]+"_"+prof[4]+"_" # everything but prof num
        prof_1 = '1_'+prof+'.arrivals.dat'
        prof_2 = '2_'+prof+'.arrivals.dat'
        # To save input for processing values
        rays_1_out = '1_'+prof+'.raypaths.dat'
        rays_2_out = '2_'+prof+'.raypaths.dat'

    if save_arrivals is True:
        # Copy output to folders that will be saved, adding the run number to name
        # to keep track of the process.
        subprocess.check_output(
                            [
                                'cp',
                                out_results+prof_1,
                                '../output/proc/arrv/'+f"{run_num-1}-"+\
                                prof_1
                            ]
                        )
        print(f"  -> File {prof_1} copied to ../output/proc/arrv/")
        subprocess.check_output(
                            [
                                'cp',
                                out_results+prof_2,
                                '../output/proc/arrv/'+f"{run_num-1}-"+\
                                prof_2
                            ]
                        )
        print(f"  -> File {prof_2} copied to ../output/proc/arrv/")
    else:
        print("--> Note: intermediate-process arrivals will not be saved.")
    if save_raypaths is True:
        subprocess.check_output(
                            [
                                'cp',
                                out_results+rays_1_out,
                                '../output/proc/rays/'+f"{run_num-1}-"+\
                                rays_1_out
                            ]
                        )
        print(f"  -> File {rays_1_out} copied to ../output/proc/rays/")
        subprocess.check_output(
                            [
                                'cp',
                                out_results+rays_2_out,
                                '../output/proc/rays/'+f"{run_num-1}-"+\
                                rays_2_out
                            ]
                        )
        print(f"  -> File {rays_2_out} copied to ../output/proc/rays/")
    else:
        print("--> Note: intermediate-process raypaths will not be saved.")

    # Filter ground intercepts near the station with a radius of 'dw'
    gridLat = grid_params[0]
    gridLon = grid_params[1]
    staLat = grid_params[2]
    staLon = grid_params[3]
    dw = grid_params[4]
    results1 = np.array([])
    results2 = np.array([])
    filt_grdInt1 = []
    filt_grdInt2 = []
    # results1 = np.loadtxt(join(out_results, prof_1))
    # print(results1)
    try:
        results1 = np.loadtxt(join(out_results, prof_1))
        print(f"  Profile {prof_1} loaded.")
        results2 = np.loadtxt(join(out_results, prof_2))
        print(f"  Profile {prof_2} loaded.")
    except FileNotFoundError:
        print(f"  Profile {join(out_results, prof_1)} not found.")
        print(f"  Profile {join(out_results, prof_2)} not found.")

    try:
        ind1 = filter_atten(
                    results1[:,10], results1[:,11], atten_thresh_db=atten_th
                )
        ind2 = filter_atten(
                    results2[:,10], results2[:,11], atten_thresh_db=atten_th
                )
        print(f"DEBUG: after filter_atten")
        print(f"        ind1={len(ind1)}, ind2={len(ind2)}")
        filt_grdInt1 = filter_grdInt(results1[ind1, 3], results1[ind1, 4],
                                     gridLat, gridLon, staLat, staLon, dw,
                                     prof_1)
        filt_grdInt2 = filter_grdInt(results2[ind2, 3], results2[ind2, 4],
                                     gridLat, gridLon, staLat, staLon, dw,
                                     prof_2)
        print(f"DEBUG: after filter_grdInt")
        print(f"len(filt_grdInt1)={len(filt_grdInt1)}")
        print(f"len(filt_grdInt2)={len(filt_grdInt2)}")
    except:
        pass
    # Calculate the mean of the filtered ground intercepts
    # and the mean of the backazimuth
    # separated by stratospheric and thermospheric arrivals
    thresh_height = 60.0 # km
    #NOTE: better declare this outside for more flexibility
    strato_tup1, thermo_tup1 = filter_stratoThermo(
            results1, filt_grdInt1, thresh_height, prof_1
            )
    strato_tup2, thermo_tup2 = filter_stratoThermo(
            results2, filt_grdInt2, thresh_height, prof_2
            )
    return run_num, strato_tup1, thermo_tup1, strato_tup2, thermo_tup2

def calculate_profiles(my_profiles, arcade_conf, atmo_type, profInd=0, rngdep=False):
    '''
    # Load, and copy duplicates of profiles (.met files)
    # Obtained profile is called 'profile_00001.met', then copied
    # into two files:  'profile_DOY_XXXXX_XXXX.met_1', and
    #                  'profile_DOY_XXXXX_XXXX.met_2'.
    # This is to be able to save two different results files (different names)
    # with infraga-sph.
    '''

    # Test profile name
    prof_name = ''
    out_file_name = ''
    if rngdep is False:
        prof_name = my_profiles[profInd][5]
        out_file_name = f"{prof_name.split('.')[0]}_out.txt"
    else:
        prof = my_profiles[profInd][5].split('_')
        prof_name = prof[0]+"_"+prof[1]+"_"+prof[2]+"_"+prof[3]+"_"+prof[4] # everything but prof num
        out_file_name = f"{prof_name}out.txt"
    # Out variables
    bphi1_s, bphi1_t = 0, 0
    bphi2_s, bphi2_t = 0, 0
    num_s_tot, num_t_tot = 0, 0
    # separated by stratospheric and thermospheric arrivals
    num1_s, num1_t = 0, 0
    num2_s, num2_t = 0, 0
    with open("../output/proc/"+out_file_name, 'w') as fout:
        with redirect_stdout(fout):
            print(f"\n[calculate_profiles] for {prof_name}")
            print("  Profile number: ", profInd)
            staLat = my_profiles[profInd][0]
            staLon = my_profiles[profInd][1]
            # Test grid node
            gridLat = my_profiles[profInd][2]
            gridLon = my_profiles[profInd][3]

            sta_name = my_profiles[profInd][4]

            # ==================================================================
            # Find distance (meters), azimuth, and back azimuth
            # from grid node to current station
            # ==================================================================
            d, az, baz = gps2dist_azimuth(gridLat, gridLon, staLat, staLon)
            # Pass distance to km
            d /= 1000.0
            # Print to be sure
            print("\nSetup source and receiver parameters")
            print("====================================\n")
            print("Source            : ({0}, {1})".format(gridLat, gridLon))
            print("Receiver          : ({0}, {1})".format(staLat, staLon))
            print("Name              : {}".format(sta_name))
            print("Distance [km]     : {}".format(d))
            print("Azimuth [deg]     : {}".format(az))
            print("Back Azimuth [deg]: {}".format(baz))

            # ==================================================================
            # Define first launch parameters
            # Here the launch parameters to be used by infraga-sph are defined.
            # This will be the first run parameters.
            # ==================================================================
            print("\nLaunch parameters")
            print("=================\n")
            dphi = config['launch']['dphi']
            phi1 = az-dphi
            phi2 = az+dphi
            launch_parameters = {
                'incl_min'  : config['launch']['incl_min'],
                'incl_max'  : config['launch']['incl_max'],
                'incl_step' : config['launch']['incl_step'],
                'azimuth'   : [phi1, phi2],
                'prof_inds' :[1,2],
                'prof_name' : prof_name,
                'rng_max'   : d+config['launch']['drng'],
                'bounces'   : config['launch']['bounces'],
                'src_lat'   : gridLat,
                'src_lon'   : gridLon,
                'src_alt'   : config['discretization']['sou_alt'][0],  # NOTE: one source
                'write_atmo': config['launch']['write_atmo'],
                'calc_amp'  : config['launch']['calc_amp'],
                'freq'      : config['launch']['freq']
                }
            for key in launch_parameters:
                print("{0:11}:{1}".format(key, launch_parameters[key]))

            # ==================================================================
            # Main Loop, where the search is performed
            # ==================================================================
            bisect = [True, True]  # The first element refers to Phi1,
                                   # the second to Phi2.
                                   # "True" means "bisect it"
            # dw = 0.5*111.19        # Threshold distance to consider ground intercepts
                                   # around station. Units in km
            dw = config['launch']['min_dist_arrv']
            grid_params = [gridLat, gridLon, staLat, staLon, dw]  # First profile
            daz = config['launch']['daz']              # Threshold distance in degrees.
                                   # Used to determine if average of filtered ground
                                   # intercepts is near enough the station.
                                   # Set up as 0.5, but could be <= dw in degrees
                                   # NOTE: I think I could reduce this variable and
                                   #       make it be the same as 'dw'
                                   #       (Fri Jul 12 15:52:57 PDT 2019)
            run_num = 1            # Auxiliary variable to keep track of the run num.
            max_run = config['launch']['max_run']          # Maximum possible number of runs.
                                   # If exceeded, end process without converging.
                                   # NOTE: I should append a "special" line to the
                                   #       table in this case, or create another
                                   #       "ill" table for reference.
                                   #       (Fri Jul 12 15:53:15 PDT 2019)
            phi_min, phi_max = 0.0, 0.0  # auxiliary variables to keep track of
                                         # the min and max launch azimuths
            ill_flag = False       # Flag to know if case falls into the two
                                   # possible ill cases
            use_thermo = config['launch']['use_thermo']
                                    # By default, do not consider thermospheric
                                    # arrivals to calculate the backazimuth
            atten_th = config['launch']['atten_th']

            save_raypaths = config['launch']['save_raypaths']
            save_arrivals = config['launch']['save_arrivals']
            while True in bisect:
                # Obtain ground intercepts
                run_num, strato_tup1, thermo_tup1, strato_tup2, thermo_tup2 =\
                    run_launch(
                        launch_parameters,
                        grid_params,
                        run_num,
                        profInd,
                        atmo_type,
                        rngdep,
                        atten_th,
                        save_arrivals,
                        save_raypaths
                        )
                num_s, num_t = strato_tup1[0], thermo_tup1[0]
                num_s2, num_t2 = strato_tup2[0], thermo_tup2[0]
                lat_av1, lon_av1, lat_av2, lon_av2 = 0, 0, 0, 0
                if use_thermo is False:
                    print("  Using only stratospheric arrivals.")
                    if num_s>0:
                        lat_av1 = strato_tup1[1]
                        lon_av1 = strato_tup1[2]
                    else:
                        lat_av1 = np.nan
                        lon_av1 = np.nan
                    if num_s2>0:
                        lat_av2 = strato_tup2[1]
                        lon_av2 = strato_tup2[2]
                    else:
                        lat_av2 = np.nan
                        lon_av2 = np.nan
                elif use_thermo is True:
                    print("  Using stratospheric and thermospheric arrivals.")
                    if num_s>0 and num_t>0:
                        lat_av1 = (num_s*strato_tup1[1]+num_t*thermo_tup1[1])/(num_s+num_t)
                        lon_av1 = (num_s*strato_tup1[2]+num_t*thermo_tup1[2])/(num_s+num_t)
                    elif num_s>0:
                        lat_av1 = strato_tup1[1]
                        lon_av1 = strato_tup1[2]
                    elif num_t>0:
                        lat_av1 = thermo_tup1[1]
                        lon_av1 = thermo_tup1[2]
                    else:
                        lat_av1 = np.nan
                        lon_av1 = np.nan
                    if num_s2>0 and num_t2>0:
                        lat_av2 = (num_s2*strato_tup2[1]+num_t2*thermo_tup2[1])/(num_s2+num_t2)
                        lon_av2 = (num_s2*strato_tup2[2]+num_t2*thermo_tup2[2])/(num_s2+num_t2)
                    elif num_s2>0:
                        lat_av2 = strato_tup2[1]
                        lon_av2 = strato_tup2[2]
                    elif num_t2>0:
                        lat_av2 = thermo_tup2[1]
                        lon_av2 = thermo_tup2[2]
                    else:
                        lat_av2 = np.nan
                        lon_av2 = np.nan
                print("    Ray 1: "
                      f"\n      lat_av1={lat_av1:8.3f}"
                      f"\n      lon_av1={lon_av1:8.3f}")
                print(f"    Ray 2: "
                      f"\n      lat_av2={lat_av2:8.3f}"
                      f"\n      lon_av2={lon_av2:8.3f}")

                if run_num >= max_run:
                    print("  WARNING: The maximum number of runs has been "
                                     "exceeded, stopping calculations.")
                    phi_min = 0.0
                    phi_max = 0.0
                    bisect = [False, False]
                    ill_flag = True
                    print("\n Final azimuths:")
                    print("================")
                    print(f"min: {phi_min:7.3f}")
                    print(f"max: {phi_max:7.3f}\n")
                    break
                elif True in np.isnan([lat_av1, lon_av1, lat_av2, lon_av2]):
                    # Case to catch possible shadow zones
                    # Re-run is there are no arrivals, with bigger dw threshold
                    # until it's possible to see something
                    print("  WARNING: No ground intercepts for profile")
                    phi_min = float('nan')  # phi_min and phi_max are
                    phi_max = float('nan')  # set up as python 'nan' values
                    bisect = [False, False]
                    ill_flag = True
                    break

                # ==============================================================
                # Decision scheme
                # ==============================================================
                # Obtain apparent azimuth with ground intercepts that fall in area
                _, az1, _ = gps2dist_azimuth(gridLat, gridLon, lat_av1, lon_av1)
                _, az2, _ = gps2dist_azimuth(gridLat, gridLon, lat_av2, lon_av2)
                print("  Apparent azimuths are: ")
                print(f"    az1={az1:7.3f}")
                print(f"    az2={az2:7.3f}")
                # Now, take care of the cases near 0 degrees
                # - First, if az1 > az + 180 means that the true azimuth is near 0, and
                #   az1 is less than 0, or near 360 degrees. In this case is better to
                #   take out 360 degrees from az1 to make it slightly negative
                if az1-az>180:
                    az1 = az1-360.0
                # - If az2 < az - 270 means that az2 is slightly bigger than zero,
                #   while az is slightly less than zero, or near 360 degrees. In this
                #   case I add 360 to az2 to make it slightly bigger than 360 degrees
                elif az-az2>180:
                    az2 = 360.0+az2
                # After taking care of those ill cases, check if az1 and az2 contain
                # az like az1<az<az2
                # - This case means that az<az1<az2, so we need to make the launch
                #   angle phi1 smaller to make az1 smaller
                if (az1>az):
                    print("  Not enclosing station (phi1>true_phi)")
                    # NOTE: this criterion is really intuitive and should be
                    #       established with more rigurosity.
                    # Criterion 0: Subtract twice the difference of az1 and az to
                    #              the launch angle phi1
                    phi1 = phi1-2*(az1-az)
                    print(f"    New phi1={phi1:7.3f}.")
                    # Update the launch parameters to consider the new phi1
                    launch_parameters['azimuth'] = [phi1, phi2]
                    launch_parameters['prof_inds'] = [1, 2]
                # - This case means that az1<az2<az, so we need to make the launch
                #   angle phi2 bigger to make az2 bigger
                elif (az2<az):
                    print("  Not enclosing station (phi2<true_phi)")
                    phi2 = phi2+2*(az-az2)
                    print(f"    New phi2={phi2:7.3f}")
                    # Update the launch parameters to consider the new phi2
                    launch_parameters['azimuth'] = [phi1, phi2]
                    launch_parameters['prof_inds'] = [1, 2]
                # - Case when az1<az<az2 as expected, so proceed
                else:
                    print("  Enclosing station (phi1<az<phi2), proceeding...")
                    # Get phi1 and phi2 from launch_parameters, keeping the one
                    # that has not been changed
                    if bisect == [True, True]:
                        print("    Bisecting phi1 and phi2.")
                        phi1 = launch_parameters['azimuth'][0]
                        phi2 = launch_parameters['azimuth'][1]
                    elif bisect == [True, False]:
                        print("    Bisecting phi1.")
                        phi1 = launch_parameters['azimuth'][0]
                    elif bisect == [False, True]:
                        print("    Bisecting phi2.")
                        phi2 = launch_parameters['azimuth'][0]
                    print("  Current launch angles are:")
                    print(f"    phi1={phi1:7.3f}")
                    print(f"    phi2={phi2:7.3f}")
                    # Calculate distance to average to decide which one needs
                    # to be bisected
                    az_dist1 = np.sqrt((lat_av1-staLat)**2+(lon_av1-staLon)**2)
                    az_dist2 = np.sqrt((lat_av2-staLat)**2+(lon_av2-staLon)**2)
                    print("  Distances from ground intercepts to station:")
                    print(f"    az_dist1={az_dist1:5.2f}")
                    print(f"    az_dist2={az_dist2:5.2f}")
                    print(f"  while daz={daz}.")
                    # Set up bisect as False, adding True in case of bisection
                    bisect = [False, False]
                    dphi = phi2-phi1  # angular distance between apparent azimuths
                    scale = config['launch']['scale']         # parameter to weight dphi. NOTE: arbitrary
                    # - If distance is bigger than threshold for phi1 ground
                    #   intercepts, "bisect" it, which in reality is just "reduce"
                    #   by a certain value. In this case is slow on purpose, to be
                    #   sure it will converge. NOTE: this should be studied to make
                    #   a more rigorous rule
                    if az_dist1 > daz:
                        phi1 = phi1+dphi/scale
                        bisect[0] = True
                        print(f"  New phi1={phi1:7.3f}")
                    # - The same for phi1 ground intercepts
                    if az_dist2 > daz:
                        phi2 = phi2-dphi/scale
                        bisect[1] = True
                        print(f"  New phi2={phi2:7.3f}")
                    # If both True is not in bisect, means that the distances are
                    # now less than the threshold, so we found a first approximation
                    # for the launch angles that will produce arrivals near the station.
                    # Now, we want to be as near as the station as possible, that means
                    # "optimizing" the launch angles until we reach a better maximum
                    # distance.
                    if True not in bisect:
                        print("  Both launching azimuths have arrivals!")
                        # Optimized list was set up outside as False, False
                        # unless there are no ground intercepts or the number of tries
                        # exceeds the threshold.
                        thresh = config['launch']['thresh']  # this parameter is the new threshold to decide
                                       # when to declare that the ground intercepts
                                       # are near enough to the station
                        k = config['launch']['k']        # this parameter serves as the "1/scale" above,
                                       # that is, it will be the weight to reduce phi1/phi2
                        dphi_opt = (phi2-phi1)*0.5*k  # this will augment or diminish
                                                      # phi1/phi2.
                                                      # The value 0.5 is the ideal case
                                                      # of bisection.
                        # We will try to make the ground intercepts symmetric around
                        # the station before saving the average. The idea is to avoid
                        # bias for the backazimuth as we use phi1 and phi2 as the
                        # extreme values to stablish that the launch azimuth is
                        # (phi1+phi2)/2
                        while np.max([az_dist1, az_dist2])> thresh:
                            print("    Distances from ground intercepts to station:")
                            print(f"     az_dist1={az_dist1:5.2f}")
                            print(f"     az_dist2={az_dist2:5.2f}")
                            print(f"   while thresh={thresh}.")
                            if az_dist1>az_dist2:
                                phi1 = phi1 + dphi_opt
                                print(f"    phi1 has been augmented by {dphi_opt:8.3f}:")
                                print(f"      phi1={phi1:7.3f}")
                                print(f"      phi2={phi2:7.3f}")
                                # Below we just save phi1 and set up the prof int for the first profile
                                # to be calculated
                                launch_parameters['azimuth'] = [phi1]
                                launch_parameters['prof_inds'] = [1]
                                # Calculate new ground intercepts! Note that lat_av2 and lon_av2 will
                                # be the same. There will not be a new run for phi2.
                                run_num, strato_tup1, thermo_tup1, strato_tup2, thermo_tup2 =\
                                        run_launch(
                                            launch_parameters,
                                            grid_params,
                                            run_num,
                                            profInd,
                                            atmo_type,
                                            rngdep,
                                            atten_th,
                                            save_arrivals,
                                            save_raypaths
                                            )
                                if run_num >= max_run:
                                    print(f"    WARNING: max_run({max_run}) reached while distance is bigger than threshold")
                                    print("              Stopping calculations and saving...")
                                    break
                                num_s, num_t = strato_tup1[0], thermo_tup1[0]
                                lat_av1, lon_av1 = 0, 0
                                if use_thermo is False:
                                    print("    Using only stratospheric arrivals.")
                                    if num_s>0:
                                        lat_av1 = strato_tup1[1]
                                        lon_av1 = strato_tup1[2]
                                    else:
                                        lat_av1 = np.nan
                                        lon_av1 = np.nan
                                else:
                                    print("    Using stratospheric and thermospheric arrivals.")
                                    if num_s>0 and num_t>0:
                                        lat_av1 = (num_s*strato_tup1[1]+num_t*thermo_tup1[1])/(num_s+num_t)
                                        lon_av1 = (num_s*strato_tup1[2]+num_t*thermo_tup1[2])/(num_s+num_t)
                                    elif num_s>0:
                                        lat_av1 = strato_tup1[1]
                                        lon_av1 = strato_tup1[2]
                                    elif num_t>0:
                                        lat_av1 = thermo_tup1[1]
                                        lon_av1 = thermo_tup1[2]
                                # Check if distance fulfills our criterion
                                az_dist1 = np.sqrt((lat_av1-staLat)**2+(lon_av1-staLon)**2)
                            else:
                                phi2 = phi2 - dphi_opt
                                print(f"    phi2 has been reduced by {dphi_opt:8.3f}")
                                print(f"      phi1={phi1:7.3f}")
                                print(f"      phi2={phi2:7.3f}")
                                launch_parameters['azimuth'] = [phi2]
                                launch_parameters['prof_inds'] = [2]
                                run_num, strato_tup1, thermo_tup1, strato_tup2, thermo_tup2 = \
                                        run_launch(
                                            launch_parameters,
                                            grid_params,
                                            run_num,
                                            profInd,
                                            atmo_type,
                                            rngdep,
                                            atten_th,
                                            save_arrivals,
                                            save_raypaths
                                            )
                                if run_num >= max_run:
                                    print(f"    WARNING: max_run({max_run}) reached while distance is bigger than threshold")
                                    print("              Stopping calculations and saving...")
                                    break
                                num_s, num_t = strato_tup2[0], thermo_tup2[0]
                                lat_av2, lon_av2 = 0, 0
                                if use_thermo is False:
                                    print("    Using only stratospheric arrivals.")
                                    if num_s>0:
                                        lat_av2 = strato_tup2[1]
                                        lon_av2 = strato_tup2[2]
                                    else:
                                        lat_av2 = np.nan
                                        lon_av2 = np.nan
                                else:
                                    print("    Using stratospheric and thermospheric arrivals.")
                                    if num_s>0 and num_t>0:
                                        lat_av2 = (num_s*strato_tup2[1]+num_t*thermo_tup2[1])/(num_s+num_t)
                                        lon_av2 = (num_s*strato_tup2[2]+num_t*thermo_tup2[2])/(num_s+num_t)
                                    elif num_s>0:
                                        lat_av2 = strato_tup2[1]
                                        lon_av2 = strato_tup2[2]
                                    elif num_t>0:
                                        lat_av2 = thermo_tup2[1]
                                        lon_av2 = thermo_tup2[2]
                                az_dist2 = np.sqrt((lat_av2-staLat)**2+(lon_av2-staLon)**2)
                            dphi_opt = (phi2-phi1)*0.5*k  # update perturbation

                        # Final launch azimuths below:
                        phi_min = phi1
                        phi_max = phi2
                        # Final backazimuth
                        # separated by stratospheric and thermospheric arrivals
                        num1_s, num1_t = strato_tup1[0], thermo_tup1[0]
                        num2_s, num2_t = strato_tup2[0], thermo_tup2[0]
                        # So the idea below is to have a flag in the table as NaN
                        # to know if there are strato/thermo arrivals or not
                        bphi1_s = strato_tup1[3] if num1_s != 0 else np.nan
                        bphi1_t = thermo_tup1[3] if num1_t != 0 else np.nan
                        bphi2_s = strato_tup2[3] if num2_s != 0 else np.nan
                        bphi2_t = thermo_tup2[3] if num2_t != 0 else np.nan
                        num_s_tot = num1_s+num2_s
                        num_t_tot = num1_t+num2_t
                        print("\nFinal azimuths:")
                        print("===============")
                        print("phi1 (min): {0:7.3f}".format(phi_min))
                        print("phi2 (max): {0:7.3f}".format(phi_max))
                        print("back phi1 strato (min): {0:7.3f}".format(bphi1_s))
                        print("back phi2 strato (max): {0:7.3f}".format(bphi2_s))
                        print("back phi1 thermo (min): {0:7.3f}".format(bphi1_t))
                        print("back phi2 thermo (max): {0:7.3f}".format(bphi2_t))
                    # Below: other possible cases, until both launch azimuths generate arrivals
                    elif bisect == [True, True]:
                        print("  Both azimuths have been bisected.")
                        launch_parameters['azimuth'] = [phi1, phi2]
                        launch_parameters['prof_inds'] = [1, 2]
                        print(f"    phi1={phi1:7.3f}")
                        print(f"    phi2={phi2:7.3f}")
                    elif bisect == [True, False]:
                        print("  phi1 has been bisected")
                        launch_parameters['azimuth'] = [phi1]
                        launch_parameters['prof_inds'] = [1]
                        print(f"    phi1={phi1:7.3f}")
                        print(f"    phi2={phi2:7.3f}")
                    elif bisect == [False, True]:
                        print("  phi2 has been bisected.")
                        launch_parameters['azimuth'] = [phi2]
                        launch_parameters['prof_inds'] = [2]
                        print(f"    phi1={phi1:7.3f}")
                        print(f"    phi2={phi2:7.3f}")

            # ==================================================================
            # Save the table
            # ==================================================================
            if True not in bisect:
                print("\n========================")
                print("Saving azimuth deviation")
                print("========================")
                # deviated_az = (phi_min+phi_max)*0.5
                deviated_baz_s = midpoint(bphi1_s, bphi2_s)
                deviated_baz_t = midpoint(bphi1_t, bphi2_t)
                # true baz - observed (deviated) baz
                baz_dev_s = baz_dev(baz, deviated_baz_s)
                baz_dev_t = baz_dev(baz, deviated_baz_t)
                # Standard deviations
                std_av_s_1 = strato_tup1[4]
                std_av_s_2 = strato_tup2[4]
                std_av_t_1 = thermo_tup1[4]
                std_av_t_2 = thermo_tup1[4]
                #s_arr_flag = num_s_tot>0 and not np.isnan(baz_dev_s)
                #t_arr_flag = num_t_tot>0 and not np.isnan(baz_dev_t)
                s_arr_flag = num_s_tot>0
                t_arr_flag = num_t_tot>0
                if s_arr_flag and t_arr_flag:
                    tot_baz_dev = (num_s_tot*baz_dev_s+num_t_tot*baz_dev_t)/(num_s_tot+num_t_tot)
                    std_av_s = np.sqrt((num1_s*std_av_s_1**2+num2_s*std_av_s_2**2)/(num1_s+num2_s))
                    std_av_t = np.sqrt((num1_t*std_av_t_1**2+num2_t*std_av_t_2**2)/(num1_t+num2_t))
                    std_av = np.sqrt((num_s_tot*std_av_s**2+num_t_tot*std_av_t**2)/(num_s_tot+num_t_tot))
                elif s_arr_flag:
                    tot_baz_dev = baz_dev_s
                    std_av = np.sqrt((num1_s*std_av_s_1**2+num2_s*std_av_s_2**2)/(num1_s+num2_s))
                elif t_arr_flag:
                    tot_baz_dev = baz_dev_t
                    std_av = np.sqrt((num1_t*std_av_t_1**2+num2_t*std_av_t_2**2)/(num1_t+num2_t))
                else:
                    tot_baz_dev = np.nan
                    std_av = np.nan
                table_name = '../output/azimuth_deviation_table.txt'

                year = config['discretization']['year'][0]  # NOTE: one year

                # E.g., prof_21600_155_00001_0001
                #            5   9 1 3 5   9 1  4
                #       0123456789012345678901234
                sec     = int(prof_name[5:10])
                doy     = int(prof_name[11:14])
                gridNum = int(prof_name[15:20])
                staNum  = int(prof_name[21:25])
                if not os.path.isfile(table_name):
                    with open(table_name, "w") as f:
                        header_lst = ["Year", "DOY", "Seconds", "SouNum",
                            "StaNum", "StaNam", "TrueBaz", "BazDevS", "#BazDS",
                            "BazDevT", "#BazDT", "BazDevA", "StdBDA", "Ill"]
                        header_str = ""
                        for i, head in enumerate(header_lst):
                            if i == 0:
                                header_str += f"{head:>7}"
                            elif i == 11:
                                # scientific notation
                                header_str += f" {head:>8}"
                            elif i == len(header_lst)-1:
                                header_str += f" {head:>7}\n"
                            else:
                                header_str += f" {head:>7}"
                        print("\n",header_str,"\n")
                        f.write(header_str)

                with open(table_name, "a") as f:
                    f.write(f"{year:>7} "
                            f"{doy:>7d} "
                            f"{sec:>7d} "
                            f"{gridNum:>7d} "
                            f"{staNum:>7d} "
                            f"{sta_name:>9}"
                            f"{baz:>7.2f} "
                            f"{float(baz_dev_s):>7.2f} "
                            f"{num1_s+num2_s:>7d} "
                            f"{float(baz_dev_t):>7.2f} "
                            f"{num1_t+num2_t:>7d} "
                            f"{float(tot_baz_dev):>7.2f} "
                            f"{float(std_av):>8.2e} "
                            f"{str(ill_flag):>7}\n"
                            )

def wrap_and_save(config):
    project_name = config['project_info']['name']
    if project_name == 'auto':
        sou_name = config['discretization']['sou_name']
        if len(sou_name)>1:
            sou_name = sou_name[0]+'etal'
        else:
            sou_name = sou_name[0]

        sta_name = config['discretization']['sta_name']
        if len(sta_name)>=3:
            sta_name = sta_name[0]+"_"+sta_name[1]+"_"+sta_name[2]+"etal"
        elif len(sta_name) == 2:
            sta_name = sta_name[0]+"_"+sta_name[1]
        elif len(sta_name) == 1:
            sta_name = sta_name[0]

        year = config['discretization']['year']
        if len(year)>2:
            year = f"{year[0]}_to_{year[-1]}"
        elif len(year)==2:
            year = f"{year[0]}_{year[1]}"
        else:
            year = f"{year[0]}"

        doys = config['discretization']['doys']
        if len(doys)>3:
            doys = f"{doys[0]}_to_{doys[-1]}"
        elif len(doys) == 3:
            doys = f"{doys[0]}_{doys[1]}_{doys[2]}"
        elif len(doys) == 2:
            doys = f"{doys[0]}_{doys[1]}"
        else:
            doys = f"{doys[0]}"

        secs = ''
        for s in config['discretization']['sec']:
            secs += f"{s}_"
        secs = secs[:-1]

        atmo_model = config['atmospheric_model']['type']
        prop_model = config['atmospheric_model']['prop_model']
        use_pert = config['atmospheric_model']['use_pert']

        project_name = sou_name+"-"+sta_name+"-"+year+"-"+doys+"-"+secs+"-"+\
            atmo_model+"-"+prop_model

        if use_pert is True:
            project_name = project_name + "-" + 'pert'

    print("\nRunning: wrap_project.sh")
    print("========================")
    proj_folder = os.path.join(config['project_info']['path'], project_name)
    ind = 1
    project_name_aux = project_name
    while isdir(proj_folder):
        print(f"[wrap_project] Project folder name \n{proj_folder}\nexists!")
        project_name_aux = project_name+f"-v{ind:d}"
        proj_folder = os.path.join(config['project_info']['path'], project_name_aux)
        print(f"[wrap_project] Trying to save as \n{proj_folder}...")
        ind += 1
    project_name = project_name_aux

    subprocess.check_output(
        [
            'bash',
            'wrap_project.sh',
            os.path.join(
                    config['project_info']['path'],
                    project_name
                    )
        ]
    )


if __name__ == '__main__':
    '''
    Default case, when calling the script from terminal
    '''
    print("\n============")
    print("   ARCADE   ")
    print("============")

    # Set working path
    os.chdir("./bin")
    print("-> Changed to working directory ./bin")

    list_of_folders = [
        "../output/proc",
        "../output/proc/arrv",
        "../output/proc/rays"
        ]

    for folder in list_of_folders:
        try:
            mkdir(folder)
            print(f"Folder {folder} created")
        except FileExistsError:
            print(f"Folder {folder} already there")

    # Check for paths that should exist
    list_of_paths = [
        "../input/secs_doys_sources_stations.txt",
        "../input/config.toml"
        ]

    for path in list_of_paths:
        try:
            open(path)
        except FileNotFoundError:
            print(f"ERROR: required file {path} doesn't exist")
            sys.exit()

    # Read TOML config file ===============================
    config = toml.load("../input/config.toml")
    # =====================================================

    # model type: 'clim', 'hybrid', or 'ncpag2s'
    atmo_type       = config['atmospheric_model']['type']
    use_pert        = config['atmospheric_model']['use_pert']
    use_rng_dep     = config['atmospheric_model']['prop_model'] == 'range_dep'

    perc_cpu        = config['launch']['perc_cpu']

    doys = config['discretization']['doys']
    doy_step = config['discretization']['doy_step']
    # Option of writing a [start, stop] and setting a doy step to create a
    # list as in [1, 2,..., 365] for whole year or long time intervals
    if doy_step > 0:
        if len(doys) == 2:
            doys = np.arange(int(float(doys[0])),
                             int(float(doys[1]))+doy_step,
                             doy_step)

    if use_rng_dep is True:
        profiles = []
        lons = None
        lats = None
        if config['atmospheric_model']['type'] == 'clim':
            lons = np.loadtxt("../output/profiles/nodes-lon.loc", dtype='float', ndmin=1)
            lats = np.loadtxt("../output/profiles/nodes-lat.loc", dtype='float', ndmin=1)
        elif config['atmospheric_model']['type'] == 'ncpag2s':
            lons = np.loadtxt("../output/profiles/ncpag2s/lons.dat")
            lats = np.loadtxt("../output/profiles/ncpag2s/lats.dat")
        profiles = get_profiles_rngdep(doys)
        # Run with multiprocessing =============================================
        num_cpu = int(mp.cpu_count()*perc_cpu)
        pool = mp.Pool(num_cpu)
        results = []
        print(f"- {len(profiles)} profiles fill be calculated with {num_cpu} cores")
        results = pool.starmap_async(
            calculate_profiles,
            [(profiles, config, atmo_type, i, use_rng_dep) for i in range(len(profiles))]).get()
        pool.close()
    else:
        profiles = []
        profiles = get_profiles(
            doys,
            pert_flag = use_pert,
            mix_flag  = config['atmospheric_model']['type'] == 'hybrid',
        )

        # Run with multiprocessing =================================================
        num_cpu = int(mp.cpu_count()*perc_cpu)
        pool = mp.Pool(num_cpu)
        results = []
        print(f"- {len(profiles)} profiles fill be calculated with {num_cpu} cores")
        results = pool.starmap_async(
            calculate_profiles,
            [(profiles, config, atmo_type, i, use_rng_dep) for i in range(len(profiles))]).get()
        pool.close()

    os.chdir("../src")

    # =========================================================
    # Run plot_results.py
    # =========================================================
    print("\nRunning: plot_arrivals.py")
    print("========================")
    import plot_arrivals
    plot_arrivals.main()

    wrap_and_save(config)
