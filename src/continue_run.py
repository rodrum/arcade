"""
Continue a run that was stopped.

The file input/secs_doys_sources_stations.txt should exist.

This file is created in the discretization process as an input for the
Fortran subroutines needed for the climatologies.

The main idea of this script is to compare `input/secs_doys_sources_stations.txt`
with the combinations in `output/azimuth_deviation_table.txt`

by RDN @ SFU
2026-07-13
"""

import sys
import multiprocessing as mp
from os.path import join, isdir
from os import listdir, chdir, mkdir
import toml
from ARCADE_main import *
import shutil

def check_file_exists(file_path):
    """Check if file exists. If not does sys.exit()"""
    try:
        f = open(file_path, 'r')
    except FileNotFoundError:
        print("[read_secs_doys_sources_stations] ERROR: file not found!")
        print(f" -> Check that <{file_path}> exists.")
        print("    Bye!")
        sys.exit()


def read_secs_doys_sources_stations(folder_path=''):
    """
    Reads `input/secs_doys_sources_stations.txt` and returns combinations.
    """

    input_path =  './input/secs_doys_sources_stations.txt'
    if folder_path != '':
        input_path = join(folder_path, 'secs_doys_sources_stations.txt')
    check_file_exists(input_path)

    combs = []
    with open(input_path, 'r') as f:
        for line in f.readlines():
            line_split = line.split()
            combs.append([
                int(float(line_split[0])),  # sec
                int(float(line_split[1])),  # doy
                int(float(line_split[2])),  # sou
                int(float(line_split[3])),  # sta
                ])
    return combs


def read_azimuth_deviation_table(folder_path=''):
    """
    Reads `output/azimuth_deviation_table.txt` and returns combinations.
    """

    input_path =  './output/azimuth_deviation_table.txt'
    if folder_path != '':
        input_path = join(folder_path, 'azimuth_deviation_table.txt')
    check_file_exists(input_path)
    combs = []
    with open(input_path, 'r') as f:
        _ = f.readline()  # header
        for line in f.readlines():
            line_split = line.split()
            combs.append([
                int(float(line_split[2])),  # Seconds
                int(float(line_split[1])),  # DOY
                int(float(line_split[3])),  # SouNum
                int(float(line_split[4])),  # StaNum
                ])
    return combs


def subtract_combs(comb_target, comb_done):
    """
    Find combinations that are in comb_target, but not in comb_done.
    """
    combs_to_do = []
    for comb_i in comb_target:
        if comb_i not in comb_done:
            combs_to_do.append(comb_i)
    return combs_to_do


def count_profiles(combs_to_do, folder_path=''):
    """
    Count profiles to see if we need to discretize again.
    """
    print('')
    input_path = './output/profiles/'
    if folder_path != '':
        input_path = join(folder_path)

    all_files = listdir(input_path)
    met_files = list(filter(lambda x: x[-4:] == '.met', all_files))

    combs_discretized = []
    for met_file in met_files:
        met_file = met_file[:-4].split('_')
        comb_i = [
            int(float(met_file[2])),  # sec
            int(float(met_file[3])),  # doy
            int(float(met_file[4])),  # sou
            int(float(met_file[5]))   # sta
            ]
        combs_discretized.append(comb_i)

    combs = []
    for comb_i in combs_to_do:
        if comb_i not in combs_discretized:
            combs.append(comb_i)

    if len(combs) == 0:
        print('[count_profiles] All combinations are discretized.')
        return combs_to_do
    else:
        print('[count_profiles] The following combinations need to be discretized:')
        for comb_i in combs:
            print(*comb_i)
        return combs



if __name__ == '__main__':
    print('---------------------------------------------------------------')
    print('                       CONTINUE RUN                            ')
    print('---------------------------------------------------------------')


    input_project = join(
        '/home/rodrigo/Desktop/arcade',
        'leftraru_runs',
        'Etna-IS26_IS42_IS48etal-2016-1_365-0_21600_43200_64800-clim-range_ind'
        )

    input_path_secs_doys_stats = 'leftraru_runs/Etna-IS26_IS42_IS48etal-2016-1_365-0_21600_43200_64800-clim-range_ind/input/'
    combs_target = read_secs_doys_sources_stations(input_path_secs_doys_stats)
    input_path_azimuth_table = 'leftraru_runs/Etna-IS26_IS42_IS48etal-2016-1_365-0_21600_43200_64800-clim-range_ind/output/'
    combs_done = read_azimuth_deviation_table(input_path_azimuth_table)
    combs_to_do = subtract_combs(combs_target, combs_done)

    print('Combinations not calculated:\n ')
    print('# Seconds DOY SouNum StaNum')
    for comb_i in combs_to_do:
        print(*comb_i)

    combs_to_discretize = count_profiles(
        combs_to_do,
        'leftraru_runs/Etna-IS26_IS42_IS48etal-2016-1_365-0_21600_43200_64800-clim-range_ind/output/profiles'
        )

    config = toml.load(join(input_project, 'input', 'config.toml'))

    sta_net = None
    try:
        sta_net = config['discretization']['sta_net']
    except KeyError:
        pass

    sta_locs = []
    if not sta_net:
        sta_locs = config['discretization']['sta_pos']
    elif sta_net == 'IMS':
        ims_sta_locs = pd.read_csv(join(input_project, 'input', 'ims_stations.csv'))
        sta_locs = []
        for sta_name in config['discretization']['sta_name']:
            sta_locs.append([
                ims_sta_locs[sta_name]['lat (deg)'].values[0],
                ims_sta_locs[sta_name]['lon (deg)'].values[0]
                ])

    sou_locs = config['discretization']['sou_pos']
    sta_names = config['discretization']['sta_name']

    pert_flag = config['atmospheric_model']['use_pert']
    mix_flag  = config['atmospheric_model']['type'] == 'hybrid'
    end_str = '.met'
    if pert_flag is True and mix_flag is False:
        end_str = '_pert.met'
    elif pert_flag is False and mix_flag is True:
        end_str = '_mix.met'
    elif pert_flag is True and mix_flag is True:
        end_str = '_pert.met'

    # Set working path ---------------------------------------------------------
    chdir("./bin")
    print("-> Changed to working directory ./bin")

    # Check that output folder and its subfolders exists or create them
    if not isdir("../output"):
        mkdir("../output")
    out_path = "../output/nodes"
    if not isdir(out_path):
        mkdir(out_path)
    out_path_prof = "../output/profiles/"
    if not isdir(out_path_prof):
        mkdir(out_path_prof)

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
    # ------------------------------------------------------------------------//

    profiles = []
    if len(combs_to_discretize) == len(combs_to_do):
        for sec, doy, sou_num, sta_num in combs_to_do:
            sta_lat = sta_locs[sta_num-1][0]
            sta_lon = sta_locs[sta_num-1][1]
            sou_lat = sou_locs[sou_num-1][0]
            sou_lon = sou_locs[sou_num-1][1]
            sta_nam = sta_names[sta_num-1]
            prof_name = f"prof_{sec:05d}_{doy:03d}_{sou_num:05d}_{sta_num:04d}"+\
                f"{end_str}"
            profiles.append((
                sta_lat, sta_lon, sou_lat, sou_lon, sta_nam, prof_name
                ))
            for i in range(2):
                prof_i = f"{i+1:d}_{prof_name}"
                shutil.copyfile(
                    join(input_project, 'output', 'profiles', prof_i),
                    join('..', 'output', 'profiles', prof_i)
                    )
    else:
        print("Not implemented yet, please finish")
        sys.exit()


    use_rng_dep = config['atmospheric_model']['prop_model'] == 'range_dep'
    atmo_type = config['atmospheric_model']['type']
    perc_cpu = config['launch']['perc_cpu']
    if use_rng_dep is True:
        lons = None
        lats = None
        if config['atmospheric_model']['type'] == 'clim':
            lons = np.loadtxt(
                join(input_project, "output/profiles/nodes-lon.loc"),
                dtype='float',
                ndmin=1
                )
            lats = np.loadtxt(
                join(input_project,"output/profiles/nodes-lat.loc"),
                dtype='float',
                ndmin=1
                )
        elif config['atmospheric_model']['type'] == 'ncpag2s':
            lons = np.loadtxt(join(input_project,"output/profiles/ncpag2s/lons.dat"))
            lats = np.loadtxt(join(input_project,"output/profiles/ncpag2s/lats.dat"))

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
        # Run with multiprocessing =================================================
        num_cpu = int(mp.cpu_count()*perc_cpu)
        pool = mp.Pool(num_cpu)
        results = []
        print(f"- {len(profiles)} profiles fill be calculated with {num_cpu} cores")
        results = pool.starmap_async(
            calculate_profiles,
            [(profiles, config, atmo_type, i, use_rng_dep) for i in range(len(profiles))]).get()
        pool.close()

    chdir("../src")

    wrap_and_save(config, '/home/rodrigo/Desktop/arcade/test_runs/')
