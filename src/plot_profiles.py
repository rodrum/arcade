"""
Code to plot profiles automatically from '.met' files and output from 
'discretize.py'
"""

import numpy as np
import matplotlib.pyplot as plt
from os.path import join
from os import mkdir, stat
import cartopy.crs as ccrs
from geographiclib.geodesic import Geodesic
import matplotlib.image as mpimg
import toml
from obspy.geodetics.base import gps2dist_azimuth
from matplotlib import colors
from itertools import product

config = toml.load("./input/config.toml")

plot_arrivals = config['launch']['plot_arrivals']
plot_profiles = config['launch']['plot_profiles']

# NOTE: added this to skip case
is_ncpag2s_model = config['atmospheric_model']['type'] == 'ncpag2s'

geod = Geodesic.WGS84

path_input      = "./input"
path_profiles   = "./output/profiles"
path_figures    = "./output/figures"

secs        = config['discretization']['sec']
doys        = config['discretization']['doys']
sources     = config['discretization']['sou_pos']
stations    = config['discretization']['sta_pos']

all_comb = [combi for combi in
            product(secs, doys, range(len(sources)), range(len(stations)))]

try:
    mkdir(path_figures)
    print(f"Folder {path_figures} created")
except FileExistsError:
    print(f"Folder {path_figures} already there")

if plot_profiles is True and is_ncpag2s_model is False:
    profile_num = 0
    for sec, doy, isou, ista in all_comb:
        sou_lat, sou_lon = sources[isou][0], sources[isou][1]
        sta_lat, sta_lon = stations[ista][0], stations[ista][1]
        # world map
        n = 20
        l = geod.InverseLine(sou_lat, sou_lon, sta_lat, sta_lon)
        # along arc discretization
        arc = np.linspace(0, l.s13, n)
        arc_lons = []
        arc_lats = []
        for s in arc:
            g = l.Position(s, Geodesic.STANDARD | Geodesic.LONG_UNROLL)
            arc_lons.append(g['lon2'])
            arc_lats.append(g['lat2'])
        cen_pos = l.Position(
            arc[int(n/2)], Geodesic.STANDARD | Geodesic.LONG_UNROLL
            )
        central_latitude = cen_pos['lat2']
        central_longitude = cen_pos['lon2']

        plt.figure(figsize=(6,3))
        ax0 = plt.subplot(121, 
            projection=ccrs.Orthographic(
                central_longitude=central_longitude,
                central_latitude=central_latitude)
            )

        ax0.stock_img()
        ax0.plot(arc_lons, arc_lats,
                    color='blue', linewidth=2, linestyle='--',
                    transform=ccrs.Geodetic(),
                    )

        ax0.plot(sou_lon, sou_lat, marker='^', color='r', 
                    transform=ccrs.Geodetic())
        ax0.plot(sta_lon, sta_lat, marker='v', color='k', 
                    transform=ccrs.Geodetic())


        ax1 = plt.subplot(122)
        ax1.text(0.5, 0.9, 
            f"Source : ({sou_lat:.1f}$^\circ$,{sou_lon:.1f}$^\circ$)",
            fontsize=12, horizontalalignment='left', 
            verticalalignment='center', transform=ax1.transAxes
            )
        ax1.plot(0.3, 0.9, 
            marker='^', markersize=10, color='r', 
            transform=ax1.transAxes
            )
        ax1.text(0.5, 0.7, 
            f"Station: ({sta_lat:.1f}$^\circ$,{sta_lon:.1f}$^\circ$)",
            fontsize=12, horizontalalignment='left',
            verticalalignment='center', transform=ax1.transAxes
            )
        ax1.plot(0.3, 0.7, 
            marker='v', markersize=10, color='k', 
            transform=ax1.transAxes
            )
        ax1.text(0.5, 0.5, 
            "Geodesic arc", fontsize=12, horizontalalignment='left',
            verticalalignment='center', 
            transform=ax1.transAxes
            )
        ax1.plot([0.2,0.4], [0.5,0.5], color='b', linewidth=2, 
            linestyle='--', 
            transform=ax1.transAxes
            )
        ax1.text(0.5, 0.4, 
            f"-  Distance: {l.s13/1000:.1f} km", fontsize=10, 
            horizontalalignment='left', verticalalignment='center', 
            transform=ax1.transAxes
            )
        ax1.text(0.5, 0.3, 
            f"-  Azimuth : {l.azi1:.1f}$^\circ$", fontsize=10, 
            horizontalalignment='left', verticalalignment='center', 
            transform=ax1.transAxes
            )

        ax1.axis('off')

        # Save fig
        namefig = f"glob_{sec:05d}_{doy:03d}_{isou+1:05d}_{ista+1:04d}.png"
        out_fig_path = join(path_figures, namefig)
        plt.savefig(out_fig_path, dpi=300, bbox_inches='tight')

        plt.close() 
        profile_num += 1
        print(f"   > Figure {out_fig_path} saved")

    if config['atmospheric_model']['prop_model'] == 'range_ind':
        end_str = '.met'
        if config['atmospheric_model']['type'] == 'hybrid':
            end_str =  '_mix.met'
        profile_num = 0
        for sec, doy, isou, ista in all_comb:
            sou_lat, sou_lon = sources[isou][0], sources[isou][1]
            sta_lat, sta_lon = stations[ista][0], stations[ista][1]
            prof_name = f"1_prof_{sec:05d}_{doy:03d}_{isou+1:05d}_{ista+1:04d}{end_str}"
            prof_name_path = join(path_profiles, prof_name)
            print(f"Processing {prof_name_path}")
            prof = np.loadtxt(prof_name_path)
            h = prof[:, 0]
            temp = prof[:, 1]
            zonw = prof[:, 2]
            merw = prof[:, 3]
            dens = prof[:, 4]
            pres = prof[:, 5]
            # Perturbed profiles
            pert_prof_name = join(
                path_profiles,
                f"1_prof_{sec:05d}_{doy:03d}_{isou+1:05d}_{ista+1:04d}_pert.met"
            )
            if stat(pert_prof_name).st_size > 0:
                pert_prof = np.loadtxt(pert_prof_name)
                zonw_pert = pert_prof[:, 2]
                merw_pert = pert_prof[:, 3]
                # Get azimuth source -> station
                l = geod.InverseLine(sou_lat, sou_lon, sta_lat, sta_lon)
                a12 = np.deg2rad(l.azi1)
                # Effective winds
                kx = np.sin(a12)
                ky = np.cos(a12)
                # Along winds
                c_along = zonw*kx + merw*ky
                c_pert_along = zonw_pert*kx + merw_pert*ky
                # Across winds
                #   rotate profile direction vector +90 degrees
                kx = -ky
                ky = kx
                c_pert_across = zonw_pert*kx + merw_pert*ky
                c_across = zonw*kx + merw*ky

                #=== create fig
                fig, axs = plt.subplots(
                num=profile_num,
                nrows = 3,
                ncols = 5,
                figsize = (8, 10),
                gridspec_kw = {"width_ratios" : [1, 1, 1, 1, 1]},
                )

                # Merge upper five cols into one plot
                gs = axs[0, 0].get_gridspec()
                for ax in axs[0, 0:5]:
                    ax.remove()
                ax_globe = fig.add_subplot(gs[0, 0:5])
                ax0 = axs[1,0]  # Temp [K]
                ax1 = axs[1,1]  # Zonal w. [m/s]
                ax2 = axs[1,2]  # Mer w. [m/s]
                ax3 = axs[1,3]  # Dens. [g/cm^3]
                ax4 = axs[1,4]  # Pres. [mbar] (Pa*10^-2)

                for ax in axs[2,:]:
                    ax.remove()
                ax5 = fig.add_subplot(gs[2, 0:2])  # Along winds [m/s]
                ax6 = fig.add_subplot(gs[2, 3:5])  # Across winds [m/s]

                # Make ticks go inwards, and background gray for each plot
                ax0.tick_params(direction='in')
                ax0.set_facecolor('xkcd:black')
                ax1.tick_params(direction='in')
                ax1.set_facecolor('xkcd:black')
                ax2.tick_params(direction='in')
                ax2.set_facecolor('xkcd:black')
                ax3.tick_params(direction='in')
                ax3.set_facecolor('xkcd:black')
                ax4.tick_params(direction='in')
                ax4.set_facecolor('xkcd:black')
                ax5.tick_params(direction='in')
                ax5.set_facecolor('xkcd:black')
                ax6.tick_params(direction='in')
                ax6.set_facecolor('xkcd:black')

                ax_globe.axis('off')

                # take y-axis out for last 4 panels
                ax1.set_yticklabels([])
                ax2.set_yticklabels([])
                ax3.set_yticklabels([])
                ax4.set_yticklabels([])
                ax6.set_yticklabels([])

                # set dens and press x-scales as log
                ax3.set_xscale('log')
                ax4.set_xscale('log')

                # add axis labels
                ax0.set_xlabel('Temp. [K]')
                ax1.set_xlabel('Zonal w. [m/s]')
                ax2.set_xlabel('Merid. w. [m/s]')
                ax3.set_xlabel('Dens. [g/cm^3]')
                ax4.set_xlabel('Pres. [mbar]')
                ax5.set_xlabel('Along winds [m/s]')
                ax6.set_xlabel('Across winds [m/s]')

                # add yaxis label
                ax0.set_ylabel("Height [km]")
                ax5.set_ylabel("Height [km]")


                # title
                title = "Profile "\
                    +f"prof_{sec:05d}_{doy:03d}_{isou+1:05d}_{ista+1:04d}{end_str}"\
                    +f" - DOY {doy:03d} "
                plt.suptitle(title, fontsize=14)

                # grid lines
                ax0.grid()
                ax1.grid()
                ax2.grid()
                ax3.grid()
                ax4.grid()
                ax5.grid()
                ax6.grid()

                # plot
                ax0.plot(temp, h)
                ax1.plot(zonw, h)
                ax2.plot(merw, h)
                ax3.plot(dens, h)
                ax4.plot(pres, h)
                ax5.plot(c_along, h, label="HWM14")
                ax5.plot(c_pert_along, h, linestyle='--', label="HWM14 pert.")
                ax6.plot(c_across, h, label="HWM14")
                ax6.plot(c_pert_across, h, linestyle='--', label="HWM14 pert.")

                # add labels for along and across winds
                ax5.legend(loc=1)
                ax6.legend(loc=1)

                # import
                namefig = f"glob_{sec:05d}_{doy:03d}_{isou+1:05d}_{ista+1:04d}.png"
                in_fig_path = join(path_figures, namefig)
                img = mpimg.imread(in_fig_path)
                imgplot = ax_globe.imshow(img)
                print(f"   > Figure {in_fig_path} loaded") 

                # save
                namefig = f"prof_{sec:05d}_{doy:03d}_{isou+1:05d}_{ista+1:04d}.png"
                out_fig_path = join(path_figures, namefig)
                plt.savefig(out_fig_path, dpi=300, bbox_inches='tight')

                plt.close(fig) 
                profile_num += 1
                print(f"   > Figure {out_fig_path} saved")
            else:
                #=== create fig
                fig, axs = plt.subplots(
                    num=profile_num,
                    nrows = 2,
                    ncols = 5,
                    figsize = (8, 8),
                    gridspec_kw = {"width_ratios" : [1, 1, 1, 1, 1]},
                )
                # Merge upper five cols into one plot
                gs = axs[0, 0].get_gridspec()
                for ax in axs[0, 0:5]:
                    ax.remove()
                ax_globe = fig.add_subplot(gs[0, 0:5])
                ax0 = axs[1,0]  # Temp [K]
                ax1 = axs[1,1]  # Zonal w. [m/s]
                ax2 = axs[1,2]  # Mer w. [m/s]
                ax3 = axs[1,3]  # Dens. [g/cm^3]
                ax4 = axs[1,4]  # Pres. prof_name[mbar] (Pa*10^-2)

                # Make ticks go inwards, and background gray for each plot
                ax0.tick_params(direction='in')
                ax0.set_facecolor('xkcd:black')
                ax1.tick_params(direction='in')
                ax1.set_facecolor('xkcd:black')
                ax2.tick_params(direction='in')
                ax2.set_facecolor('xkcd:black')
                ax3.tick_params(direction='in')
                ax3.set_facecolor('xkcd:black')
                ax4.tick_params(direction='in')
                ax4.set_facecolor('xkcd:black')

                ax_globe.axis('off')

                # take y-axis out for last 4 panels
                ax1.set_yticklabels([])
                ax2.set_yticklabels([])
                ax3.set_yticklabels([])
                ax4.set_yticklabels([])

                # set dens and press x-scales as log
                ax3.set_xscale('log')
                ax4.set_xscale('log')

                # add axis labels
                ax0.set_xlabel('Temp. [K]')
                ax1.set_xlabel('Zonal w. [m/s]')
                ax2.set_xlabel('Merid. w. [m/s]')
                ax3.set_xlabel('Dens. [g/cm^3]')
                ax4.set_xlabel('Pres. [mbar]')

                # add yaxis label
                ax0.set_ylabel("Height [km]")

                # title
                title = "Profile "\
                    +f"prof_{sec:05d}_{doy:03d}_{isou+1:05d}_{ista+1:04d}.met"\
                    +f" - DOY {doy:03d} "
                plt.suptitle(title, fontsize=14)

                # grid lines
                ax0.grid()
                ax1.grid()
                ax2.grid()
                ax3.grid()
                ax4.grid()

                # plot
                ax0.plot(temp, h)
                ax1.plot(zonw, h)
                ax2.plot(merw, h)
                ax3.plot(dens, h)
                ax4.plot(pres, h)

                # import
                namefig = f"glob_{sec:05d}_{doy:03d}_{isou+1:05d}_{ista+1:04d}.png"
                in_fig_path = join(path_figures, namefig)
                img = mpimg.imread(in_fig_path)
                imgplot = ax_globe.imshow(img)
                print(f"   > Figure {in_fig_path} loaded") 

                # save
                namefig = f"prof_{sec:05d}_{doy:03d}_{isou+1:05d}_{ista+1:04d}.png"
                out_fig_path = join(path_figures, namefig)
                plt.savefig(out_fig_path, dpi=300, bbox_inches='tight')

                plt.close(fig) 
                profile_num += 1
                print(f"   > Figure {out_fig_path} saved")
    elif config['atmospheric_model']['prop_model'] == 'range_dep':
        print("\nNote: this plot type has been implemented only for a case of one ")
        print("      station, one source, and one day of the year. Trying more")
        print("      than one of each could create unexpected output plots.\n")

        def get_wind_values(x_points, y_points, nsec, ndoy, nsou, nsta, num_h=50):
            zonw_values = np.zeros(shape=(x_points.shape[0], y_points.shape[0]))
            merw_values = np.zeros(shape=(x_points.shape[0], y_points.shape[0]))
            prof_num = 0
            for xi in range(x_points.shape[0]):
                for yi in range(y_points.shape[0]):
                    prof_name = f"1_prof_{nsec:05d}_{ndoy:03d}_{nsou+1:05d}"\
                                f"_{nsta+1:04d}_{prof_num}.met"
                    prof_vals = np.loadtxt(f"./output/profiles/{prof_name}")
                    zonw_values[xi, yi] = prof_vals[num_h, 2]
                    merw_values[xi, yi] = prof_vals[num_h, 3] 
                    prof_num += 1
            return zonw_values, merw_values

        x_points = np.loadtxt('./output/profiles/nodes-lon.loc')
        y_points = np.loadtxt('./output/profiles/nodes-lat.loc')
        grid_x, grid_y = np.meshgrid(x_points, y_points)

        profile_num = 0
        nr, nc = 2, 3
        alts = [25, 50, 75, 100, 125, 150]
        for nsec, ndoy, nsou, nsta in all_comb:
            sou_lat, sou_lon = sources[nsou][0], sources[nsou][1]
            sta_lat, sta_lon = stations[nsta][0], stations[nsta][1]
            # Azimuth from source to station
            _, a12, _ = gps2dist_azimuth(sou_lat, sou_lon, sta_lat, sta_lon)
            print(f"a12={a12:.2f} deg.")
            a12 = np.deg2rad(a12) # pass to radians
            # Along winds
            kx = np.sin(a12)
            ky = np.cos(a12)

            # find max value
            max_wind = -1.0e10
            zonw_vals = []
            merw_vals = []
            alti = 0 
            for i in range(nr):
                for j in range(nc):
                    zonw, merw = get_wind_values(
                        x_points, y_points, nsec, ndoy, nsou, nsta, alts[alti]*2
                        )
                    zonw_vals.append(zonw)
                    merw_vals.append(merw)
                    wind = np.sqrt(zonw**2 + merw**2)
                    wind_max = np.amax(wind)
                    if wind_max > max_wind:
                        max_wind = wind_max
                    alti+=1
            print(f"Max. wind: {max_wind:.2f} m/s")

            norm=colors.Normalize(vmin=-1, vmax=1)

            along_coef = []
            k = 0
            for i in range(nr):
                for j in range(nc):
                    zonw, merw = zonw_vals[k], merw_vals[k]
                    along_wind = zonw*kx+merw*ky
                    along_wind /= np.abs(along_wind)
                    along_coef.append(along_wind)
                    k+=1

            #=== create fig
            fig, axs = plt.subplots(
                num=profile_num,
                nrows = 3,
                ncols = 4,
                figsize = (12,12),
                gridspec_kw = {"width_ratios" : [1, 1, 1, 0.05]},
                )
                # Merge upper five cols into one plot
            gs = axs[0, 0].get_gridspec()
            for ax in axs[0, 0:4]:
                ax.remove()
            ax_globe = fig.add_subplot(gs[0, 0:4])
            for ax in axs[1:3, 3]:
                ax.remove()
            ax_cbar = fig.add_subplot(gs[1:3, 3])
            ax00 = axs[1,0]   
            ax01 = axs[1,1]   
            ax02 = axs[1,2]   
            ax10 = axs[2,0]   
            ax11 = axs[2,1]   
            ax12 = axs[2,2]   

            # Make ticks go inwards, and background gray for each plot
            ax00.tick_params(direction='in')
            ax00.set_facecolor('xkcd:gray')
            ax00.grid()
            ax01.tick_params(direction='in')
            ax01.set_facecolor('xkcd:gray')
            ax01.grid()
            ax02.tick_params(direction='in')
            ax02.set_facecolor('xkcd:gray')
            ax02.grid()
            ax10.tick_params(direction='in')
            ax10.set_facecolor('xkcd:gray')
            ax10.grid()
            ax11.tick_params(direction='in')
            ax11.set_facecolor('xkcd:gray')
            ax11.grid()
            ax12.tick_params(direction='in')
            ax12.set_facecolor('xkcd:gray')
            ax12.grid()

            ax_globe.axis('off')

            # take y-axis out for last 4 panels
            ax01.set_yticklabels([])
            ax02.set_yticklabels([])
            ax11.set_yticklabels([])
            ax12.set_yticklabels([])

            # take x-axis out for last 4 panels
            ax00.set_xticklabels([])
            ax01.set_xticklabels([])
            ax02.set_xticklabels([])

            # add axis labels
            ax10.set_xlabel('Longitude [$^{\circ}$]')
            ax11.set_xlabel('Longitude [$^{\circ}$]')
            ax12.set_xlabel('Longitude [$^{\circ}$]')
            ax00.set_ylabel('Latitude [$^{\circ}$]')
            ax10.set_ylabel('Latitude [$^{\circ}$]')

            # title
            title = f"Profile prof_{nsec:05d}_{ndoy:03d}_{nsou+1:05d}_{nsta+1:04d}.met"\
                    +f" - DOY {ndoy:03d} "
            plt.suptitle(
                title,
                fontsize=14
                )

            # grid lines
            ax00.grid()
            ax01.grid()
            ax02.grid()
            ax10.grid()
            ax11.grid()
            ax12.grid()

            alti = 0
            k = 0
            for i in range(nr):
                i = i+1
                for j in range(nc):
                    zonw, merw = zonw_vals[k], merw_vals[k]
                    # Generate data with a range that varies from one plot to the next.
                    q = axs[i,j].quiver(grid_x, grid_y, zonw, merw, along_coef[k],
                                pivot='tail',
                                units='xy',
                                angles='xy',
                                scale_units='xy',
                                scale=max_wind/1,
                                cmap="viridis",
                                norm=norm
                                )
                    axs[i,j].plot(sou_lon, sou_lat, '^r')
                    axs[i,j].plot(sta_lon, sta_lat, 'vk')
                    #axs[i,j].plot([sou_lon, sta_lon], [sou_lat, sta_lat], '--b')
                    head_length = 0.4
                    vec_magnitude = np.sqrt((sta_lon-sou_lon)**2+(sta_lat-sou_lat)**2)
                    dx = (sta_lon-sou_lon)/vec_magnitude
                    dy = (sta_lat-sou_lat)/vec_magnitude
                    vec_ab_mag = vec_magnitude-head_length
                    axs[i,j].arrow(
                        sou_lon, sou_lat, 
                        vec_ab_mag*dx, vec_ab_mag*dy, 
                        head_width=0.4, head_length=head_length, 
                        fc='white', ec='black'
                        )
                    axs[i,j].set_title(f"Altitude: {alts[alti]:.2f} km")
                    axs[i,j].grid()
                    if i==1 and j==2:
                        axs[1,1].quiverkey(q, X=1.2, Y=1.1, U=100,
                                    label="Scale: 100 m/s", labelpos='N')
                    alti+=1
                    k+=1

            SMALL_SIZE = 12 
            MEDIUM_SIZE = 14
            BIGGER_SIZE = 16

            plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
            plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
            plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
            plt.rc('xtick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
            plt.rc('ytick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
            plt.rc('legend', fontsize=MEDIUM_SIZE)    # legend fontsize
            plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

            fig.colorbar(
                q, cax=ax_cbar, orientation='vertical', 
                label=r'Along coeff. ($\cos(\theta)$)'
                )

            # import
            namefig = f"glob_{nsec:05d}_{ndoy:03d}_{nsou+1:05d}_{nsta+1:04d}.png"
            in_fig_path = join(path_figures, namefig)
            img = mpimg.imread(in_fig_path)
            imgplot = ax_globe.imshow(img)
            print(f"   > Figure {in_fig_path} loaded") 

            # save
            namefig = f"prof_{nsec:05d}_{ndoy:03d}_{nsou+1:05d}_{nsta+1:04d}.png"
            out_fig_path = join(path_figures, namefig)

            plt.savefig(
                out_fig_path,
                dpi=300, bbox_inches='tight'
                )

            profile_num += 1
else:
    print("NOTE: Skipping profile plots...")
