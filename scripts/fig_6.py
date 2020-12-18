#!/usr/bin/env python3

# Import utils- and plot-specific modules
from modules_plot import *

#====================================================================
def plot_atmosphere( output_dir, sub_dirs, output_times ):

    width   = 6.00 #* 3.0/2.0
    height  = 12.0
    fig_o   = FigureData( 2, 1, width, height, '../figures/fig_6', units='kyr' )
    fig_o.fig.subplots_adjust(wspace=0.15, hspace=0.0)

    ax0 = fig_o.ax[0]
    ax1 = fig_o.ax[1]
    handle_l = [] # handles for legend

    ymax_atm_pressure   = 0
    ymin_atm_pressure   = 1000
    ymax_atm_z          = 0
    ymin_atm_z          = 0
    ymax_sp_flux        = 0
    xmin                = 100
    xmax                = 3000

    fs_legend   = 11
    fs_label    = 11

    # Initiate legends
    legend_ax0_1_handles = []
    legend_ax0_2_handles = []

    for subdir in sub_dirs:

        data_dir = output_dir+"/"+subdir

        # Get data times for subdir
        data_times = get_all_output_pkl_times(data_dir)

        # Find solidus and liquidus lines
        time        = data_times[0]                    
        myjson_o    = MyJSON( data_dir+'/{}.json'.format(time) )
        xx_pres     = myjson_o.get_dict_values(['data','pressure_b'])*1.0E-9
        xx_pres_s   = myjson_o.get_dict_values(['data','pressure_s'])*1.0E-9
        xx_radius   = myjson_o.get_dict_values(['data','radius_b'])*1.0E-3
        xx_depth    = xx_radius[0] - xx_radius
        xx_radius_s = myjson_o.get_dict_values(['data','radius_s'])*1.0E-3
        xx_depth_s  = xx_radius_s[0] - xx_radius_s
        r_planet    = np.max(xx_radius*1e3) # m

        # Find planet mass
        core_mass   = myjson_o.get_dict_values(['atmosphere','mass_core'])
        mantle_mass = myjson_o.get_dict_values(['atmosphere','mass_mantle'])
        planet_mass = core_mass + mantle_mass

        color_idx   = 5

        for nn, output_time in enumerate( output_times ):

            if nn == 0:
                lw = 2.0
                ls = "-"
            if nn == 1:
                lw = 2.0
                ls = "--"
            if nn == 3:
                lw = 2.0
                ls = ":"

            print(subdir, output_time)

            # Find data_time closest to output_time wanted
            time, time_idx = find_nearest(data_times, output_time)

            label = vol_latex[subdir]#+", "+latex_float(time)+" yr"

            atm_file = data_dir+"/"+str(int(time))+"_atm.pkl"

            # Read pickle file
            atm_file_stream = open(atm_file,'rb')
            atm = pkl.load(atm_file_stream)
            atm_file_stream.close()

            # color = fig_o.get_color( nn )
            color   = vol_colors[subdir][color_idx]
          
            # Atmosphere T-Z
            z_profile = AtmosphericHeight(atm, planet_mass, r_planet) # m
            z_profile = z_profile*1e-3 # km

            # Connect interior and atmosphere lines
            atm.tmp[-1] = atm.ts

            # Find pressure level closest to 1 Pa
            prs_cut, prs_cut_idx = find_nearest(atm.p, 1)
            if np.amax(atm.p) < 1: prs_cut_idx = 0

            # Plot height coordinates
            l1, = ax0.plot( atm.tmp[prs_cut_idx:], z_profile[prs_cut_idx:], lw=lw, color=color, label=label, ls=ls)

            # Fill color and P_surf legends
            if nn == 0:
                legend_ax0_1_handles.append(l1)

            if subdir == sub_dirs[0]:
                if output_time == 5e+2: 
                    time_label = "100"
                elif output_time == 1e+7: 
                    time_label = "10"
                else: 
                    time_label = latex_float(output_time)#+" yr"
                l2, = ax0.plot( 0, 0, lw=lw, color=qgray_dark, label=time_label, ls=ls)
                legend_ax0_2_handles.append(l2)

            ## INTERIOR

            myjson_o = MyJSON( data_dir+"/"+'{}.json'.format(time) )

            # T-P mantle
            mantle_prs  = myjson_o.get_dict_values(['data','pressure_b'])*1.0E-9 # GPa
            mantle_tmp  = myjson_o.get_dict_values(['data','temp_b'])

            ax1.plot( mantle_tmp, mantle_prs, ls=ls, lw=lw, color=color )

            # Get rheological front parameters
            rf_idx = int(myjson_o.get_dict_values(['rheological_front_dynamic', 'mesh_index']))
            rf_depth = round(myjson_o.get_dict_values(['rheological_front_dynamic','depth']))/1e+3 # km
            rf_prs = round(myjson_o.get_dict_values(['rheological_front_dynamic','pressure']))/1e+9 # GPa

            # Reset y-axis boundaries
            if np.min(atm.p) > 0:
                ymax_atm_pressure = np.max([ymax_atm_pressure, np.max(atm.p/1e5)])
                ymin_atm_pressure = np.min([ymin_atm_pressure, np.min(atm.p/1e5)])
                ymax_atm_z = np.max([ymax_atm_z, np.max(z_profile)])
                ymin_atm_z = np.min([ymin_atm_z, np.min(z_profile)])

            # Reset x-axis boundaries
            xmin = np.min([ xmin, np.min(atm.tmp) ])
            xmax = np.max([ xmax, np.max(mantle_tmp) ])

    # Print liquidus and solidus
    xx_pres = myjson_o.get_dict_values(['data','pressure_b'])*1.0E-9
    yy_liqt = myjson_o.get_dict_values(['data','liquidus_temp_b'])
    yy_solt = myjson_o.get_dict_values(['data','solidus_temp_b'])
    ax1.fill_betweenx( xx_pres, yy_liqt, yy_solt, facecolor=qmagenta_light, alpha=0.3, linewidth=0 )

    ### Figure settings
    
    # xmin    = 0
    xticks  = [xmin, 500, 1000, 1500, 2000, 2500, 3000, 3500, 4000, xmax]
    title_xcoord = -0.09
    title_ycoord = 0.5

    #####  T-Z
    ax0.set_ylabel("Atmosphere height, $z_\mathrm{atm}$ (km)", fontsize=fs_label)
    ax0.yaxis.set_label_coords(title_xcoord,title_ycoord)
    ax0.set_xlim( xmin, xmax )
    ax0.set_xticks(xticks)
    ax0.set_xticklabels([])
    ax0.set_yscale("log")
    ax0.set_yscale("symlog", linthreshy=50)
    ax0.set_ylim( 0, ymax_atm_z )
    yticks = [ 10, 20, 30, 40, 50, 100, 200, 300, 500, 1000, ymax_atm_z]
    ax0.set_yticks( yticks )
    ax0.set_yticklabels( [ str(int(i)) for i in yticks ] )
    
    ax1.set_ylim([0, int(mantle_prs[-1])])
    ax1.set_yticks([])
    ax1.set_yticklabels([])
    ax1.set_xlim( xmin, xmax )
    ax1.set_xticks(xticks)
    ax1.set_xlabel( r'Temperature, $T$ (K)', fontsize=fs_label )
    ax1.invert_yaxis()


    ax1b = ax1.twinx()
    ax1b.plot( mantle_tmp, xx_depth, alpha=0.0)
    ax1b.set_xlim( right=xmax, left=xmin )
    ax1b.set_xticks(xticks)
    ax1b.set_ylim(top=xx_depth[-1], bottom=xx_depth[0])
    yticks = [ 0, 250, 500, 1000, 1500, 2000, 2500, xx_depth[-1] ]
    ax1b.set_yticks( yticks )
    ax1b.set_yticklabels( [ str(int(i)) for i in yticks ] )
    ax1b.invert_yaxis()
    ax1b.yaxis.tick_left()
    ax1.set_ylabel( 'Mantle depth, $d_\mathrm{mantle}$ (km)', fontsize=fs_label )
    ax1.yaxis.set_label_coords(title_xcoord,title_ycoord)
     
    ax1.text(0.55, 0.48, '100%\nsolid', color=qblue_dark, rotation=0, ha="center", va="center", fontsize=fs_label, transform=ax1.transAxes, bbox=dict(fc='white', ec="white", alpha=0.01, pad=0.1, boxstyle='round'))
    
    ax1.text(0.72, 0.48, 'Mixed phase', color=qmagenta_dark, rotation=-62, ha="center", va="center", fontsize=fs_label-1.5, transform=ax1.transAxes, bbox=dict(fc='white', ec="white", alpha=0.01, pad=0.1, boxstyle='round'))
    
    ax1.text(0.95, 0.48, '100%\nmelt', color=qred_dark, rotation=0, ha="center", va="center", fontsize=fs_label, transform=ax1.transAxes, bbox=dict(fc='white', ec="white", alpha=0.01, pad=0.1, boxstyle='round'))

    ax0.arrow(3000, 30, 0.0, 18, head_width=50, head_length=1.5, fc=qgray_light, ec=qgray_light, lw=1.0, alpha=0.7) # , transform=ax0.transAxes
    ax0.arrow(3000, 52, 0.0, 80, head_width=50, head_length=11, fc=qgray_light, ec=qgray_light, lw=1.0, alpha=0.7) # , transform=ax0.transAxes
    ax0.text(3050, 40, 'linear\nscale', color=qgray_light, rotation=0, ha="left", va="center", fontsize=fs_label-1, alpha=0.99, bbox=dict(fc='white', ec="white", alpha=0.01, pad=0.1, boxstyle='round')) # , transform=ax0.transAxes
    ax0.text(3050, 80, 'log\nscale', color=qgray_light, rotation=0, ha="left", va="center", fontsize=fs_label-1, alpha=0.99, bbox=dict(fc='white', ec="white", alpha=0.01, pad=0.1, boxstyle='round')) # , transform=ax0.transAxes

    try:
        sns.set_style("ticks")
        sns.despine()
    except:
        print("No seaborn.")

    ax0.text(0.02, 0.985, 'A', color="k", rotation=0, ha="left", va="top", fontsize=fs_label+5, transform=ax0.transAxes, bbox=dict(fc='white', ec="white", alpha=0.1, pad=0.2, boxstyle='round'))
    ax1.text(0.02, 0.985, 'B', color="k", rotation=0, ha="left", va="top", fontsize=fs_label+5, transform=ax1.transAxes, bbox=dict(fc='white', ec="white", alpha=0.1, pad=0.2, boxstyle='round'))

    # Legend(s)
    legend_ax0_1 = ax0.legend(handles=legend_ax0_1_handles, loc=1, ncol=1, fontsize=fs_legend, framealpha=0.3, title="Volatiles")
    ax0.add_artist(legend_ax0_1)
    legend_ax0_2 = ax1.legend(handles=legend_ax0_2_handles, loc=3, ncol=1, fontsize=fs_legend, framealpha=0.3, title="Time (yr)")
    ax1.add_artist(legend_ax0_2)

    fig_o.savefig(1)
    plt.close()

#====================================================================
def main():

    # Optional command line arguments for running from the terminal
    # Usage: $ python plot_atmosphere.py -t 0,718259
    parser = argparse.ArgumentParser(description='COUPLER plotting script')
    parser.add_argument('-odir', '--output_dir', type=str, help='Full path to output directory');
    parser.add_argument('-t', '--times', type=str, help='Comma-separated (no spaces) list of times');
    args = parser.parse_args()

    # Define output directory for plots
    if args.output_dir:
        output_dir = args.output_dir
        print("Output directory:", output_dir)
        
    else:
        output_dir = os.getcwd()
        print("Output directory:", output_dir)

    # Define specific one
    output_dir  = "../data/int_atm/coupler_runs/"
    
    # Define volatiles considered
    sub_dirs    = [ "H2", "CH4", "H2O", "CO2", "N2", "CO", "O2" ]

    # Times at which to plot
    output_times = [ 1e+2, 1e+6 ]

    print("Host directory:", output_dir)

    # Plot fixed set from above
    plot_atmosphere( output_dir=output_dir, sub_dirs=sub_dirs, output_times=output_times )

#====================================================================

if __name__ == "__main__":

    # Import utils- and plot-specific modules
    from modules_plot import *

    main()
