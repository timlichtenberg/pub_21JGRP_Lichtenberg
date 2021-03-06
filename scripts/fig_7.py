#!/usr/bin/env python3

# Import utils- and plot-specific modules
from modules_plot import *

#====================================================================
def plot_atmosphere( output_dir, sub_dirs ):


    width   = 12.00 #* 3.0/2.0
    height  = 6.0

    fig = plt.figure(tight_layout=True, constrained_layout=False, figsize=[width, height])

    gs = fig.add_gridspec(nrows=2, ncols=4, wspace=0.2, hspace=0.2, left=0.055, right=0.98, top=0.98, bottom=0.08)

    sns.set_style("ticks")
    sns.despine()

    ax0 = fig.add_subplot(gs[0, 0])
    ax1 = fig.add_subplot(gs[0, 1])
    ax2 = fig.add_subplot(gs[0, 2])
    ax3 = fig.add_subplot(gs[0, 3])
    ax4 = fig.add_subplot(gs[1, 0])
    ax5 = fig.add_subplot(gs[1, 1])
    ax6 = fig.add_subplot(gs[1, 2])

    sns.set_style("ticks", { 'xtick.bottom': False, 'xtick.top': False, 'ytick.left': False, 'ytick.right': False,})
    ax7 = fig.add_subplot(gs[1, 3])


    fs_legend   = 11
    fs_label    = 11

    # Initiate legends
    legend_ax7_handles = []

    output_times = [ 1e+3, 1e+4, 1e+6 ]

    for subdir in sub_dirs:

        color_idx   = 2

        if subdir == "H2":
            ax = ax0
        if subdir == "H2O":
            ax = ax1
        if subdir == "CO2":
            ax = ax2
        if subdir == "CH4":
            ax = ax3
        if subdir == "O2":
            ax = ax4
        if subdir == "CO":
            ax = ax5
        if subdir == "N2":
            ax = ax6

        ymax_sp_flux = 0
        xmax_wavenumber = 0

        legend_ax_handles = []

        for nn, output_time in enumerate( output_times ):

            if nn == 0:
                lw = 1.5
                ls = "-"
            if nn == 1:
                lw = 1.5
                ls = "--"
            if nn == 2:
                lw = 1.5
                ls = "-."
            if nn == 3:
                lw = 1.5
                ls = ":"
            if nn == 4:
                lw = 1.5
                ls = (0, (2, 1))
            ls = "-"

            data_dir = output_dir+"/"+subdir

            # Find data_time closest to output_time wanted
            atm_data_times = get_all_output_pkl_times(data_dir)
            time, time_idx = find_nearest(atm_data_times, output_time)

            print(subdir, output_time)

            time_label = latex_float(output_time)


            atm_file = data_dir+"/"+str(int(time))+"_atm.pkl"

            # Read pickle file
            atm_file_stream = open(atm_file,'rb')
            atm = pkl.load(atm_file_stream)
            atm_file_stream.close()

            Ts_label = str(int(round(atm.ts)))

            color   = vol_colors[subdir][color_idx]

            sp_flux = atm.net_spectral_flux[:,0]/atm.band_widths
          
            l1, = ax.plot( atm.band_centres, sp_flux, ls=ls, color=color, label=Ts_label, lw=lw)

            legend_ax_handles.append(l1)
            
            # Blackbody curves
            ax.plot(atm.band_centres,surf_Planck_nu(atm)/atm.band_widths, color=color, ls="--", lw=0.5, alpha=0.5, label=str(round(atm.ts))+' K blackbody')

            ymax_sp_flux = np.max( [ ymax_sp_flux, np.max(sp_flux)] )
            xmax_wavenumber = np.max( [ xmax_wavenumber, np.max(np.where(sp_flux > ymax_sp_flux/100., atm.band_centres, 0.)) ] )
            
            # Time legend    
            if ax == ax0:
                l2, = ax7.plot( 0, 0, ls=ls, color=vol_colors["greys"][color_idx], label=time_label, lw=lw)
                legend_ax7_handles.append(l2)
     
            color_idx   += 3

        ax.set_xlim( left=0, right=xmax_wavenumber )
        ax.set_ylim( bottom=0, top=ymax_sp_flux*10 )
        ax.set_yscale("symlog", linthreshy=ymax_sp_flux/1e2)

        if ax == ax0:
            ax.set_ylabel( 'Spectral flux density (W m$^{-2}$ cm)', fontsize=fs_label)
            ax.yaxis.set_label_coords(-0.18, -0.12)

        if ax == ax5:
            ax.set_xlabel( 'Wavenumber (cm$^{-1}$)', fontsize=fs_label)

        # Individual legends
        ax.legend(handles=legend_ax_handles, loc=1, ncol=1, fontsize=fs_legend, framealpha=0.3, title="$T_\mathrm{surf}$ (K)")
        ax.text(0.02, 0.985, vol_latex[subdir], color=vol_colors[subdir][5], rotation=0, ha="left", va="top", fontsize=fs_label+3, transform=ax.transAxes, bbox=dict(fc='white', ec="white", alpha=0.01, pad=0.1, boxstyle='round'))

    sns.set_style("ticks")
    sns.despine()
    sns.despine(ax=ax7, top=True, right=True, left=True, bottom=True)
    ax7.set_yticklabels([])
    ax7.set_xticklabels([])

    ax0.text(0.42, 0.45, 'Blackbody', color=qgray_dark, rotation=0, ha="right", va="top", fontsize=fs_label, transform=ax0.transAxes, bbox=dict(fc='white', ec="white", alpha=0.01, pad=0.1, boxstyle='round'))
    ax0.arrow(0.35, 0.46, 0.07, 0.1, head_width=0.0, head_length=0.0, fc=qgray_dark, ec=qgray_dark, transform=ax0.transAxes, ls="--", lw=0.5)

    legend_ax7 = ax7.legend(handles=legend_ax7_handles, loc="center", ncol=1, fontsize=fs_legend, framealpha=0.3, title=r"Time (yr)")
    ax7.add_artist(legend_ax7)
    plt.savefig('../figures/fig_7.pdf', tight_layout=True)

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

    vols    = [ "H2", "H2O", "CO2", "CH4", "O2", "N2", "CO" ]

    output_dir  = "../data/int_atm/coupler_runs/"

    print("Host directory:", output_dir)

    # Plot fixed set from above
    plot_atmosphere( output_dir=output_dir, sub_dirs=vols )

#====================================================================

if __name__ == "__main__":

    # Import utils- and plot-specific modules
    from modules_plot import *

    main()
