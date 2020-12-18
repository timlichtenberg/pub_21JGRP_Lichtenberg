#!/usr/bin/env python3

# Import utils- and plot-specific modules
from modules_plot import *
from utils.modules_utils import *

import matplotlib.gridspec as gridspec
import json

#====================================================================
def plot_atmosphere( output_dir, sub_dirs ):

    # article class text width is 4.7747 inches
    # http://tex.stackexchange.com/questions/39383/determine-text-width
    # logger.info( 'building stacked interior atmosphere' )

    width   = 12.00 #* 3.0/2.0
    height  = 5.0
    
    # fig_o   = su.FigureData( 1, 3, width, height, output_dir+'/compare_sfd+melt', units='kyr' )
    # fig_o.fig.subplots_adjust(wspace=0.15, hspace=0.0)

    

    # ax0 = fig_o.ax[0]#[0]
    # ax1 = fig_o.ax[1]#[1]
    # ax2 = fig_o.ax[2]#[0]
    # # ax3 = fig_o.ax[1][1]

    # Define figure with gridspec
    # https://matplotlib.org/3.2.1/tutorials/intermediate/gridspec.html
    fig = plt.figure(tight_layout=True, constrained_layout=False, figsize=[width, height])

    gs = fig.add_gridspec(nrows=1, ncols=2, wspace=0.1, hspace=0.25, left=0.055, right=0.98, top=0.98, bottom=0.09)


    ax0 = fig.add_subplot(gs[0, 0])
    ax1 = fig.add_subplot(gs[0, 1])
    # ax2 = fig.add_subplot(gs[1, 1])

    # ax0b = fig.add_subplot(gs[0, 2:4])
    # ax1b = fig.add_subplot(gs[1, 2])
    # ax2b = fig.add_subplot(gs[1, 3])

    # sns.set_style("ticks")
    # sns.despine()

    
    handle_l = [] # handles for legend

    ymax_atm_pressure   = 0
    ymin_atm_pressure   = 1000
    ymax_atm_z          = 0
    ymin_atm_z          = 0
    ymax_sp_flux        = 0
    xmin                = 0
    xmax                = 0.1

    fs_legend   = 11
    fs_label    = 11

    # Initiate legends
    legend_ax0_handles = []
    legend_ax1_handles = []
    legend_ax1_dummy_handles = []

    # Show wavelenght or wavenumber
    print_wavelength = True

    # Define smoothing length
    nsmooth = 1

    # # Load runtime helpfile
    # runtime_helpfile = pd.read_csv(output_dir+"runtime_helpfile.csv")

    
    data_file_photosphere = '../data/fig8_cff/photosphere.json'
    with open(data_file_photosphere) as json_file:
        photosphere = json.load(json_file)


    for ni, subdir in enumerate(sub_dirs):

        settings = [ "" ]

        color_idx   = 5

        if ni == 0:
            lw = 1.5
            ls = "-"
        if ni == 1:
            lw = 1.5
            ls = "--"
        if ni == 2:
            lw = 1.5
            ls = ":"
        if ni == 3:
            lw = 1.5
            ls = "-."
        if ni == 4:
            lw = 1.5
            ls = (0, (2, 1))

        ls="-"

        for nn, setting in enumerate( settings ):

            color   = vol_colors[subdir][color_idx]

            # Define data file
            data_file_cff = '../data/fig8_cff/'+str(subdir)+"_cff.txt"

            prs         = []
            cff_tot     = []
            cff_2       = []
            cff_6       = []
            cff_10      = []
            cff_lambda  = {}

            print("Read:", data_file_cff)

            with open(data_file_cff, "r") as filestream:
                next(filestream) # skip header
                filestream = [line.rstrip('\n') for line in filestream]
                for line in filestream:
                    print(line)
                    data_values = line.split(" ")
                    cff_tot.append(float(data_values[0]))
                    prs.append(float(data_values[1]))
                    cff_2.append(float(data_values[2]))
                    cff_6.append(float(data_values[3]))
                    cff_10.append(float(data_values[4]))

            cff_lambda["2"]     = cff_2
            cff_lambda["6"]     = cff_6
            cff_lambda["10"]    = cff_10    

            print(type(photosphere))
            print(photosphere["H2_z"])

            label_a = vol_latex[subdir]+": "+str(photosphere[subdir+"_z"])+" km, "+str(photosphere[subdir+"_prs_bar"])+" bar"  

            l1, = ax0.semilogy(cff_tot, prs, ls=ls, lw=lw, color=color, label=label_a)
            legend_ax0_handles.append(l1)

            ax0.plot([0.081, 1], [photosphere[subdir+"_prs_norm"], photosphere[subdir+"_prs_norm"]], ls="--", lw=lw/2., color=color, alpha=0.8)
            
            wavelength_bands = [ 2, 6, 10 ] # microns
            for w_idx, wavelength in enumerate(wavelength_bands):

                if w_idx == 0:
                    lw = 1.5
                    ls = "-"
                if w_idx == 1:
                    lw = 1.5
                    ls = "--"
                if w_idx == 2:
                    lw = 1.5
                    ls = ":"
                if w_idx == 3:
                    lw = 1.5
                    ls = "-."
                if w_idx == 4:
                    lw = 1.5
                    ls = (0, (2, 1))

                wavenumber = (1./wavelength)*1e+4 # cm

                label_b = vol_latex[subdir]+", "+str(round(wavelength, 1))+" $\mu$m"#+band_name

                cff_band = cff_lambda[str(wavelength)]

                l2, = ax1.semilogy(cff_band, prs, ls=ls, lw=lw, color=color, label=label_b) # 
                legend_ax1_handles.append(l2)

                # Dummy legend
                if subdir == sub_dirs[0]:
                    dummy, = ax1.semilogy([0], [0], ls=ls, lw=lw, color="k", label=str(round(wavelength, 1))+" $\mu$m") 
                    legend_ax1_dummy_handles.append(dummy)

    ax0.set_xlabel( r'Normalized flux contribution, $\mathcal{CF}_\mathrm{F}$ (non-dim.)', fontsize=fs_label )
    ax0.invert_yaxis()
    ax0.set_ylabel( 'Atmospheric pressure, $P/P_{\mathrm{surf}}$ (non-dim.)', fontsize=fs_label )
    ax0.set_ylim(bottom=1, top=1e-5) # , top=1e-5
    ax0.set_xlim(left=0, right=xmax)

    ax1.set_xlabel( r'Normalized flux contribution per band, $\mathcal{CF}_\mathrm{F}^{\nu}$ (non-dim.)', fontsize=fs_label )
    ax1.invert_yaxis()
    # ax1.set_yticklabels([])
    ax1.set_ylim(bottom=1, top=1e-5) # , top=1e-5
    
    try:
        sns.set_style("ticks")
        sns.despine()
    except:
        print("No seaborn.")

    # # Legend(s)
    legend_ax0 = ax0.legend(handles=legend_ax0_handles, loc=1, ncol=1, fontsize=fs_legend, framealpha=0.99, title=r"Volatile: $z$, $p$ weighted by $\mathcal{CF}_\mathrm{F}$") # , time $t$"
    ax0.add_artist(legend_ax0)
    
    # # Detailed legend
    # legend_ax1 = ax1.legend(handles=legend_ax1_handles, loc=[0.47, 0.1], ncol=2, fontsize=fs_legend, framealpha=0.3, title=r"Volatile species, wavelength $\lambda_\mathrm{c}$" )
    # ax1.add_artist(legend_ax1)

    # Only wavelengths legend
    legend_dummy = ax1.legend(handles=legend_ax1_dummy_handles, loc=1, ncol=1, fontsize=fs_legend, framealpha=0.3, title=r"Wavelength $\lambda$" )
    ax1.add_artist(legend_dummy)

    # ax2.text(0.6, 0.28, 'Mush', color=qmagenta_light, rotation=0, ha="left", va="top", fontsize=fs_label, transform=ax2.transAxes, bbox=dict(fc='white', ec="white", alpha=0.01, pad=0.1, boxstyle='round'))
    # ax2b.text(0.6, 0.28, 'Mush', color=qmagenta_light, rotation=0, ha="left", va="top", fontsize=fs_label, transform=ax2b.transAxes, bbox=dict(fc='white', ec="white", alpha=0.01, pad=0.1, boxstyle='round'))
    
    ax0.text(0.98, 0.015, 'A', color="k", rotation=0, ha="right", va="bottom", fontsize=fs_label+4, transform=ax0.transAxes, bbox=dict(fc='white', ec="white", alpha=0.01, pad=0.1, boxstyle='round'))
    ax1.text(0.98, 0.015, 'B', color="k", rotation=0, ha="right", va="bottom", fontsize=fs_label+4, transform=ax1.transAxes, bbox=dict(fc='white', ec="white", alpha=0.01, pad=0.1, boxstyle='round'))

    ax0.text(0.995, 0.31, r'$\mathcal{CF}_\mathrm{F}$–weighted'+'\n'+'photosphere', color="k", rotation=0, ha="right", va="bottom", fontsize=fs_legend, transform=ax0.transAxes, bbox=dict(fc='white', ec="white", alpha=0.01, pad=0.1, boxstyle='round'))

    # ax1.invert_xaxis()
    # ax1.yaxis.tick_right()
    # sns.despine(left=True, right=False)

    # ax0.legend( fancybox=True, framealpha=0.5, ncol=1, fontsize=fs_legend)
    # ax2.legend( fontsize=8, fancybox=True, framealpha=0.5 )
    plt.savefig('../figures/fig_8.pdf', tight_layout=True)

    # fig_o.savefig(1)
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

    # # Define which times are plotted
    # if args.times:
    #     data_times = [ int(time) for time in args.times.split(',') ]
    #     print("Snapshots:", output_times)
    # else:
    #     data_times = su.get_all_output_times(output_dir)
    #     print("Snapshots:", output_times)

    vols    = [ "H2", "H2O", "CO2", "CH4", "O2", "N2", "CO" ]
    # vols    = [ "CH4" ]

    output_dir  = "../data/int_atm/coupler_runs/"
    
    print("Host directory:", output_dir)

    # Plot fixed set from above
    plot_atmosphere( output_dir=output_dir, sub_dirs=vols )

#====================================================================

if __name__ == "__main__":

    # Import utils- and plot-specific modules
    from modules_utils import *
    from modules_plot import *
    import utils_coupler as cu
    import utils_spider as su

    main()
