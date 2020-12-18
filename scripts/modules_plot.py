#!/usr/bin/env python3

import argparse
import logging
import matplotlib
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import numpy as np
import os, sys, glob, shutil, re
import matplotlib.ticker as ticker
import pandas as pd
import matplotlib.transforms as transforms
import mmap
import pathlib
import errno
import json
import subprocess
import fileinput # https://kaijento.github.io/2017/05/28/python-replacing-lines-in-file/
import math
import importlib.util
import pickle as pkl
from datetime import datetime
from scipy import interpolate
from decimal import Decimal
from scipy.integrate import solve_ivp
from scipy import stats
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import numpy as np
import seaborn as sns
from matplotlib import cm

# Color definitions, https://chrisalbon.com/python/seaborn_color_palettes.html
qgray       = "#768E95"
qblue       = "#4283A9" # http://www.color-hex.com/color/4283a9
qgreen      = "#62B4A9" # http://www.color-hex.com/color/62b4a9
qred        = "#E6767A"
qturq       = "#2EC0D1"
qmagenta    = "#9A607F"
qyellow     = "#EBB434"
qgray_dark  = "#465559"
qblue_dark  = "#274e65"
qgreen_dark = "#3a6c65"
qred_dark   = "#b85e61"
qturq_dark  = "#2499a7"
qmagenta_dark = "#4d303f"
qyellow_dark  = "#a47d24"
qgray_light  = "#acbbbf"
qblue_light  = "#8db4cb"
qgreen_light = "#a0d2cb"
qred_light   = "#eb9194"
qturq_light  = "#57ccda"
qmagenta_light = "#c29fb2"
qyellow_light = "#f1ca70"

from matplotlib import cm

try:
    import seaborn as sns
except:
    print("Seaborn package not installed.")

# https://matplotlib.org/tutorials/colors/colormaps.html

no_colors   = 9

vol_zorder  = {
    "H2O"            : 11,
    "CO2"            : 10,
    "H2"             : 9,
    "CH4"            : 8,
    "N2"             : 7,
    "N2_reduced"     : 7,
    "O2"             : 5,
    "CO"             : 4,
    "S"              : 3,
    "He"             : 2,
    "NH3"            : 1,
}

vol_colors  = {
    "H2O"            : cm.get_cmap('PuBu', no_colors)(range(no_colors)),
    "CO2"            : cm.get_cmap("Reds", no_colors)(range(no_colors)),
    "H2"             : cm.get_cmap("Greens", no_colors)(range(no_colors)),
    "N2"             : cm.get_cmap("Purples", no_colors)(range(no_colors)),
    "N2_reduced"     : cm.get_cmap("Purples", no_colors)(range(no_colors)),
    "O2"             : cm.get_cmap("Wistia", no_colors+2)(range(no_colors+2)),
    "CH4"            : cm.get_cmap("RdPu", no_colors)(range(no_colors)),
    "CO"             : cm.get_cmap("pink_r", no_colors)(range(no_colors)),
    "S"              : cm.get_cmap("YlOrBr", no_colors)(range(no_colors)),
    "He"             : cm.get_cmap("Greys", no_colors)(range(no_colors)),
    "NH3"            : cm.get_cmap("cool", no_colors)(range(no_colors)),
    "greys"          : cm.get_cmap("Greys", no_colors)(range(no_colors)),
    "mixtures"       : cm.get_cmap("Set3", 9)(range(no_colors)),
    "H2O-CO2"        : cm.get_cmap("Set3", 9)(range(no_colors))[1],
    "CO2-H2O"        : cm.get_cmap("Set3", 9)(range(no_colors))[1],
    "H2O-H2"         : cm.get_cmap("Set3", 9)(range(no_colors))[2],
    "H2-H2O"         : cm.get_cmap("Set3", 9)(range(no_colors))[2],
    "H2-CO"          : cm.get_cmap("Set3", 9)(range(no_colors))[3],
    "CO-H2"          : cm.get_cmap("Set3", 9)(range(no_colors))[3],
    "H2-CO2"         : cm.get_cmap("Set3", 9)(range(no_colors))[4],
    "CO2-H2"         : cm.get_cmap("Set3", 9)(range(no_colors))[4],
    "H2-CH4"         : cm.get_cmap("Set3", 9)(range(no_colors))[5],
    "CH4-H2"         : cm.get_cmap("Set3", 9)(range(no_colors))[5],
    "H2-N2"          : cm.get_cmap("Set2", 9)(range(no_colors))[0],
    "N2-H2"          : cm.get_cmap("Set2", 9)(range(no_colors))[0],
    "CO2-N2"         : cm.get_cmap("Set2", 9)(range(no_colors))[1],
    "N2-CO2"         : cm.get_cmap("Set2", 9)(range(no_colors))[1],
    "black_1"        : "#000000",
    "black_2"        : "#323232",
    "black_3"        : "#7f7f7f",
    "H2O_1"          : "#8db4cb",
    "H2O_2"          : "#4283A9",
    "H2O_3"          : "#274e65",
    "CO2_1"          : "#811111",
    "CO2_2"          : "#B91919",
    "CO2_3"          : "#ce5e5e",
    "H2_1"           : "#a0d2cb",
    "H2_2"           : "#62B4A9",
    "H2_3"           : "#3a6c65",
    "CH4_1"          : "#eb9194",
    "CH4_2"          : "#E6767A",
    "CH4_3"          : "#b85e61",
    "CO_1"           : "#eab597",
    "CO_2"           : "#DD8452",
    "CO_3"           : "#844f31",
    "N2_1"           : "#c29fb2",
    "N2_2"           : "#9A607F",
    "N2_3"           : "#4d303f",  
    "S_1"            : "#f1ca70",
    "S_2"            : "#EBB434",
    "S_3"            : "#a47d24",    
    "O2_1"           : "#57ccda",
    "O2_2"           : "#2EC0D1",
    "O2_3"           : "#2499a7",
    "He_1"           : "#acbbbf",
    "He_2"           : "#768E95",
    "He_3"           : "#465559",
    "qgray"          : "#768E95",
    "qgray2"         : "#888888",
    "qblue"          : "#4283A9", # http://www.color-hex.com/color/4283a9
    "qgreen"         : "#62B4A9", # http://www.color-hex.com/color/62b4a9
    "qred"           : "#E6767A",
    "qturq"          : "#2EC0D1",
    "qorange"        : "#ff7f0e",
    "qmagenta"       : "#9A607F",
    "qyellow"        : "#EBB434",
    "qgray_dark"     : "#465559",
    "qblue_dark"     : "#274e65",
    "qgreen_dark"    : "#3a6c65",
    "qred_dark"      : "#b85e61",
    "qturq_dark"     : "#2499a7",
    "qmagenta_dark"  : "#4d303f",
    "qyellow_dark"   : "#a47d24",
    "qgray_light"    : "#acbbbf",
    "qblue_light"    : "#8db4cb",
    "qgreen_light"   : "#a0d2cb",
    "qred_light"     : "#eb9194",
    "qturq_light"    : "#57ccda",
    "qmagenta_light" : "#c29fb2",
    "qyellow_light"  : "#f1ca70",
}

# Volatile Latex names
vol_latex = {
    "H2O"     : r"H$_2$O",
    "CO2"     : r"CO$_2$",
    "H2"      : r"H$_2$" ,
    "CH4"     : r"CH$_4$",
    "CO"      : r"CO",
    "N2"      : r"N$_2$",
    "N2_reduced" : r"N$_2^{-}$",
    "S"       : r"S",
    "O2"      : r"O$_2$",
    "He"      : r"He",
    "NH3"     : r"NH$_3$",
    "H2O-CO2" : r"H$_2$O–CO$_2$",
    "H2O-H2"  : r"H$_2$O–H$_2$",
    "H2O-CO"  : r"H$_2$O–CO",
    "H2O-CH4" : r"H$_2$O–CH$_4$",
    "H2O-N2"  : r"H$_2$O–N$_2$",
    "H2O-O2"  : r"H$_2$O–O$_2$",
    "H2-H2O"  : r"H$_2$–H$_2$O",
    "H2-CO"   : r"H$_2$–CO",
    "H2-CH4"  : r"H$_2$–CH$_4$",
    "H2-CO2"  : r"H$_2$–CO$_2$",
    "H2-N2"   : r"H$_2$–N$_2$",
    "H2-O2"   : r"H$_2$-O$_2$",
    "CO2-N2"  : r"CO$_2$–N$_2$",
    "CO2-H2O" : r"CO$_2$–H$_2$O",
    "CO2-CO"  : r"CO$_2$–CO",
    "CO2-CH4"  : r"CO$_2$–CH$_4$",
    "CO2-O2"  : r"CO$_2$–O$_2$",
    "CO2-H2"  : r"CO$_2$–H$_2$",
    "CO-H2O" : r"CO–H$_2$O",
    "CO-CO2" : r"CO–CO$_2$",
    "CO-H2"  : r"CO–H$_2$",
    "CO-CH4" : r"CO–CH$_4$",
    "CO-N2"  : r"CO–N$_2$",
    "CO-O2"  : r"CO–O$_2$",
    "CH4-H2O" : r"CH$_4$–H$_2$O",
    "CH4-CO2" : r"CH$_4$–CO$_2$",
    "CH4-H2"  : r"CH$_4$–H$_2$",
    "CH4-CO"  : r"CH$_4$–CO",
    "CH4-CH4" : r"CH$_4$–CH$_4$",
    "CH4-N2"  : r"CH$_4$–N$_2$",
    "CH4-O2"  : r"CH$_4$–O$_2$",
    "N2-H2O" : r"N$_2$–H$_2$O",
    "N2-CO2" : r"N$_2$–CO$_2$",
    "N2-H2"  : r"N$_2$–H$_2$",
    "N2-CO"  : r"N$_2$–CO",
    "N2-CH4" : r"N$_2$–CH$_4$",
    "N2-N2"  : r"N$_2$–N$_2$",
    "N2-O2"  : r"N$_2$–O$_2$",
    "O2-H2O" : r"O$_2$–H$_2$O",
    "O2-CO2" : r"O$_2$–CO$_2$",
    "O2-H2"  : r"O$_2$–H$_2$",
    "O2-CO"  : r"O$_2$–CO",
    "O2-CH4" : r"O$_2$–CH$_4$",
    "O2-N2"  : r"O$_2$–N$_2$",
    "O2-O2"  : r"O$_2$–O$_2$",
}

molar_mass      = {
          "H2O" : 0.01801528,           # kg mol−1
          "CO2" : 0.04401,              # kg mol−1
          "H2"  : 0.00201588,           # kg mol−1
          "CH4" : 0.01604,              # kg mol−1
          "CO"  : 0.02801,              # kg mol−1
          "N2"  : 0.028014,             # kg mol−1
          "O2"  : 0.031999,             # kg mol−1
          "SO2" : 0.064066,             # kg mol−1
          "H2S" : 0.0341,               # kg mol−1 
          "H"   : 0.001008,             # kg mol−1 
          "C"   : 0.012011,             # kg mol−1 
          "O"   : 0.015999,             # kg mol−1 
          "N"   : 0.014007,             # kg mol−1 
          "S"   : 0.03206,              # kg mol−1 
          "He"  : 0.0040026,            # kg mol−1 
          "NH3" : 0.017031,             # kg mol−1 
        }

volatile_species = [ "H2O", "CO2", "H2", "CH4", "CO", "N2", "O2", "S", "He" ]


# https://stackoverflow.com/questions/13490292/format-number-using-latex-notation-in-python
def latex_float(f):
    float_str = "{0:.2g}".format(f)
    if "e" in float_str:
        base, exponent = float_str.split("e")
        return r"${0} \times 10^{{{1}}}$".format(base, int(exponent))
    else:
        return float_str

def gravity( m, r ):

    g = 6.67408E-11*m/r**2
    return g

def AtmosphericHeight(atm, m_planet, r_planet):

    z_profile       = np.zeros(len(atm.p))
    P_s             = np.max(atm.p)
    grav_s          = gravity( m_planet, r_planet )

    # Reverse arrays to go from high to low pressure
    atm.p   = atm.p[::-1]
    atm.tmp = atm.tmp[::-1]
    for vol in atm.vol_list.keys():
        atm.mr_gas[vol] = atm.mr_gas[vol][::-1]

    # print(atm.p)

    for n in range(0, len(z_profile)-1):

        # Gravity with height
        grav_z = grav_s * ((r_planet)**2) / ((r_planet + z_profile[n])**2)

        # Mean molar mass depending on mixing ratio
        mean_molar_mass = 0
        for vol in atm.vol_list.keys():
            mean_molar_mass += molar_mass[vol]*atm.mr_gas[vol][n]

        # Temperature below present height
        T_mean_below    = np.mean(atm.tmp[n:])

        # Integration
        dz = - R_gas * T_mean_below * np.log(atm.p[n+1]/atm.p[n]) / (mean_molar_mass*grav_z)
        
        # Next height
        z_profile[n+1] = z_profile[n] + dz

    # Reverse arrays again back to normal
    atm.p   = atm.p[::-1]
    atm.tmp = atm.tmp[::-1]
    for vol in atm.vol_list.keys():
        atm.mr_gas[vol] = atm.mr_gas[vol][::-1]
    z_profile = z_profile[::-1]

    return z_profile


# Constants
R_gas           = 8.31446261815324      # J K−1 mol−1

def find_nearest(array, value):
    array   = np.asarray(array)
    idx     = (np.abs(array - value)).argmin()
    return array[idx], idx

### Constants ###

# Astronomical constants
L_sun           = 3.828e+26             # W, IAU definition
AU              = 1.495978707e+11       # m
R_gas           = 8.31446261815324      # J K−1 mol−1
M_earth         = 5.972E24              # kg
R_core_earth    = 3485000.0             # m
M_core_earth    = 1.94E24               # kg
mol             = 6.02214076e+23        # mol definition

# Elements
H_mol_mass      = 0.001008              # kg mol−1
C_mol_mass      = 0.012011              # kg mol−1
O_mol_mass      = 0.015999              # kg mol−1
N_mol_mass      = 0.014007              # kg mol−1
S_mol_mass      = 0.03206               # kg mol−1
He_mol_mass     = 0.0040026             # kg mol−1
Ar_mol_mass     = 0.039948              # kg mol−1
Ne_mol_mass     = 0.020180              # kg mol−1
Kr_mol_mass     = 0.083798              # kg mol−1
Xe_mol_mass     = 0.131293              # kg mol−1

# Volatile molar masses
H2O_mol_mass    = 0.01801528            # kg mol−1
CO2_mol_mass    = 0.04401               # kg mol−1
H2_mol_mass     = 0.00201588            # kg mol−1
CH4_mol_mass    = 0.01604               # kg mol−1
CO_mol_mass     = 0.02801               # kg mol−1
N2_mol_mass     = 0.028014              # kg mol−1
O2_mol_mass     = 0.031999              # kg mol−1
SO2_mol_mass    = 0.064066              # kg mol−1
H2S_mol_mass    = 0.0341                # kg mol−1

molar_mass      = {
          "H2O" : 0.01801528,           # kg mol−1
          "CO2" : 0.04401,              # kg mol−1
          "H2"  : 0.00201588,           # kg mol−1
          "CH4" : 0.01604,              # kg mol−1
          "CO"  : 0.02801,              # kg mol−1
          "N2"  : 0.028014,             # kg mol−1
          "O2"  : 0.031999,             # kg mol−1
          "SO2" : 0.064066,             # kg mol−1
          "H2S" : 0.0341,               # kg mol−1 
          "H"   : 0.001008,             # kg mol−1 
          "C"   : 0.012011,             # kg mol−1 
          "O"   : 0.015999,             # kg mol−1 
          "N"   : 0.014007,             # kg mol−1 
          "S"   : 0.03206,              # kg mol−1 
          "He"  : 0.0040026,            # kg mol−1 
          "NH3" : 0.017031,             # kg mol−1 
        }

volatile_species = [ "H2O", "CO2", "H2", "CH4", "CO", "N2", "O2", "S", "He" ]
element_list     = [ "H", "O", "C", "N", "S", "He" ] 

# Disable and enable print: https://stackoverflow.com/questions/8391411/suppress-calls-to-print-python
def blockPrint():
    sys.stdout = open(os.devnull, 'w')
def enablePrint():
    sys.stdout = sys.__stdout__

def CleanOutputDir( output_dir ):
    types = ("*.json", "*.log", "*.csv", "*.pkl", "current??.????", "profile.*") 
    files_to_delete = []
    for files in types:
        files_to_delete.extend(glob.glob(output_dir+"/"+files))
    for file in natural_sort(files_to_delete):
        os.remove(file)

dirs = {"output": "/Users/tim/bitbucket/pcd_couple-interior-atmosphere/atm_rad_conv/output", "data_dir": "/Users/tim/bitbucket/pcd_couple-interior-atmosphere/atm_rad_conv/output/radiation_limits_data", "rad_conv": "/Users/tim/bitbucket/pcd_couple-interior-atmosphere/atm_rad_conv"}

def surf_Planck_nu(atm):
    h   = 6.63e-34
    c   = 3.0e8
    kb  = 1.38e-23
    B   = np.zeros(len(atm.band_centres))
    c1  = 1.191042e-5
    c2  = 1.4387752
    for i in range(len(atm.band_centres)):
        nu      = atm.band_centres[i]
        B[i]    = (c1*nu**3 / (np.exp(c2*nu/atm.ts)-1))
    B   = (1.-atm.albedo_s) * np.pi * B * atm.band_widths/1000.0
    return B

def get_all_output_times( odir='output' ):

    '''get all times (in Myrs) from the json files located in the
       output directory'''

    # locate times to process based on files located in odir/
    file_l = [f for f in os.listdir(odir) if os.path.isfile(os.path.join(odir,f))]
    if not file_l:
        logger.critical('output directory contains no files')
        sys.exit(0)

    time_l = [fname for fname in file_l]
    time_l = list(filter(lambda a: a.endswith('json'), time_l))

    # Filter out original/non-hacked jsons
    time_l = [ file for file in time_l if not file.startswith("orig_")]

    time_l = [int(time.split('.json')[0]) for time in time_l]
    
    # ascending order
    time_l = sorted( time_l, key=int)
    time_a = np.array( time_l )

    return time_a

def get_dict_surface_values_for_times( keys_t, time_l, indir='output'):

    '''Similar to above, but only loop over all times once and get
       all requested (surface / zero index) data in one go'''

    data_l = []

    for time in time_l:
        filename = os.path.join( indir, '{}.json'.format(time) )
        myjson_o = MyJSON( filename )
        keydata_l = []
        for key in keys_t:
            values_a = myjson_o.get_dict_values( key )
            try:
                value = values_a[0]
            except TypeError:
                value = values_a
            keydata_l.append( value )
        data_l.append( keydata_l )

    data_a = np.array( data_l )

    # rows time, cols data
    data_a.reshape( (len(time_l),-1 ) )
    # rows data, cols time
    data_a = data_a.transpose()

    return data_a

class MyJSON( object ):

    '''load and access json data'''

    def __init__( self, filename ):
        self.filename = filename
        self._load()

    def _load( self ):
        '''load and store json data from file'''
        try:
            json_data  = open( self.filename )
        except FileNotFoundError:
            logger.critical('cannot find file: %s', self.filename )
            logger.critical('please specify times for which data exists')
            sys.exit(1)
        self.data_d = json.load( json_data )
        json_data.close()

    # was get_field_data
    def get_dict( self, keys ):
        '''get all data relating to a particular field'''
        try:
            dict_d = recursive_get( self.data_d, keys )
            return dict_d
        except NameError:
            logger.critical('dictionary for %s does not exist', keys )
            sys.exit(1)

    # was get_field_units
    def get_dict_units( self, keys ):
        '''get the units (SI) of a particular field'''
        dict_d = recursive_get( self.data_d, keys )
        units = dict_d['units']
        units = None if units == 'None' else units
        return units

    # was get_scaled_field_values
    def get_dict_values( self, keys, fmt_o='' ):
        '''get the scaled values for a particular quantity'''
        dict_d = recursive_get( self.data_d, keys )
        scaling = float(dict_d['scaling'])
        if len( dict_d['values'] ) == 1:
            values_a = float( dict_d['values'][0] )
        else:
            values_a = np.array( [float(value) for value in dict_d['values']] )
        scaled_values_a = scaling * values_a
        if fmt_o:
            scaled_values_a = fmt_o.ascale( scaled_values_a )
        return scaled_values_a

    # was get_scaled_field_value_internal
    def get_dict_values_internal( self, keys, fmt_o='' ):
        '''get the scaled values for the internal nodes (ignore top
           and bottom nodes)'''
        scaled_values_a = self.get_dict_values( keys, fmt_o )
        return scaled_values_a[1:-1]

    def get_mixed_phase_boolean_array( self, nodes='basic' ):
        '''this array enables us to plot different linestyles for
           mixed phase versus single phase quantities'''
        if nodes == 'basic':
            phi = self.get_dict_values( ['data','phi_b'] )
        elif nodes == 'basic_internal':
            phi = self.get_dict_values_internal( ['data','phi_b'] )
        elif nodes == 'staggered':
            phi = self.get_dict_values( ['data','phi_s'] )
        # define mixed phase by these threshold values
        MIX = (phi<0.999) & (phi>0.001)
        MIX = MIX * 1.0 # convert to float array
        # set single phase region to nan to prevent plotting
        MIX[MIX==0] = np.nan
        return MIX

    def get_rho_interp1d( self ):
        '''return interp1d object for determining density as a
           function of pressure for static structure calculations'''
        pressure_a = self.get_dict_values( ['data','pressure_s'] )
        density_a = self.get_dict_values( ['data','rho_s'] )
        rho_interp1d = interp1d( pressure_a, density_a, kind='linear',
            fill_value='extrapolate' )
        return rho_interp1d

    def get_temp_interp1d( self ):
        '''return interp1d object for determining temperature as a
           function of pressure for static structure calculations'''
        pressure_a = self.get_dict_values( ['data','pressure_b'] )
        temp_a = self.get_dict_values( ['data','temp_b'] )
        temp_interp1d = interp1d( pressure_a, temp_a, kind='linear',
            fill_value='extrapolate' )
        return temp_interp1d

    def get_atm_struct_depth_interp1d( self ):
        '''return interp1d object for determining atmospheric height
           as a function of pressure for static structure calculations'''
        apressure_a = self.get_dict_values( ['atmosphere', 'atm_struct_pressure'] )
        adepth_a = self.get_dict_values( ['atmosphere', 'atm_struct_depth'] )
        atm_interp1d = interp1d( apressure_a, adepth_a, kind='linear' )
        return atm_interp1d

    def get_atm_struct_temp_interp1d( self ):
        '''return interp1d object for determining atmospheric temperature
           as a function of pressure'''
        apressure_a = self.get_dict_values( ['atmosphere', 'atm_struct_pressure'] )
        atemp_a = self.get_dict_values( ['atmosphere', 'atm_struct_temp'] )
        atm_interp1d = interp1d( apressure_a, atemp_a, kind='linear' )
        return atm_interp1d

def recursive_get(d, keys):

    '''function to access nested dictionaries'''

    if len(keys) == 1:
        return d[keys[0]]
    return recursive_get(d[keys[0]], keys[1:])

def get_all_output_pkl_times( odir='output' ):

    '''get all times (in Myrs) from the pkl files located in the
       output directory'''

    # locate times to process based on files located in odir/
    file_l = [f for f in os.listdir(odir) if os.path.isfile(os.path.join(odir,f))]
    if not file_l:
        logger.critical('output directory contains no PKL files')
        sys.exit(0)

    time_l = [fname for fname in file_l]
    time_l = list(filter(lambda a: a.endswith('pkl'), time_l))

    # print(time_l)

    # Filter and split files
    time_l = [ file for file in time_l if not file.startswith("orig_")]
    time_l = [ time.split('.pkl')[0] for time in time_l ]
    time_l = [ int(time.split('_atm')[0]) for time in time_l ]
    
    # ascending order
    time_l = sorted( time_l, key=int)
    time_a = np.array( time_l )

    return time_a

class FigureData( object ):

    def __init__( self, nrows, ncols, width, height, outname='fig',
        times=None, units='kyr' ):
        dd = {}
        self.data_d = dd
        if times:
            dd['time_l'] = times
            self.process_time_list()
        if units:
            dd['time_units'] = units
            dd['time_decimal_places'] = 2 # hard-coded
        dd['outname'] = outname
        self.set_properties( nrows, ncols, width, height )

    def get_color( self, num ):
        dd = self.data_d
        return dd['colors_l'][num]

    def get_legend_label( self, time ):
        dd = self.data_d
        units = dd['time_units']
        dp = dd['time_decimal_places']
        age = float(time)
        if units == 'yr':
            age = round( age, 0 )
            label = '%d'
        elif units == 'kyr':
            age /= 1.0E3
            label = '%0.1f'
        elif units == 'Myr':
            age /= 1.0E6
            label = '%0.2f'
        elif units == 'Byr' or units == 'Gyr':
            age /= 1.0E9
            label = '%0.2f'
        #label = '%0.{}e'.format( dp )
        #label = '%0.{}f'.format( dp )
        label = label % age
        return label

    def process_time_list( self ):
        dd = self.data_d
        time_l = dd['time_l']
        try:
            time_l = [int(time_l)]
        except ValueError:
            time_l = [int(time) for time in time_l.split(',')]
        self.time = time_l

    def make_figure( self ):
        dd = self.data_d
        nrows = dd['nrows']
        ncols = dd['ncols']
        fig, ax = plt.subplots( nrows, ncols )
        fig.subplots_adjust(wspace=0.3,hspace=0.3)
        fig.set_size_inches( dd['width'], dd['height'] )
        self.fig = fig
        self.ax = ax

    def savefig( self, num ):
        dd = self.data_d
        if dd['outname']:
            outname = dd['outname'] + '.pdf'
        else:
            outname = 'fig{}.pdf'.format( num)
        self.fig.savefig(outname, transparent=True, bbox_inches='tight',
            pad_inches=0.05, dpi=dd['dpi'])

    def set_colors( self, num=8, cmap='bkr8' ):
        dd = self.data_d
        # color scheme from Tim.  Nice reds and blues
        colors_l = ['#2364A4',
                   '#1695F9',
                   '#95D5FD',
                   '#8B0000',
                   '#CD5C5C',
                   '#FA141B',
                   '#FFA07A']
        # color scheme 'bkr8' for light background from Crameri
        # see f_Colours.m at http://www.fabiocrameri.ch/visualisation.php
        # this is actually very similar (same?) as Tim's scheme above
        # used in Bower et al. (2018)
        if cmap=='bkr8' and num==3:
            colors_l = [(0.0,0.0,0.3),
                        #(0.1,0.1,0.5),
                        #(0.2,0.2,0.7),
                        (0.4,0.4,0.8),
                        #(0.8,0.4,0.4),
                        #(0.7,0.2,0.2),
                        (0.5,0.1,0.1)]#,
                        #(0.3,0.0,0.0)]
            colors_l.reverse()
        elif cmap=='bkr8' and num==5:
            colors_l = [(0.0,0.0,0.3),
                        #(0.1,0.1,0.5),
                        (0.2,0.2,0.7),
                        #(0.4,0.4,0.8),
                        (0.8,0.4,0.4),
                        #(0.7,0.2,0.2),
                        (0.5,0.1,0.1),
                        (0.3,0.0,0.0)]
            colors_l.reverse()
        elif cmap=='bkr8' and num==6:
            colors_l = [(0.0,0.0,0.3),
                        (0.1,0.1,0.5),
                        (0.2,0.2,0.7),
                        #(0.4,0.4,0.8),
                        #(0.8,0.4,0.4),
                        (0.7,0.2,0.2),
                        (0.5,0.1,0.1),
                        (0.3,0.0,0.0)]
            colors_l.reverse()
        elif cmap=='bkr8' and num==8:
            colors_l = [(0.0,0.0,0.3),
                        (0.1,0.1,0.5),
                        (0.2,0.2,0.7),
                        (0.4,0.4,0.8),
                        (0.8,0.4,0.4),
                        (0.7,0.2,0.2),
                        (0.5,0.1,0.1),
                        (0.3,0.0,0.0)]
            colors_l.reverse()
        else:
            try:
                cmap = plt.get_cmap( cmap )
            except ValueError:
                cmap = plt.get_cmap('viridis_r')
            colors_l = [cmap(i) for i in np.linspace(0, 1, num)]
        dd['colors_l'] = colors_l

    def set_properties( self, nrows, ncols, width, height ):
        dd = self.data_d
        dd['nrows'] = nrows
        dd['ncols'] = ncols
        dd['width'] = width # inches
        dd['height'] = height # inches
        # TODO: breaks for MacOSX, since I don't think Mac comes
        # with serif font.  But whatever it decides to switch to
        # also looks OK and LaTeX-like.
        font_d = {'family' : 'sans-serif',
                  #'style': 'normal',
                  #'weight' : 'bold'
                  'serif': ['Arial'],
                  'sans-serif': ['Arial'],
                  'size'   : '10'}
        mpl.rc('font', **font_d)
        # Do NOT use TeX font for labels etc.
        plt.rc('text', usetex=False)
        dd['dpi'] = 300
        dd['extension'] = 'png'
        dd['fontsize_legend'] = 8
        dd['fontsize_title'] = 10
        dd['fontsize_xlabel'] = 10
        dd['fontsize_ylabel'] = 10
        try:
            self.set_colors( len(self.time) )
        except AttributeError:
            self.set_colors( num=8 )
        self.make_figure()

    def set_myaxes( self, ax, title='', xlabel='', xticks='',
        ylabel='', yticks='', yrotation='', fmt='', xfmt='', xmin='', xmax='', ymin='', ymax='' ):
        if title:
            self.set_mytitle( ax, title )
        if xlabel:
            self.set_myxlabel( ax, xlabel )
        if xticks:
            self.set_myxticks( ax, xticks, xmin, xmax, xfmt )
        if ylabel:
            self.set_myylabel( ax, ylabel, yrotation )
        if yticks:
            self.set_myyticks( ax, yticks, ymin, ymax, fmt )

    def set_mylegend( self, ax, handles, loc=4, ncol=1, TITLE=None, **kwargs ):
        dd = self.data_d
        fontsize = self.data_d['fontsize_legend']
        # FIXME
        if not TITLE:
            legend = ax.legend(handles=handles, loc=loc, ncol=ncol, fontsize=fontsize, **kwargs )
            #units = dd['time_units']
            #title = r'Time ({0})'.format( units )
        else:
            title = TITLE
            legend = ax.legend(title=title, handles=handles, loc=loc,
                ncol=ncol, fontsize=fontsize, **kwargs)
        plt.setp(legend.get_title(),fontsize=fontsize)

    def set_mytitle( self, ax, title ):
        dd = self.data_d
        fontsize = dd['fontsize_title']
        title = r'{}'.format( title )
        ax.set_title( title, fontsize=fontsize )

    def set_myxlabel( self, ax, label ):
        dd = self.data_d
        fontsize = dd['fontsize_xlabel']
        label = r'{}'.format( label )
        ax.set_xlabel( label, fontsize=fontsize )

    def set_myylabel( self, ax, label, yrotation ):
        dd = self.data_d
        fontsize = dd['fontsize_ylabel']
        if not yrotation:
            yrotation = 'horizontal'
        label = r'{}'.format( label )
        ax.set_ylabel( label, fontsize=fontsize, rotation=yrotation )

    def set_myxticks( self, ax, xticks, xmin, xmax, fmt ):
        dd = self.data_d
        if fmt:
            xticks = fmt.ascale( np.array(xticks) )
            ax.xaxis.set_major_formatter(
                mpl.ticker.FuncFormatter(fmt))
        ax.set_xticks( xticks)
        # set x limits to match extent of ticks
        if not xmax: xmax=xticks[-1]
        if not xmin: xmin=xticks[0]
        ax.set_xlim( xmin, xmax )

    def set_myyticks( self, ax, yticks, ymin, ymax, fmt ):
        dd = self.data_d
        if fmt:
            yticks = fmt.ascale( np.array(yticks) )
            ax.yaxis.set_major_formatter(
                mpl.ticker.FuncFormatter(fmt))
        ax.set_yticks( yticks)
        # set y limits to match extent of ticks
        if not ymax: ymax=yticks[-1]
        if not ymin: ymin=yticks[0]
        ax.set_ylim( ymin, ymax )

import phys
## Saturation vapor pressure [Pa] for given temperature T [K]. 
def p_sat(switch,T): 
    
    # Define volatile
    if switch == 'H2O':
        e = phys.satvps_function(phys.water)
    if switch == 'CH4':
        e = phys.satvps_function(phys.methane)
    if switch == 'CO2':
        e = phys.satvps_function(phys.co2)
    if switch == 'CO':
        e = phys.satvps_function(phys.co)
    if switch == 'N2':
        e = phys.satvps_function(phys.n2)
    if switch == 'O2':
        e = phys.satvps_function(phys.o2)
    if switch == 'H2':
        e = phys.satvps_function(phys.h2)
    if switch == 'He':
        e = phys.satvps_function(phys.he)
    if switch == 'NH3':
        e = phys.satvps_function(phys.nh3)
    
    # Return saturation vapor pressure
    return e(T)