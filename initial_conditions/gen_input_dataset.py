#! /usr/bin/env python

# --- load modules ---
import numpy as np
import xarray as xr
import pandas as pd
import os

# --- set names ---
sp_names = ['pp', 'mp1', 'mp2', 'mp3', 'mp4']
diat_names = ['diat1', 'diat2', 'diat3']
diaz_names = ['diaz']
#zoo_names = ['zoo1', 'zoo2', 'zoo3', 'zoo4']
zoo_names = ['zoo1', 'zoo2', 'zoo3', 'zoo4', 'zoo5', 'zoo6']

calcifiers = ['mp1','mp2','mp3']
silicifiers = diat_names

autotrophs = diaz_names + sp_names + diat_names
zooplankton = zoo_names

auto_vars = ['Chl', 'C', 'P', 'Fe']
si_vars = ['Si']
ca_vars = ['CaCO3']
zoo_vars = ['C']

# --- extract sizes ---
auto_sizes = {'ugC': {'diaz': 5.5e-5, 'pp': 1e-7, 'mp1': 2.511886e-06, 'mp2': 2.417315e-04, 'mp3': 2.326305e-02, 
                      'mp4': 2.238721, 'diat1': 2e-3, 'diat2': 4.742747e-2, 'diat3': 1.12468265}}

auto_sizes['mmolC'] = {}
for x in auto_sizes['ugC']:
    auto_sizes['mmolC'][x] = auto_sizes['ugC'][x] / 12.011 / 1000

zoo_sizes = {'mgC':{'zoo1': 1.123214e-08, 'zoo2': 1.05284e-06, 'zoo3': 1.253224e-04, 'zoo4': 0.01578,
                    'zoo5': 1.98623, 'zoo6': 250.0511}}

zoo_sizes['mmolC'] = {}
for x in zoo_sizes['mgC']:
    zoo_sizes['mmolC'][x] = zoo_sizes['mgC'][x] / 12.011

# --- set scalings ---
sp_scale_factor = 1.75
diat_scale_factor = 0.4
diaz_scale_factor = 1

n_sp = len(sp_names)
n_diat = len(diat_names)
n_zoo = len(zoo_names)
n_impcalc = len(calcifiers)

sp_scalings = np.logspace(n_sp, 1, n_sp, base=2) / np.sum(np.logspace(n_sp, 1, n_sp, base=2)) * sp_scale_factor
diat_scalings = np.logspace(n_diat, 1, n_diat, base=2) / np.sum(np.logspace(n_diat, 1, n_diat, base=2)) * diat_scale_factor
zoo_scalings = np.logspace(n_zoo, 1, n_zoo, base=2) / np.sum(np.logspace(n_zoo, 1, n_zoo, base=2)) * diaz_scale_factor
impcalc_scalings = np.logspace(n_impcalc, 1, n_impcalc, base=2) / np.sum(np.logspace(n_impcalc, 1, n_impcalc, base=2))

# --- read in dataset ---
marbl_ic = '/glade/p/cesmdata/cseg/inputdata/ocn/pop/gx1v6/ic/ecosys_jan_IC_gx1v6_20161123.nc'
if os.path.exists(marbl_ic):
    d = xr.open_dataset(marbl_ic)
else:
    subprocess.run(['wget','--no-check-certificate','--user','anonymous','--password','user@example.edu',\
    'ftp://ftp.cgd.ucar.edu/cesm/inputdata/ocn/pop/gx1v6/ic/ecosys_jan_IC_gx1v6_20161123.nc'])
    d = xr.open_dataset('ecosys_jan_IC_gx1v6_20161123.nc')

# 'small phytoplankton'
for i,names in enumerate(sp_names):
    for vars in auto_vars:
        tracer = names + vars
        print(tracer)
        d[tracer] = d['sp' + vars]
        d[tracer].values = d['sp' + vars] * sp_scalings[i]

for i,calcifier in enumerate(calcifiers):
    tracer = calcifier + ca_vars[0]
    print(tracer)
    d[tracer] = d['sp' + ca_vars[0]]
    if len(calcifiers) > 1:
        d[tracer].values = d['sp' + ca_vars[0]] * impcalc_scalings[i]


# diatoms
for i,names in enumerate(diat_names):
    for vars in auto_vars + si_vars:
        tracer = names + vars
        print(tracer)
        d[tracer] = d['diat' + vars]
        d[tracer].values = d['diat' + vars] * diat_scalings[i]

# check diazotrophs
if len(diaz_names) != 1:
    print("ERROR: Need additional initial values for diazotrophs")

# zooplankton
for i,names in enumerate(zoo_names):
    tracer = names + zoo_vars[0]
    print(tracer)
    d[tracer] = d['zoo' + zoo_vars[0]]
    d[tracer].values = d['zoo' + zoo_vars[0]] * zoo_scalings[i]


# check total biomass
ds = d # make a copy
varnames = [v + 'C' for v in autotrophs+zooplankton]
varnames_zint = [v + '_zint' for v in varnames]

for i,v in enumerate(varnames):
    vname_zint = varnames_zint[i]
    ds[vname_zint] = 10 * ds[v].sum(dim='z_t')
    ds[vname_zint].attrs['units'] = 'mmol/m^2'


pft_C = {}

for i,v in enumerate(varnames_zint):
    vals = ds[varnames_zint[i]] * (ds.TAREA / 10000) #TAREA in cm^2
    vals = vals.sum(dim=['nlat', 'nlon'])
    name = varnames[i][:-1]
    pft_C[name] = vals.values

# make dataframes, and merge together
pft_C = pd.DataFrame.from_dict(pft_C, orient='index')
pft_C.columns = ['tot_biomass_mmolC']

pft_size = pd.DataFrame.from_dict(auto_sizes['mmolC'], orient='index')
pft_size = pft_size.append(pd.DataFrame.from_dict(zoo_sizes['mmolC'], orient='index'))
pft_size.columns = ['size_mmolC']

size_df = pd.merge(pft_size.reset_index(), pft_C.reset_index())

size_df['log_mmolC'] = np.log(size_df.size_mmolC)
size_df['log_biomass_mmolC'] = np.log(size_df.tot_biomass_mmolC)


# write to netcdf
d.to_netcdf('ecosys_jan_IC_gx1v6_9p6z_20180124_3impCalc.nc', mode='w')

# note: remember to convert file to ncdf3
# ncks -3 foo4c.nc foo3.nc

