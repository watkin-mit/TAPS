#!/usr/bin/env python
# coding: utf-8
# %%

# %%


'''
This script exports gridded netCDF scaling files for global chemical transport models (CTMs). 
As scaling files, they are designed to be formatted like and multiplied by the corresponding inventory.
For more details about the inventories, see https://github.com/watkin-mit/TAPS/wiki#2-external-data-sources
'''

# import packages 
import numpy as np
import gcgridobj
import xarray
import netCDF4
import os
import io
import datetime
import time
import pandas

# choose: directories
input_dir = 'input_files'
output_dir = 'scaling_output'
netcdf_dir = 'netcdf_output'


# %%


# function to read in files
def df_to_dict(df): # num_cols = the number of category columns before the data 
    t0 = time.time()
    
    # find column index of the first data point (marked by a column title of a year)
    yrs = [int(yr_col) for yr_col in df.columns if (str.isdigit(yr_col))]
    i_yr = list(df.columns).index(str(yrs[0]))
    
    df = df.fillna(0) # some places with zeroes were exported as blanks
    nested_dict = {}
    for keys, v in df.groupby(list(df.columns[:i_yr]))[df.columns[:i_yr]]:
        d = nested_dict   # restart at root
        val = v.iloc[:,i_yr:].to_numpy().flatten() # extract the time series
        for k in keys:
            if k not in d:
                d[k] = {}   # create child if missing
            parent = d
            d = d[k]        # go down in nested level
        parent[k] = val
    
    print('seconds taken for import: ' + str(time.time() - t0))
    return nested_dict,i_yr

print('reading in the scaling files')

# read in CEDS scaling file 
CEDS_scaling_file = pandas.read_csv(os.path.join(output_dir,'CEDS_scaling.csv'))
CEDS_scaling,i_yr_CEDS_scaling = df_to_dict(CEDS_scaling_file)

# read in GFED scaling file
GFED_scaling_file = pandas.read_csv(os.path.join(output_dir,'GFED_scaling.csv'))
GFED_scaling,i_yr_GFED_scaling = df_to_dict(GFED_scaling_file)

# read in CEDS inventory dict
em_CEDS_file = pandas.read_csv(os.path.join(output_dir,'em_CEDS.csv'))
em_CEDS,i_yr_em_CEDS = df_to_dict(em_CEDS_file)

# read in GFED inventory dict. 
em_GFED_file = pandas.read_csv(os.path.join(output_dir,'em_GFED.csv'))
em_GFED,i_yr_em_GFED = df_to_dict(em_GFED_file)


# %%


# to make the figures, define the parameters

# map species 
spc_map_csv = pandas.read_csv(os.path.join(input_dir,'spc_map.csv'))
spc_map_CEDS = spc_map_csv.groupby('EPPA')['CEDS'].apply(list).to_dict() 
spc_map_CEDS['VOC'] = set(spc_map_CEDS['VOC']) # remove duplicates from VOC speciation
spc_map_CEDS_GAINS = dict(zip(spc_map_csv.dropna().EPPA,spc_map_csv.dropna().GAINS))
spc_map_GFED = spc_map_csv.dropna().groupby('EPPA')['GFED'].apply(list).to_dict() 
spc_all_CEDS = list(spc_map_csv.CEDS.unique())

# map inventory sector codes to names
inv_sec_map_file = pandas.read_csv(os.path.join(input_dir,'sectoral_mapping_EPPA7_inventories.csv'))
inv_sec_code_to_name = inv_sec_map_file.groupby('Inventory_Code')['Inventory_Name'].apply(list).to_dict()
for key, values in inv_sec_code_to_name.items():
    inv_sec_code_to_name[key] = list(set(inv_sec_code_to_name[key]))
# subset for GFED (which starts with DM), reversing the key-value order to match scale_emissions.py
GFED_sec_dict = {values[0]:key for (key,values) in inv_sec_code_to_name.items() if 'DM' in key}
    
# map inventory sectors to GAINS sectors
sec_map_CEDSGAINSEMF_file = pandas.read_csv(os.path.join(input_dir,'CEDS_GAINSEMF_sectoral_mapping.csv'))
sec_map_CEDSGAINSEMFfuels = sec_map_CEDSGAINSEMF_file.groupby(['CEDS2020','FuelCEDS'])['GAINS_EMF'].apply(list).to_dict()
sec_map_CEDSGAINSNH3_file = pandas.read_csv(os.path.join(input_dir,'CEDS_GAINSNH3_sectoral_mapping.csv'))
sec_map_CEDSGAINSNH3fuels = sec_map_CEDSGAINSNH3_file.groupby(['CEDS','CEDSFuel'])['GAINS'].apply(list).to_dict() # fix some of the names here
# point to the correct fuel mapping depending on the GAINS species (EMF or NH3)
spc_GAINSEMF = [spc for spc in spc_map_csv.dropna().GAINS.unique() if spc != 'NH3']
fuelmapdict = {spc_cat: 'sec_map_CEDSGAINSEMFfuels' for spc_cat in spc_GAINSEMF}
fuelmapdict['NH3'] = 'sec_map_CEDSGAINSNH3fuels'

# general lists from the scaling files
list_act_scen = CEDS_scaling_file.Activity_Scenario.unique()
list_reg = CEDS_scaling_file.Region.unique()
list_int_scen = CEDS_scaling_file.Intensity_Scenario.unique()
list_sec_CEDS = CEDS_scaling_file.Sector.unique()
list_fuel = CEDS_scaling_file.Fuel.unique()
list_sec_GFED = GFED_scaling_file.Sector.unique() 


# %%
# years for netCDF output
ref_yr = [int(yr_col) for yr_col in em_CEDS_file.columns if str.isdigit(yr_col)][0]

# choose: scaling years to output, and find their indices in the scaling time series 
yr_to_scale = [2050] # example
i_yr_to_scale = [list(CEDS_scaling_file.columns).index(str(yr)) 
                 - i_yr_CEDS_scaling for yr in yr_to_scale] 

yr_out = [ref_yr]+yr_to_scale


# %%


# export CEDS netCDFs

# variables that won't change during the loop
CEDS2020_data_dir = '/net/geoschem/data/gcgrid/data/ExtData/HEMCO/CEDS/v2020-08'
mask_data_CEDS = netCDF4.Dataset(os.path.join(input_dir,'CEDS_EPPA_masks_filled.nc'),'r')
n_lat = mask_data_CEDS['lat'].size
n_lon = mask_data_CEDS['lon'].size
region_names = [x for x in mask_data_CEDS.variables.keys() if x not in ['lon','lat','time']]
yr_dt = datetime.datetime(1950,1,1,0,0,0)
t_units = 'days since ' + yr_dt.strftime('%Y-%m-%d %H:%M:%S')

# for the ref_yr, multiply by 1 
reg = 'AFR' # example region for formatting
scale_grid_ref_yr = (mask_data_CEDS[reg]*np.array(0.0))+np.array(1.0)
scale_grid_ref_yr = np.tile(scale_grid_ref_yr,(12,1,1)) # convert to a monthly version

# loop for the outputs
t0 = time.time()
num_scen_yr = 1
tot_scen_yr = len(list_act_scen)*len(list_int_scen)*len(yr_out)
for act_scen in list_act_scen:
    for int_scen in list_int_scen:
        # one directory for every policy scenario (activity scenario x intensity scenario)
        scen_dir = os.path.join(netcdf_dir,'CEDS',act_scen+'_'+int_scen)
        if not os.path.isdir(scen_dir):
            os.mkdir(scen_dir)
        # Loop over each species category
        for spc_cat in spc_map_CEDS.keys():
            # loop over each CEDS species in the category
            for spc_CEDS in spc_map_CEDS[spc_cat]:
                for i_yr,yr in enumerate(yr_out): 
                    # if it's the first species, print out a progress report
                    if spc_cat == list(spc_map_CEDS.keys())[0]:
                        print('working on scen-yr',num_scen_yr,'of',tot_scen_yr)
                
                    # subdirectory for each year we scale (matching CEDS input format)
                    print('Processing {:d} {:s}'.format(yr,spc_CEDS))
                    yr_dir = os.path.join(scen_dir,'{:04d}'.format(yr))
                    if not os.path.isdir(yr_dir):
                        os.mkdir(yr_dir)
                    t_vec = [(datetime.datetime(yr,x,1,0,0,0) - yr_dt).days for x in range(1,13)] # monthly version

                    for fuel in list_fuel:
                        # read in emissions data, different for each fuel and species 
                        CEDS_path = os.path.join(CEDS2020_data_dir,'{:04d}'.format(ref_yr),'{:s}-em-{:s}_CEDS_{:04d}.nc'.format(spc_CEDS,fuel,ref_yr))
                        dsin = netCDF4.Dataset(CEDS_path,'r')

                        # prepare output file
                        dsout = netCDF4.Dataset((os.path.join(yr_dir,'{:s}-em-{:s}_CEDS_{:04d}.nc'.format(spc_CEDS,fuel,yr))), "w", format="NETCDF4_CLASSIC")
                        dsout.set_fill_off()

                        #Copy dimensions
                        for dname, the_dim in dsin.dimensions.items():
                            dsout.createDimension(dname, len(the_dim) if not the_dim.isunlimited() else None)

                        # Copy variables (sectors in this case)
                        for v_name, varin in dsin.variables.items():
                            fill_value = None
                            if hasattr(varin, "_FillValue"):
                                fill_value = varin._FillValue
                            outVar = dsout.createVariable(v_name, varin.datatype, varin.dimensions, fill_value=fill_value,
                                                          zlib=True,shuffle=True,complevel=1) # could change settings for speed-size trade-off
                            # Copy variable attributes
                            outVar.setncatts({k: varin.getncattr(k) for k in varin.ncattrs() if k not in ["_FillValue"]})

                            # Now, fill that variable with the correct value
                            if v_name == 'time':
                                outVar[:] = varin[:] - varin[:] + t_vec
                            elif v_name == 'lon' or v_name == 'lat':
                                outVar[:] = varin[:]
                            else:
                                # first, extract the sector code from the CEDS variable name format ("spc_sec")
                                sec = v_name[(len(spc_CEDS)+1):] 
                                if sec in inv_sec_code_to_name.keys():
                                    CEDS_sector = inv_sec_code_to_name[sec][0]

                                    # multiply reference year by 1, otherwise use scaling array
                                    if yr == ref_yr:
                                        outVar[:] = (varin[:]*0.0) + scale_grid_ref_yr
                                    else: 
                                        # then, get that sector's data (yr - 1 since the scaling excludes the reference year)
                                        scale_grid = 0.0
                                        if (CEDS_sector, fuel) in locals()[fuelmapdict[spc_map_CEDS_GAINS[spc_cat]]].keys():
                                            for reg in region_names:
                                                scale_grid += np.array(CEDS_scaling[act_scen][int_scen][spc_CEDS][reg][CEDS_sector][fuel][i_yr_to_scale[i_yr-1]]) * mask_data_CEDS[reg] 
                                            scale_grid = np.tile(scale_grid,(12,1,1)) # monthly version
                                        outVar[:] = (varin[:]*0.0) + scale_grid

                        dsin.close()
                        dsout.close() 
                    # estimate how much time is left: time of one scenario-year * the number of scenario-years left
                    if num_scen_yr == 1:
                        t_scen_yr = (time.time() - t0)*len(spc_all_CEDS) # time it takes to complete one scenario-year
                    if spc_cat == list(spc_map_CEDS.keys())[-1]:
                        print('estimated minutes to completion:',((t_scen_yr)*(tot_scen_yr-num_scen_yr) )/60)
                        num_scen_yr += 1
                    
mask_data_CEDS.close()


# %%


# output GFED netCDF

GFED_data_dir2014 = '/net/geoschem/data/gcgrid/data/ExtData/HEMCO/GFED4/v2015-10'

mask_data_GFED = netCDF4.Dataset(os.path.join(input_dir,'GFED_EPPA_masks_filled.nc'),'r')
yr_dt = datetime.datetime(1985,1,1,0,0,0) # note that GFED's "since year X" is different than CEDS'
days_to_hours = 24
fuel = 'process'
spc_GFED = 'BC' # GFED data is in terms of "dry matter", so all species have the same scaling

# for reference year multiplying by 1
reg = 'AFR' # example region for formatting
scale_grid_one = (mask_data_GFED[reg]*np.array(0.0))+np.array(1.0)

# loop for the outputs
t0 = time.time()
for act_scen in list_act_scen:
    for int_scen in list_int_scen:
        # one directory for every policy scenario (activity scenario x intensity scenario)
        scen_dir = os.path.join(netcdf_dir,'GFED',act_scen+'_'+int_scen)
        if not os.path.isdir(scen_dir):
            os.mkdir(scen_dir)
            
        # note that GFED trends are the same across species (as with the "dry matter" inventory)
        for i_yr,yr in enumerate(yr_out):            
            # subdirectory for each year we scale (matching GFED input format)
            print('Processing {:d}'.format(yr))
            yr_dir = os.path.join(scen_dir,'{:04d}'.format(yr))
            if not os.path.isdir(yr_dir):
                os.mkdir(yr_dir)
                
            # prepare the scaling (yr - 1 since the scaling file doesn't include the reference year)
            if yr != ref_yr:
                scale_grid = {}
                # different scaling for each sector that we scale 
                for sec_name,sec_code in GFED_sec_dict.items():
                    scale_grid[sec_code] = 0.0
                    for reg in list_reg:
                        scale_grid[sec_code] += np.array(
                            GFED_scaling[act_scen][int_scen][spc_GFED][reg][sec_name][fuel][i_yr_to_scale[i_yr-1]]) * mask_data_GFED[reg]
                            
            # file per month
            for mo in range(1,13):   
                t_vec = [(datetime.datetime(yr,mo,1,0,0,0) - yr_dt).days] # monthly version

                # generate the right file and its metadata
                GFED_path = os.path.join(GFED_data_dir2014,str(ref_yr),'GFED4_gen.025x025.{:04d}{:02d}.nc'.format(ref_yr,mo)) 
                dsin = netCDF4.Dataset(GFED_path,'r')
                dsout = netCDF4.Dataset((os.path.join(yr_dir,'GFED4_gen.025x025.{:04d}{:02d}.nc'.format(yr,mo))), "w", format="NETCDF4_CLASSIC")

                #Copy dimensions
                for dname, the_dim in dsin.dimensions.items():
                    dsout.createDimension(dname, len(the_dim) if not the_dim.isunlimited() else None)

                # Copy variables (by sector in this case)
                for v_name, varin in dsin.variables.items():
                    fill_value = None
                    if hasattr(varin, "_FillValue"):
                        fill_value = varin._FillValue
                    outVar = dsout.createVariable(v_name, varin.datatype, varin.dimensions, fill_value=fill_value,
                                                  zlib=True,shuffle=True,complevel=1)
                    # Copy variable attributes
                    outVar.setncatts({k: varin.getncattr(k) for k in varin.ncattrs() if k not in ["_FillValue"]})

                    # and fill the values! 
                    if v_name == 'lon' or v_name == 'lat':
                        outVar[:] = varin[:]
                    elif v_name == 'time':
                        outVar[:] = varin[:] - varin[:] + np.array(t_vec)*(days_to_hours)
                    elif v_name in GFED_sec_dict.values(): # if we scale this sector 
                        # multiply 1 for reference year, otherwise multiply by scaling values
                        if yr == ref_yr:
                            outVar[:] = (varin[:]*0.0) + scale_grid_one
                        else:
                            outVar[:] = (varin[:]*0.0) + scale_grid[v_name]   
                    else: 
                        # if we don't scale this sector, just multiply by it by one 
                        outVar[:] = (varin[:]*0.0) + scale_grid_one

                dsin.close()    
                dsout.close()
mask_data_GFED.close()            

print('seconds taken: ' + str(time.time() - t0))
print('export complete')


# %%




