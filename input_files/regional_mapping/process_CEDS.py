#!/usr/bin/env python3

# Note: updated 4/2/21 for EPPA7 regions 

import time
import sdeutils
import pickle
import numpy as np
import gcgridobj
import xarray
import netCDF4
import os
import datetime
import great_circle

import pandas

def gen_region_masks(hrz_grid_GPW,hrz_grid_out,EPPA_to_GPW,
                     fill_missing_native=False,fill_missing_regridded=True,
                     fill_mix=True,verbose=False):

    # Area-conserving regridder from GPW to CEDS resolution
    regrid_GPW_to_out = gcgridobj.regrid.gen_regridder(hrz_grid_GPW,hrz_grid_out)
    
    region_names = EPPA_to_GPW.keys()
    mask_data = sdeutils.blank_nc(region_names,datetime.datetime(2000,1,1,0,0,0),
                                  time_vec=[0],time_units='hours since 2000-01-01 00:00:00',
                                  dim_names=['time','lat','lon'],hrz_grid=hrz_grid_out)
    combo_mask = None

    for idx, reg in enumerate(region_names):
        if verbose:
            print('Generating native mask for ' + reg)
        mask = np.full(cid_grid.shape,False)
        # Get the list of CIESIN codes which match this region
        code_list = EPPA_to_GPW[reg]
        for code in code_list:
            mask[cid_grid == code] = True
        if combo_mask is None:
            combo_mask = np.full(mask.shape,0.0,np.uint16)
        combo_mask[mask] = idx + 1
    
    # Now need to exclusively classify every location on Earth (!)
    if fill_missing_native:
        print('Using Voronoi fill on native data')
        voronoi_fill(combo_mask,combo_mask<1,hrz_grid_GPW,d_pcg=0.1)
    
    for idx, reg in enumerate(region_names):
        # Degrade to the CEDS grid resolution. This is then the fraction of the area which is "owned" by that region
        mask = (combo_mask == (idx+1))*1.0
        mask_data[reg][0,...] = regrid_GPW_to_out(mask)
    
    if fill_missing_regridded:
        print('Using Voronoi fill on regridded data')
        # Make a map of the largest contributor to each cell
        n_lon = hrz_grid_out['lon'].size
        n_lat = hrz_grid_out['lat'].size
        # Build a stack of the mask data
        mask_stack = np.zeros((len(region_names),n_lat,n_lon))
        for i_reg, reg in enumerate(region_names):
            mask_stack[i_reg,...] = mask_data[reg][0,...].copy()

        if fill_mix:
            # Find cells with no contributors
            missing_mask = np.max(mask_stack,axis=0) == 0
            if verbose:
                print('{:5.1%} of cells to be filled'.format(np.sum(missing_mask)/(n_lon*n_lat)))
            voronoi_fill(mask_stack,missing_mask,hrz_grid_out,d_pcg=5)
            for i_lon in range(n_lon):
                for i_lat in range(n_lat):
                    mask_stack[:,i_lat,i_lon] /= np.sum(mask_stack[...,i_lat,i_lon])
        else:
            # Get the sum. This will be used to get the "mix" in cells
            # with at least one contributor
            mask_sum = np.sum(mask_stack,axis=0)
            # Get the max. This will be used wherever 
            mask_max = np.argmax(mask_stack,axis=0)
            missing_mask = mask_sum == 0
            voronoi_fill(mask_max,missing_mask,hrz_grid_out)
            for i_lon in range(n_lon):
                for i_lat in range(n_lat):
                    if missing_mask[i_lat,i_lon]:
                        # Assign the closest dominant region
                        mask_stack[i_reg,i_lat,i_lon] = mask_max[i_lat,i_lon]
                    else:
                        mask_stack[:,i_lat,i_lon] /= mask_sum[i_lat,i_lon]
        for i_reg, reg in enumerate(region_names):
            mask_data[reg][0,...] = mask_stack[i_reg,...]

    return mask_data

def voronoi_fill(data_grid,missing_mask,hrz_grid,d_pcg=0.1):
    t0 = time.time()
    n_lon = hrz_grid['lon'].size
    n_lat = hrz_grid['lat'].size
    last_pcg = 0
    req_lat = int(np.ceil(hrz_grid['lat'].size/2.0))
    new_find = False
    for i_lat in range(req_lat):
        # Calculate distance from this [this lat,first lon] to all other points
        dist_gc_base = great_circle.gen_gc_xmat(hrz_grid,lat_idx=i_lat).astype(np.uint16)
        for lat_flip in [False, True]:
            if lat_flip:
                dist_gc = np.flip(dist_gc_base,axis=0)
                x_lat = (n_lat - 1) - i_lat
            else:
                dist_gc = dist_gc_base
                x_lat = i_lat
            if new_find:
                # SLOOOOOOW
                # Sort the elements based on the closest
                search_order = np.unravel_index(np.argsort(dist_gc, axis=None),dist_gc.shape)
                for i_lon in range(n_lon):
                    if missing_mask[x_lat,i_lon]:
                        #t0 = time.time()
                        for test_lon,test_lat in zip(np.mod(search_order[1]+i_lon,n_lon),search_order[0]):
                            #test_lon = np.mod(flat_lon + i_lon,n_lon)
                            #test_lat = flat_lat
                            if not missing_mask[test_lat,test_lon]:
                                data_grid[...,x_lat,i_lon] = data_grid[...,test_lat,test_lon]
                                break
                        #print('found in {:8.2f} s'.format(time.time()-t0))
            else:
                for i_lon in range(hrz_grid['lon'].size):
                    if missing_mask[x_lat,i_lon]:
                        #t0 = time.time()
                        # Identify the nearest VALID entry
                        # "Rotate" the distance map to the target longitude
                        #t_vec = []
                        #t_vec.append(time.time())
                        dist_temp = np.roll(dist_gc,i_lon,axis=1)
                        #t_vec.append(time.time())
                        dist_temp[missing_mask] = 60000
                        #t_vec.append(time.time()) 
                        idx_min = np.unravel_index(np.argmin(dist_temp),dist_temp.shape)
                        #t_vec.append(time.time())
                        #data_grid[x_lat,i_lon] = data_grid[idx_min]
                        # Allow for another dimension?
                        data_grid[...,x_lat,i_lon] = data_grid[...,idx_min[0],idx_min[1]]
                        #t_vec.append(time.time())
                        #print(['{:.3f}'.format(t_vec[i+1]-t_vec[i]) for i in range(len(t_vec) - 1)])
                        #print('found in {:8.2f} s'.format(time.time()-t0))
                    
        curr_pcg = (i_lat + 1) * 100.0 / req_lat
        if curr_pcg > (last_pcg + d_pcg):
            last_pcg += d_pcg
            t_elapsed = time.time() - t0
            t_rem = (t_elapsed/curr_pcg) * (100.0 - curr_pcg)
            t_elapsed_str = sdeutils.parse_time(t_elapsed)
            t_rem_str = sdeutils.parse_time(t_rem)
            print('Completed {:5.1%} in {:s} (estimate: {:s} remaining)'.format(curr_pcg/100.0,t_elapsed_str,t_rem_str))

if __name__ == '__main__':
    ceds_data_dir = '/net/geoschem/data/gcgrid/data/ExtData/HEMCO/CEDS/v2018-08' # NOTE: changed path
    
    spc_all = ['ALD2','ALK4_butanes','ALK4_hexanes','ALK4_pentanes',
               'BC','BENZ','C2H2','C2H4','C2H6','C3H8','CH2O','CH4','CO2','CO','EOH','HCOOH',
               'MEK','NH3','NO','OC','PRPE','SO2','TOLU','XYLE']
    
    spc_map = {'NMVOC': ['ALD2','ALK4_butanes','ALK4_hexanes','ALK4_pentanes','BENZ','C2H2','C2H4','C2H6','C3H8','CH2O',
                         'CH4','EOH','HCOOH','MEK','PRPE','TOLU','XYLE'],
               'SO2':   ['SO2'],
               'CO':    ['CO'],
               'AMO':   ['NH3'],
               'BC':    ['BC'],
               'OC':    ['OC'],
               'NOx':   ['NO']}
    
    spc_missing = [spc for spc in spc_all]
    for spc_set, spc_list in spc_map.items():
        for spc in spc_list:
            if spc in spc_missing:
                spc_missing.remove(spc)
    print('The following species will not be scaled: ', spc_missing)
    
    # Read in a representative file to get the grid
    yr = 2014
    spc = 'BC'
    f_path = os.path.join(ceds_data_dir,'{:04d}'.format(yr),'{:s}-em-anthro_CMIP_CEDS_{:04d}.nc'.format(spc,yr))
    #print(os.path.isfile(f_path))
    f_nc = netCDF4.Dataset(f_path,'r')
    try:
        hrz_grid_CEDS = gcgridobj.latlontools.extract_grid(f_nc)
    finally:
        f_nc.close()
        
    d_lon = np.median(np.diff(hrz_grid_CEDS['lon']))
    d_lat = np.median(np.diff(hrz_grid_CEDS['lat']))
    print('Grid size: {:d} lat x {:d} lon, or ~{:.2f}x{:.2f}'.format(
        hrz_grid_CEDS['lat'].size,hrz_grid_CEDS['lon'].size,d_lon,d_lat))
    
    # Generate spatial mappings
    # NOTE 4/2/21: added IDZ and KOR regions for EPPA7
    region_names = ['USA','CAN','MEX','JPN','ANZ','EUR','ROE','RUS','ASI','CHN','IND','BRA','AFR','MES','LAM','REA','IDZ','KOR']
    
    # Read in the lists of which GPW IDs are contained in each region
    mapping_file = 'EPPA_to_GPW.pkl'
    if os.path.isfile(mapping_file):
        print('Loading mappings from pickled file')
        stored_EtoG = pickle.load(open(mapping_file,'rb'))
    else:
        raise ValueError('Need ' + mapping_file + '. Run gen_GPW_to_EPPA_mappings first')
    
    # Now get the GPW data
    ds_GPW = netCDF4.Dataset('/home/watkin/sebeppa/Python/gpw_v4_population_density_adjusted_rev11_2pt5_min.nc','r') # NOTE: changed path
    # GPW data is flipped in latitude, so create the grid manually
    hrz_grid_GPW = gcgridobj.latlontools.gen_grid(lon_stride=360/8640,lat_stride=180/4320,half_polar=False,center_180=False)
    var_name = 'UN WPP-Adjusted Population Density, v4.11 (2000, 2005, 2010, 2015, 2020): 2.5 arc-minutes'
    cid_grid = (np.flip(ds_GPW[var_name][10,:,:],axis=0)).astype(np.int)
    mask_data = gen_region_masks(hrz_grid_GPW,hrz_grid_CEDS,stored_EtoG,fill_missing_native=False,
                                 fill_missing_regridded=True,fill_mix=True,verbose=True)
    print('say hi again')
    mask_nc_path = '/home/watkin/sebeppa/aux_data2/CEDS_EPPA_masks_filled.nc' # NOTE: changed path
    if os.path.isfile(mask_nc_path):
        print('Clobbering ' + mask_nc_path)
        os.remove(mask_nc_path)
    mask_data.to_netcdf(mask_nc_path)
