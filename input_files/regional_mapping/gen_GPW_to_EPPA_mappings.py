#!/usr/bin/env python3

# 3/14/21: this code is adapted from Seb's to reflect EPPA7 regions (which separate Indonesia-IDZ and Korea-KOR as separate regions)

import pickle
import numpy as np
import gcgridobj
import xarray
import netCDF4
import os
import datetime

import pandas

# Generate spatial mappings
region_names = ['USA','CAN','MEX','JPN','ANZ','EUR','ROE','RUS','ASI','CHN','IND','BRA','AFR','MES','LAM','REA','IDZ','KOR']
# added IDZ and KOR for EPPA7

# Load in the EPPA and GPWv4 mappings
# Need to drop the NaN column resulting from Excel adding a trailing comma
eppa_table = pandas.read_csv('/home/watkin/sebeppa/Python/EPPA7_region_map.csv').dropna(how='all',axis=1)
# Now get the GPWv4 spatial mappings
cid_table = pandas.read_csv('/home/watkin/sebeppa/Python/gpw_v4_national_identifier_grid_rev11_lookup.txt',
                            delimiter='\t')

mapping_file = 'EPPA_to_GPW.pkl'
if not os.path.isfile(mapping_file):                                 # note 3/17/21: had to change from "is_file" to "isfile"
    # We need every GPWv4 territory name to map to ONE EPPA region
    eppa_countries = list(eppa_table['Countries'].values)
    missing_countries = {}
    EPPA_to_GPW = {}
    for reg in region_names:
        # List of the country identifier which is gridded by GPW
        EPPA_to_GPW[reg] = []
    GPW_to_EPPA = {}
    for idx, row in cid_table.iterrows():
        c_name_gpw = row['NAME0']
        GPW_code = row['CIESINCODE']
        # Assume failure
        match_found = False
        if c_name_gpw in eppa_countries:
            c_name_eppa = c_name_gpw
            match_found = True
        # Try replacing 'and' with ampersand
        if (not match_found) and 'and' in c_name_gpw:
            c_name_test = c_name_gpw.replace('and','&')
            if c_name_test in eppa_countries:
                c_name_eppa = c_name_test
                match_found = True
        if (not match_found) and '&' in c_name_gpw:
            c_name_test = c_name_gpw.replace('&','and')
            if c_name_test in eppa_countries:
                c_name_eppa = c_name_test
                match_found = True
        if match_found:
            reg = eppa_table['EPPA7 regions'][eppa_countries.index(c_name_eppa)].upper() # adjusted to say "EPPA7 regions"
            EPPA_to_GPW[reg].append(GPW_code)
            GPW_to_EPPA[GPW_code] = reg
        else:
            missing_countries[c_name_gpw] = GPW_code
    n_countries = len(cid_table)
    print('Auto-matched {:d} of {:d}'.format(n_countries - len(missing_countries),n_countries))
    # 3/17/21: there were 63 missing countries

    # this line does a "second" round of auto-matching
    with open('/home/watkin/sebeppa/Python/EPPA_to_GPW_mapping.dat') as f: 
        for line in f:
            ls = line.replace('\n','').split('#')
            c_name_gpw = ls[0]
            reg = ls[1].upper()
            GPW_idx = list(cid_table['NAME0'].values).index(c_name_gpw)
            GPW_code = cid_table['CIESINCODE'].values[GPW_idx]
            EPPA_to_GPW[reg].append(GPW_code)
            GPW_to_EPPA[GPW_code] = reg
            del missing_countries[c_name_gpw]

    # 3/17/21: there were no more missing countries after that automatch
    if len(missing_countries) > 0: 
        raise ValueError('Missing countries must be filled in by hand')
        # Code to get it from the user, if necessary
        for c_name_gpw, GPW_code in missing_countries.items():
            print('Enter region code for ', c_name_gpw)
            is_ok = False
            while not is_ok:
                print('Enter a region code')
                reg = input().upper()
                is_ok = reg in region_names
            GPW_to_EPPA[GPW_code] = reg
            EPPA_to_GPW[reg].append(GPW_code)
    
    pickle.dump(GPW_to_EPPA, open( "GPW_to_EPPA.pkl", "wb" ))
    pickle.dump(EPPA_to_GPW, open( mapping_file, "wb" ))
