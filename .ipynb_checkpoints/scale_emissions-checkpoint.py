#!/usr/bin/env python
# coding: utf-8
# %%

# %%


'''
Welcome! This script is organized as follows:

0.	Inputs
a.	Input parameters with activities
b.	Input parameters with intensities
c.  Input parameters with inventories

1.	Inventories
a.	Read CEDS inventory
b.	Read GFED inventory
c.	Aggregations and exports

2.	Scaling
a.	Calculate activity trends
b.	Calculate intensity trends
c.	Create intensity scenarios 

3.	Outputs
a.	Export intensity outputs
b.	Export scaling of CEDS inventory
c.	Export scaling of GFED inventory

Key user inputs are marked by the keyword "choose:"
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


# %%


#~# 0: Inputs

#~# 0(a): Input parameters with activities ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



#~# input_location ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# folder for input data and where the output is to go
input_dir = 'input_files'
output_dir = 'scaling_output'

#~# input_activity_scenarios ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Read in EPPA model energy outputs, which will be inputs for activity scaling
EPPA_energy = pandas.read_csv(os.path.join(input_dir,'EPPA7_Sectoral_energy_20210903.csv'))

# choose: scenarios considered for emitting activities (default is all in file; based on climate policy extent)
act_scen = EPPA_energy.scenario.unique()
print('climate policy scenarios considered:',act_scen)

#~# input_years ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# choose: year inputs (default is all in file; could change based on research question)
yr_list = list(np.sort(EPPA_energy.year.unique()))
ref_yr = yr_list[0]
end_yr = yr_list[-1]

yr_list = [yr for yr in yr_list if yr <= end_yr] # limit year list and EPPA datasets to end year 
yrs_from_ref = list(np.array(yr_list)-ref_yr) 

#~# format_activity_scenarios ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# EPPA regions and gridding -- produced via gen_GPW_to_EPPA_mappings.py and process_CEDS.py in input_files/regional_mapping  
region_names = [var for var in EPPA_energy.region.unique()]
xmask_data_CEDS = xarray.open_dataset(os.path.join(input_dir,'CEDS_EPPA_masks_filled.nc'))
xmask_data_GFED = xarray.open_dataset(os.path.join(input_dir,'GFED_EPPA_masks_filled.nc'))

# set the index as the year so you can fill missing values later on
EPPA_energy = EPPA_energy[EPPA_energy.year <= end_yr]
EPPA_energy = EPPA_energy.set_index('year')
# certain zero values are included as 'Eps' ... change those to zero
EPPA_energy['seuse (EJ)'] = pandas.to_numeric(EPPA_energy['seuse (EJ)'],errors = 'coerce').astype(float)
EPPA_energy['seuse (EJ)'] = EPPA_energy['seuse (EJ)'].fillna(0)

# EPPA 'energy' categories to scale CEDS 'fuels'
ene_to_fuels = {'total-coal':['COAL'],
            'solid-biofuel':['COAL','GAS','ROIL','bio','hydro','nuclear','renewables'],
            'liquid-fuel-plus-natural-gas':['GAS','ROIL'],
              'process':['COAL','GAS','ROIL','bio','hydro','nuclear','renewables']
           }

# now, read in EPPA non-energy output for population and land use scaling
EPPA_other = pandas.read_csv(os.path.join(input_dir,'EPPA7_nonsector_results.csv'))

# reformat and extract the variables you need 
EPPA_other = EPPA_other[EPPA_other.year <= end_yr]
EPPA_other.set_index('year',inplace=True)
pop = '04_Population (million people)' # for population scaling
landdict = {'CROP':'31_landuse_Cropland', # for land use scaling
           'LIVE':'33_landuse_Pasture',
           'FORS':'34_landuse_Managed forest'
           }
# certain zero values are included as 'Eps' ... change those to zero
EPPA_other.Outlook = pandas.to_numeric(EPPA_other.Outlook,errors = 'coerce').astype(float).fillna(0)
# map scenarios from EPPA_energy to EPPA_other format
scen_map_EPPA_file = pandas.read_csv(os.path.join(input_dir,'scen_map_EPPA.csv'))
scen_map_EPPA = dict(zip(scen_map_EPPA_file.EPPA_energy,scen_map_EPPA_file.EPPA_other))


# %%


#~# 0(b): Input parameters with intensities ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~# input_GAINSEMF ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# files 
GAINSEMF_em_file = 'G102021export_EMISS_LIMITS_2.csv'   # emissions
GAINSEMF_act_file = 'G102021ACTIVITY_LIMITS_2.csv'      # activities

# sectors (GAINS activities to GAINS emissions)
sec_map_GAINSGAINS_file = pandas.read_csv(os.path.join(input_dir,'GAINSEMF_sectoral_mapping.csv'))
sec_map_GAINSGAINS = dict(zip(sec_map_GAINSGAINS_file.Emissions,sec_map_GAINSGAINS_file.Activities))

# sectors (GAINS emissions to inventory)
sec_map_inv_GAINS_file = pandas.read_csv(os.path.join(input_dir,'CEDS_GAINSEMF_sectoral_mapping.csv'))
sec_map_inv_GAINS = sec_map_inv_GAINS_file.groupby('CEDS2020')['GAINS_EMF'].apply(list).to_dict()
sec_map_inv_GAINSfuels = sec_map_inv_GAINS_file.groupby(['CEDS2020','FuelCEDS'])['GAINS_EMF'].apply(list).to_dict()

# regions (from other mapping doc)
reg_map_EPPAGAINSEMF_file = pandas.read_csv(os.path.join(input_dir,'EPPA7_GAINSEMF_regional_mapping.csv'),usecols=['EPPA7','EMF30REGION'],nrows=25)
reg_GAINSEMF = reg_map_EPPAGAINSEMF_file['EMF30REGION'].unique() # list of all GAINSEMF regions
reg_map_EPPAGAINSEMF = reg_map_EPPAGAINSEMF_file.groupby('EPPA7')['EMF30REGION'].apply(list).to_dict()


#~# input_GAINSNH3 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# NH3 files (n)
fn = 'G102021export_NH3_by_source_G20.csv'
n = pandas.read_csv(os.path.join(input_dir,fn))
n = n.set_index('IDYEARS')
# correct any negative numbers. No NA fill needed (no blanks)
n.iloc[:,-3:] = n.iloc[:,-3:].abs()

# sectors (G20 to CEDS)
sec_map_CEDSGAINSNH3_file = pandas.read_csv(os.path.join(input_dir,'CEDS_GAINSNH3_sectoral_mapping.csv'))
sec_map_CEDSGAINSNH3fuels = sec_map_CEDSGAINSNH3_file.groupby(['CEDS','CEDSFuel'])['GAINS'].apply(list).to_dict()

# regions
reg_map_EPPAGAINSNH3_file = pandas.read_csv(os.path.join(input_dir,'EPPA7_GAINSG20_FAO_regional_mapping.csv'))
reg_GAINSNH3 = reg_map_EPPAGAINSNH3_file['GAINSG20'].unique() # list of all GAINSNH3 regions
reg_map_EPPAGAINSNH3 = reg_map_EPPAGAINSNH3_file.groupby('EPPA7')['GAINSG20'].apply(list).to_dict()
# subset of EPPA regions that are well-covered by G20
reg_G20 = ['CAN','USA','MEX','BRA','EUR','RUS','IND','KOR','CHN','JPN','ANZ']
# so the others would be ['LAM','AFR','ROE','MES','ASI','IDZ','REA']


#~# input_FAO ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# agriculture (FAO to G20) sectors and regions
sec_map_NH3FAO = sec_map_CEDSGAINSNH3_file.groupby(['GAINS'])[['Item','Element','Units']].apply(lambda g: g.values.tolist()).to_dict()
reg_map_NH3FAO = reg_map_EPPAGAINSNH3_file.groupby(['GAINSG20'])['FAO'].apply(list).to_dict()
reg_map_EPPAFAO = reg_map_EPPAGAINSNH3_file.groupby(['EPPA7'])['FAO'].apply(list).to_dict()
reg_FAO = (reg_map_EPPAGAINSNH3_file.FAO.unique()) # list of FAO global regions

# FAO data, scenarios and years
FAOyr = [2012,2030,2050]
FAO = pandas.read_csv(os.path.join(input_dir,'FOFA2050RegionsData_all.csv')).set_index('Year').loc[FAOyr]
FAOscen = ['Business As Usual','Toward Sustainability'] # matching GAINS CLE / MFR
scen_map_EPPAFAO = dict(zip(scen_map_EPPA_file.EPPA_energy,scen_map_EPPA_file.FAO))
# multiply normalizations by this to normalize to our base year, not the FAO base year (2012)
FAOscale = [(fyr-yr_list[0])/(fyr-FAOyr[0]) for fyr in FAOyr[1:] ]
# Separate csv for total domestic production (to avoid having to do the math in subsequent loops)
FAOtotal = pandas.read_csv(os.path.join(input_dir,'FOFA2050RegionsData_all_total.csv'))
FAOtotal = FAOtotal.set_index(FAOtotal.columns[0])[list(map(str, FAOyr))].reset_index()
# Separate csv for total yields by scenario (to avoid having to do the math in subsequent loops)
FAOyield = pandas.read_csv(os.path.join(input_dir,'FOFA2050RegionsData_all_yield.csv'))
# prepare the conversion/extension from FAO years to EPPA years
FAOyieldn = FAOyield.T.iloc[1:,:]
FAOyieldn.index = pandas.to_datetime(FAOyieldn.index,format="%Y")
FAOyieldn2 = FAOyieldn.reindex(pandas.to_datetime(list(np.sort(list(set(yr_list+FAOyr)))),format="%Y"),).astype(float)
# linearly interpolate and extrapolate (based on the final timepoints)
FAOyieldint = FAOyieldn2.interpolate(method="slinear", axis=0, fill_value="extrapolate", limit_direction="both")
# and convert back into a dataframe with a readable format of EPPA years
FAOyieldint.index = FAOyieldint.index.year
FAOEPPA = pandas.concat([FAOyield.T.iloc[0:1,:],FAOyieldint]).T.drop(columns=FAOyr[0])


# %%


#~# 0(c): Input parameters with inventories ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~# Prepare CEDS_GBD-MAPS and GFED4.1s inputs ~~~~~~~~~~~~~~~~~

# CEDS inputs (McDuffie et al., 2020)
CEDS2020_data_dir = '/net/geoschem/data/gcgrid/data/ExtData/HEMCO/CEDS/v2020-08' # choose: specify your folder
fuels = ['total-coal','solid-biofuel','liquid-fuel-plus-natural-gas','process']

# species mapping
spc_map_csv = pandas.read_csv(os.path.join(input_dir,'spc_map.csv'))
spc_map_CEDS = spc_map_csv.groupby('EPPA')['CEDS'].apply(list).to_dict() 
spc_map_CEDS['VOC'] = set(spc_map_CEDS['VOC']) # remove duplicates
spc_map_GFED = spc_map_csv.dropna().groupby('EPPA')['GFED'].apply(list).to_dict() 
spc_map_CEDSGFED = spc_map_csv.groupby('CEDS')['GFED'].apply(list).to_dict()
spc_map_GAINS = dict(zip(spc_map_csv.dropna().EPPA,spc_map_csv.dropna().GAINS))

# also keep a list of all the species in each for convenience
spc_all_CEDS = list(spc_map_csv.CEDS.unique())
spc_all_GFED = list(spc_map_csv.GFED.unique())
spc_all_GFED.remove(np.nan)
spc_GAINSNH3 = ['NH3']
spc_GAINSEMF = [spc for spc in spc_map_csv.dropna().GAINS.unique() if spc not in spc_GAINSNH3]

# GFED inputs (van Marle et al., 2017). This reflects the latest source for the 2014 base year used here
GFED_data_dir = '/net/geoschem/data/gcgrid/data/ExtData/HEMCO/GFED4/v2015-10' # choose: your download folder

# compile monthly base year data 
tmplst = []
for mth in range(1,13):
    file = f'GFED4_gen.025x025.{str(ref_yr).zfill(4)}{str(mth).zfill(2)}.nc'
    tmplst.append(xarray.open_dataset(os.path.join(GFED_data_dir,str(ref_yr),file))) 
GFED = xarray.concat(tmplst,dim='time') 

# read in GFED emission factors (EF), from https://www.geo.vu.nl/~gwerf/GFED/GFED4/ancill/
GFED4_EF = pandas.read_csv(input_dir+'/'+"GFED4_Emission_Factors.txt", sep="\s+",lineterminator='\n', header=15,skiprows=[1],
                           names=['SPECIE','SAVA','BORF','TEMF','DEFO','PEAT','AGRI'])


#~# prepare the CEDS and GFED grids ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# CEDS: extract grid cell area from CEDS 0.5x0.5 example file (with example species and fuel)
spc = 'BC'
f = 'process'
CEDS_path = os.path.join(CEDS2020_data_dir,'{:04d}'.format(ref_yr),'{:s}-em-{:s}_CEDS_{:04d}.nc'.format(spc,f,ref_yr))
f_nc = netCDF4.Dataset(CEDS_path,'r')
try:
    hrz_grid_CEDS = gcgridobj.latlontools.extract_grid(f_nc)
finally:
    f_nc.close()

# GFED: extract grid cell area from GFED 0.25x0.25 example file
GFED_path = os.path.join(GFED_data_dir,str(ref_yr),'GFED4_gen.025x025.201401.nc')
f_nc = netCDF4.Dataset(GFED_path,'r')
try:
    hrz_grid_GFED = gcgridobj.latlontools.extract_grid(f_nc)
finally:
    f_nc.close()


#~# Prepare the sectoral secaling with the activity data ~~~~~~~~~~~~~~~

inv_sec_map_file = pandas.read_csv(os.path.join(input_dir,'sectoral_mapping_EPPA7_inventories.csv'))

# map inventory sector names to codes
inv_sec_map = inv_sec_map_file.groupby('Inventory_Name')['Inventory_Code'].apply(list).to_dict()
for key, values in inv_sec_map.items():
    inv_sec_map[key] = list(set(inv_sec_map[key]))
    
# map inventory sector names to which inventory it comes from
inv_sec_map_source = inv_sec_map_file.groupby('Inventory_Name')['Inventory'].apply(list).to_dict()
for key, values in inv_sec_map_source.items():
    inv_sec_map_source[key] = list(set(inv_sec_map_source[key]))

# find the subset of sectors from GFED (labeled "DM" = Dry Matter)
GFED_sec_dict = {key:value for (key,value) in inv_sec_map.items() if 'DM' in value[0]}

# map inventory sector names to activity sectors
inv_act_sec_map = inv_sec_map_file.groupby('Inventory_Name')['EPPA7'].apply(list).to_dict()


# %%


#~# 1: Inventories

#~# 1(a): Process CEDS inventory ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Read in the CEDS base-year emissions (em)
t0 = time.time()
print('getting base-year CEDS emissions')
em_CEDS = {} 
em_CEDS_export = [] # export file for other scripts
for spc in spc_map_CEDSGFED.keys():
    print('working on ',spc)
    em_CEDS[spc] = {}
    for f in fuels: 
        em_CEDS[spc][f] = {}
        CEDS_path = os.path.join(CEDS2020_data_dir,'{:04d}'.format(ref_yr),'{:s}-em-{:s}_CEDS_{:04d}.nc'.format(spc,f,ref_yr))
        CEDS = xarray.open_dataset(CEDS_path)
        for sec,values in CEDS.items():
            em_CEDS[spc][f][sec] = {}
            # use annual mean (m) for annual scaling (to keep the base year's monthly distribution, as in Feng et al. (2020) for the SSPs))
            CEDSm = values.mean(dim='time')
            # get emissions (em) from species (kg m-2 s-1 * m2 * s yr-1 * Tg kg-1 = Tg yr-1), using the Sidereal year
            CEDSem = CEDSm * hrz_grid_CEDS.area * (365.25636*24*3600) * 1e-9
            for reg in region_names:
                # use masks (specifying proportions of each grid cell in each region) to get regional data
                CEDSemgrid = CEDSem * xmask_data_CEDS[reg]
                em_CEDS[spc][f][sec][reg] = float(CEDSemgrid.sum().values)
                
                # and export inventory values
                em_CEDS_export.append([spc,f,sec,reg,em_CEDS[spc][f][sec][reg]])
                            
em_CEDS_path = os.path.join(output_dir,'em_CEDS.csv')
em_CEDS_df = pandas.DataFrame(em_CEDS_export, 
                                   columns = ['Species','Fuel','Sector',
                                              'Region',str(ref_yr)])
em_CEDS_df.to_csv(em_CEDS_path, index = False, header=True)

print('seconds taken: ' + str(time.time() - t0))
print('estimated minutes to finish script: ',round((time.time() - t0)*(11/60)))


# %%


#~# 1(b): Process GFED inventory ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# now, do the same for GFED (defined as "process" emissions)
t0 = time.time()
print('getting base-year GFED emissions')
em_GFED = {}
em_GFED_export = [] # export file for other scripts

for i_sec,(sec_name,sec_code) in enumerate(GFED_sec_dict.items()):
    em_GFED[sec_name] = {}
    GFEDm = GFED[sec_code].mean(dim='time') # use annual mean (m) for annual scaling (same as above)
    sec_code_EF = sec_code[0][3:] # emission factor code doesn't have the "DM_"
    
    for spc in spc_all_GFED:
        print('working on',spc)
        em_GFED[sec_name][spc] = {}
        EF = GFED4_EF[GFED4_EF.SPECIE == spc][sec_code_EF].squeeze()
        # get emissions from species (kg m-2 s-1 * m2 * s yr-1 * Tg kg-1 = Tg yr-1)
        GFEDem = GFEDm * EF * hrz_grid_GFED.area * (365.25636*24*3600) * 1e-12
        
        for reg in region_names:
            # allot emissions to each region based on grid cell area and mask data
            GFEDemgrid = GFEDem * xmask_data_GFED[reg]
            em_GFED[sec_name][spc][reg] = float(GFEDemgrid.sum().to_array()) 

            # and export inventory values in format that matches the CEDS order (species, region, sector)
            em_GFED_export.append([spc,reg,sec_name,em_GFED[sec_name][spc][reg]])
                            
em_GFED_path = os.path.join(output_dir,'em_GFED.csv')
em_GFED_df = pandas.DataFrame(em_GFED_export, 
                                   columns = ['Species','Region','Sector',  
                                              str(ref_yr)])
em_GFED_df.to_csv(em_GFED_path, index = False, header=True)        
        
print('seconds taken: ' + str(time.time() - t0))
                


# %%


#~# 1(c): Aggregations and exports ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#~# inventory by fuel across pollutants, to distinguish the nonzero categories for activity scaling
emsum = {}
for fuel in fuels:
    emsum[fuel] = {}
    for inv_sector in inv_act_sec_map.keys():
        emsum[fuel][inv_sector] = {}
        for reg in region_names:
            emsum[fuel][inv_sector][reg] = 0.0
            for s, (spc_cat,spcs_CEDS) in enumerate(spc_map_CEDS.items()):
                for spc_CEDS in spcs_CEDS:
                    if inv_sector in GFED_sec_dict.keys():
                        # define GFED fuel as "process"; GFED sum is 0 for other fuels
                        if fuel == 'process':
                            for spc_GFED in spc_map_CEDSGFED[spc_CEDS]:
                                if type(spc_GFED) == str:
                                    emsum[fuel][inv_sector][reg] += em_GFED[inv_sector][spc][reg]
                        else:
                            emsum[fuel][inv_sector][reg] = 0.0
                    else: 
                        sec = spc_CEDS + '_' + inv_sec_map[inv_sector][0]
                        emsum[fuel][inv_sector][reg] += em_CEDS[spc_CEDS][fuel][sec][reg]

#~# calculate proportion of species/sector/region's emissions from each CEDS fuel (for activity scaling)
fpct = {}
fpctexport = []
for spc in spc_map_CEDSGFED.keys():
    fpct[spc] = {}
    for inv_sector,sec_code in inv_sec_map.items():
        if inv_sector not in GFED_sec_dict.keys():
            sec = spc + '_' + sec_code[0]
            fpct[spc][sec] = {}
            for reg in region_names:
                fpct[spc][sec][reg] = {}
                # first, find the sum over all fuels 
                fuelsum = 0.0
                for f in fuels:
                    fuelsum += em_CEDS[spc][f][sec][reg]    
                # then, find the proportion from each fuel
                for f in fuels:
                    if fuelsum == 0.0:
                        fpct[spc][sec][reg][f] = 0.0
                    else: 
                        fpct[spc][sec][reg][f] = em_CEDS[spc][f][sec][reg] / fuelsum

                # the rest is export code if desired:
                fpctexport.append([spc,sec,reg,f,fpct[spc][sec][reg][f]])

fpctpath = os.path.join(output_dir,'fuel_proportions.csv')
fpctdf = pandas.DataFrame(fpctexport, columns = ['Species','Sector','Region','Fuel','2014_Proportion']) 
fpctdf.to_csv(fpctpath, index = False, header=True)

#~# Table A1: CEDS global fuel proportions by sector and species~~~~~~~~~~~~~~~~~~

# Goal: row per sector-fuel, column per species
# globally: sum the emissions for each sector/species, and then find the % from each fuel
TableA1_dict = {}
spc_table = ['SO2','NO','NH3','OC','BC','CO','C2H4'] # NO for NOx, and C2H4 as example VOC species
for spc in spc_table:
    TableA1_dict[spc] = {}
    for inv_sector,sec_code in inv_sec_map.items():
        if inv_sector not in GFED_sec_dict.keys():
            sec = spc + '_' + sec_code[0]
            TableA1_dict[spc][inv_sector] = {}
            # first, find the global emissions across fuels
            fuelsum = 0.0
            for reg in region_names:
                for f in fuels:
                    fuelsum += em_CEDS[spc][f][sec][reg]  
            # then, find the proportion from each fuel, if any
            for f in fuels:
                if fuelsum == 0.0:
                    TableA1_dict[spc][inv_sector][f] = 0.0
                else:
                    regionsum = 0.0
                    for reg in region_names:
                        regionsum += em_CEDS[spc][f][sec][reg]  
                    TableA1_dict[spc][inv_sector][f] = regionsum/fuelsum
                    
# then, format the table
TableA1 = []
for inv_sector,sec_code in inv_sec_map.items():
    if inv_sector not in GFED_sec_dict.keys():    
        for nf,f in enumerate(fuels):
            # start a row with the sector and fuel (naming the sector if it's the first fuel)
            if nf == 0:
                row = [inv_sector,f]
            else:
                row = ['',f]   
            # continue the row with a column per species, expressed as a whole-number percentage
            for spc in spc_table: 
                row.append(round(TableA1_dict[spc][inv_sector][f]*100))
            # append the row to the overall table
            TableA1.append(row)

# and export the table
TableA1path = os.path.join(output_dir,'TableA1.csv')
TableA1df = pandas.DataFrame(TableA1, columns = ['Sector','Fuel']+spc_table) 
TableA1df.to_csv(TableA1path, index = False, header=True)


# %%


#~# 2: Scaling

#~# 2(a): Calculate activity trends ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

t0 = time.time()
print('working on activity trends')
act_scaling = {} # EPPA activity trends  
# 3 scenarios
for i_scen, scen in enumerate(act_scen):
    print('working on ' + scen)
    act_scaling[scen] = {}
    # for agriculture production calculation
    FAOi = (FAOEPPA[FAOEPPA['Yield (tonnes/ha)'] == scen_map_EPPAFAO[scen]].index[0])
    FAOyields = FAOEPPA.loc[(FAOi+1):(FAOi+8)]
    
    # now, work with the activity trend data
    # for each CEDS sector, sum the mapped EPPA sector trends and divide by their summed base value 
    for i,(inv_sector,EPPA_sectors) in enumerate(inv_act_sec_map.items()):
        act_scaling[scen][inv_sector] = {}
        print('working on ' + inv_sector)
        for reg in region_names:
            act_scaling[scen][inv_sector][reg] = {}
            for nf,fuel in enumerate(fuels):
                act_scaling[scen][inv_sector][reg][nf] = 0.0

            # first, scale the activities mapped to population
            if 'Population' in EPPA_sectors: 
                poptrend = EPPA_other[(EPPA_other.Scenario == scen_map_EPPA[scen]) 
                                      & (EPPA_other.Parameter == pop) & 
                                      (EPPA_other.region == reg)].Outlook.astype(float)
                poptrend.reset_index(drop=True, inplace=True) 
                for nf,fuel in enumerate(fuels):   
                    # only scale if we have emissions in the inventory (saves a bit of time)
                    if emsum[fuel][inv_sector][reg] > 0:
                        act_scaling[scen][inv_sector][reg][nf] = np.array(poptrend/poptrend[0])
                    else:                
                        act_scaling[scen][inv_sector][reg][nf] = np.ones(len(yr_list)) 
                continue
            
            # next, scale agricultural (non-energy) process activities by land use
            if 'Agricult' in inv_sector: 
                nf = fuels.index('process')  # agriculture emissions are all 'process'
                # only scale if we have emissions in the inventory (saves a bit of time)
                if (emsum[fuels[nf]][inv_sector][reg] > 0):
                    for esec in EPPA_sectors:
                        trend = EPPA_other[(EPPA_other.Scenario == scen_map_EPPA[scen])
                                          &(EPPA_other.Parameter == landdict[esec])
                                          &(EPPA_other.region == reg)].Outlook.astype(float)
                        # if there are values here: 
                        if len(trend) > 0:
                            # fill in zero values (which are not included in the spreadsheet)
                            if len(trend) < len(yr_list):
                                for yr in yr_list:
                                    if yr not in trend.index:
                                        trend[yr] = 0.0        
                                trend = trend.astype(float).sort_index()
                            # sum for the trend (aggregating by total land use). Will be normalized below
                            act_scaling[scen][inv_sector][reg][nf] += np.array(trend)       

            # finally, scale other sectors by EPPA energy trends
            for nf,(fuel,energies) in enumerate(ene_to_fuels.items()):
                # only scale if we have emissions in the inventory (saves a bit of time)
                if (emsum[fuel][inv_sector][reg] > 0) and ('Agricult' not in inv_sector):
                    # save a bit of time by 'autofilling' sectors whose scaling (TRAN) was already found
                    if inv_sector == 'Transport' or inv_sector == 'Shipping':
                        act_scaling[scen][inv_sector][reg][nf] = act_scaling[scen]['Non-road transport'][reg][nf]
                        continue
                    # for each EPPA energy associated with that CEDS fuel 
                    for ene in energies: 
                        # for each EPPA sector associated with this CEDS sector
                        for esec in EPPA_sectors:
                            trend = EPPA_energy[(EPPA_energy.scenario == act_scen[i_scen])
                                               &(EPPA_energy.sector == esec)&(EPPA_energy.energy == ene)
                                               &(EPPA_energy.region == reg)]['seuse (EJ)']
                            # if there are values here: 
                            if len(trend) > 0:
                                # fill in zero values (which are not included in the spreadsheet)
                                if len(trend) < len(yr_list):
                                    for yr in yr_list:
                                        if yr not in trend.index:
                                            trend[yr] = 0.0        
                                    trend = trend.astype(float).sort_index()
                                # aggregates based on total energy (up-weighting higher-energy sectors)
                                act_scaling[scen][inv_sector][reg][nf] += np.array(trend)
                 
                # and while we're looping across fuels, finish the processing for energy and ag. scaling
                # if no emissions in inventory or no activity data from EPPA, just put a constant trend
                if type(act_scaling[scen][inv_sector][reg][nf]) is float:
                    act_scaling[scen][inv_sector][reg][nf] = np.ones(len(yr_list)) 
                
                # normalize to base year
                act_scaling[scen][inv_sector][reg][nf] /= act_scaling[scen][inv_sector][reg][nf][0]  
                
            # for agriculture, multiply land use (ha) by yield (production/ha) to match GAINS units
            if inv_sector == 'Agriculture':
                nf = fuels.index('process')
                FAOyieldtrend = FAOyields[FAOyields['Yield (tonnes/ha)'] == 
                              reg_map_EPPAFAO[reg][0]].iloc[:,1:].squeeze().to_numpy()
                act_scaling[scen][inv_sector][reg][nf] = act_scaling[scen][inv_sector][reg][nf] * (FAOyieldtrend/FAOyieldtrend[0])

print('seconds taken: ' + str(time.time() - t0))


# %%


#~# 2(b): Calculate intensity trends ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~# calculate GAINS EMF fuel-specfic emission factors

# read in sheets
EMFact = pandas.read_csv(os.path.join(input_dir,GAINSEMF_act_file))
EMFem = pandas.read_csv(os.path.join(input_dir,GAINSEMF_em_file))
# GAINS scenarios: CLE and MFR 
int_scen = EMFact.IDSCENARIOS.unique()
print('pollution policy scenarios considered:',int_scen)
int_yrs = {}           # vector per scenario: all years covered by that scenario
i_yr_EMFact = []       # item per scenario: index of activity first year in the time series 
i_yr_EMFem = []        # same for emissions file
int_scen_baseyr = []   # item per scenario: closest year to EPPA base year (for aggregations)o
for scen in int_scen:
    scen_head = EMFact[EMFact.IDSCENARIOS==scen].head().dropna(axis=1, how='all')
    int_yrs[scen] = [yr_col for yr_col in scen_head.columns if (str.isdigit(yr_col))]
    int_yrs_num = [int(yr) for yr in int_yrs[scen]]
    i_yr_EMFact.append(list(EMFact.columns).index(int_yrs[scen][0]))
    i_yr_EMFem.append(list(EMFem.columns).index(int_yrs[scen][0]))
    int_scen_baseyr.append(str(int_yrs_num[min(range(len(int_yrs_num)), 
                                               key = lambda i: abs(int_yrs_num[i]-ref_yr))]))
    
# find the year parameters to get to the base year: ref_yr = max_below * prop_to*(min_above - max_below)
max_below_ref_yr = max(int(i) for i in EMFact.columns[min(i_yr_EMFact):] if int(i) < ref_yr)
min_above_ref_yr = min(int(i) for i in EMFact.columns[min(i_yr_EMFact):] if int(i) > ref_yr)
prop_to_ref_yr = (ref_yr - max_below_ref_yr)/(min_above_ref_yr - max_below_ref_yr)

# correct a few negative numbers from the numerical columns (assuming they meant to be positive, not zero)
EMFact.iloc[:,min(i_yr_EMFact):] = EMFact.iloc[:,min(i_yr_EMFact):].abs()
EMFem.iloc[:,min(i_yr_EMFem):] = EMFem.iloc[:,min(i_yr_EMFem):].abs()

# fill the emissions NAs with zeros
EMFem.fillna(0, inplace=True)
# and linearly interpolate activity NAs (including starting/ending values), if row has 2 or more values... 
# ...first convert year columns to datetime
EMFact.columns = EMFact.columns[:min(i_yr_EMFact)].tolist() + pandas.to_datetime(
    EMFact.columns[min(i_yr_EMFact):],format="%Y").tolist() 
# and then do the interpolation
EMFact.iloc[:,min(i_yr_EMFact):] = EMFact.iloc[:,min(i_yr_EMFact):].interpolate(
    axis=1,fill_value="extrapolate",limit_direction='both') 
# and correct the negatives from interpolation to be zero
EMFact.iloc[:,min(i_yr_EMFact):][EMFact.iloc[:,min(i_yr_EMFact):] < 0] = 0
# and reset the year columns to not be datetime
EMFact.columns = EMFact.columns[:min(i_yr_EMFact)].tolist() + list(pandas.Series(
    EMFact.columns[min(i_yr_EMFact):]).dt.strftime('%Y').astype(str))

# Loop to calculate trends for GAINSEMF emissions factors (EF)
t0 = time.time()
print('working on GAINS EMF')
EF_EMF = {}       # EF trend dictionary (main output)
em_GAINSreg = {}  # emissions by GAINS region and inventory sector (to aggregate to EPPA regions in next loop)
cle_byr = {}      # CLE base year value (to normalize MFR) 
EF_EMF_raw = []   # export: list of raw EF 

# for each emissions scenario (2: CLE, MFR)
for i_scen, scen in enumerate(int_scen):
    EF_EMF[scen] = {}
    em_GAINSreg[scen] = {}
    trend_ones = (EMFem.iloc[2,(i_yr_EMFem[i_scen]):]/
                  EMFem.iloc[2,(i_yr_EMFem[i_scen]):]) # trend of ones for filling constant EF trends
    
    # for each species in GAINSEMF
    for s, spc_cat in enumerate(spc_GAINSEMF):
        print('Working on ' + spc_cat)
        EF_EMF[scen][spc_cat] = {}
        em_GAINSreg[scen][spc_cat] = {}
        if i_scen == 0:
            cle_byr[spc_cat] = {} 
        
        # for each GAINSEMF region  
        for reg in reg_GAINSEMF:
            if i_scen == 0:
                cle_byr[spc_cat][reg] = {} 
            EF_EMF[scen][spc_cat][reg] = {}
            em_GAINSreg[scen][spc_cat][reg] = {}

            # for each inventory sector-fuel combination 
            inv_sec_last = ''
            for (inv_sector,fuel),GAINS_sectors in sec_map_inv_GAINSfuels.items():
               
                if inv_sector != inv_sec_last:
                    em_GAINSreg[scen][spc_cat][reg][inv_sector] = {}
                    EF_EMF[scen][spc_cat][reg][inv_sector] = {}
                    inv_sec_last = inv_sector
                
                # first, loop over GAINS sectors in each inventory sector to make a summed emissions trend 
                em_GAINSreg[scen][spc_cat][reg][inv_sector][fuel] = trend_ones * 0.0 
                for gsec in GAINS_sectors:    
                    # find the emissions vector (ev)
                    ev = EMFem.loc[((EMFem.EMF30REGION == reg) & (EMFem.EMF30_DET == gsec) &
                                  (EMFem.IDSCENARIOS == scen) & 
                                  (EMFem.IDPOLLUTANT_FRACTIONS == spc_cat))]
                    if len(ev) > 0:
                        ev = ev.iloc[0,(i_yr_EMFem[i_scen]):] # subset the time series for simplicity
                        if sum(ev) > 0:
                            em_GAINSreg[scen][spc_cat][reg][inv_sector][fuel] += ev

                # then, loop over GAINS sectors again to calculate EF trends
                EF_EMF[scen][spc_cat][reg][inv_sector][fuel] = 0.0
                for gsec in GAINS_sectors:
                    # initialize parameters
                    gsec_EF = []                            # EF trend for this GAINS sector
                    gsecpctem = trend_ones * 0.0            # % of inv. sector emissions from GAINS sector  
                    if i_scen == 0:
                        cle_byr[spc_cat][reg][gsec] = 1.0   # CLE EF at base year

                    # find the emissions vector and % of CEDS sector emissions from that GAINS sector
                    # or set to zero if there are no emissions
                    ev = EMFem.loc[((EMFem.EMF30REGION == reg) & (EMFem.EMF30_DET == gsec) &
                                  (EMFem.IDSCENARIOS == scen) & 
                                  (EMFem.IDPOLLUTANT_FRACTIONS == spc_cat))]
                    if len(ev) > 0:
                        ev = ev.iloc[0,(i_yr_EMFem[i_scen]):] # subset the time series for simplicity
                        gsecpctem = (ev / em_GAINSreg[scen][spc_cat][reg][inv_sector][fuel]).replace(np.inf,0).fillna(0)
                                            
                    # if GAINS sector doesn't exist in activity mapping, EF is constant due to lack of info
                    # exception for Waste, which lacks activity data but is scaled further down
                    if sec_map_GAINSGAINS[gsec] == ['Blank'] and not("Waste" in gsec):
                        EF_EMF[scen][spc_cat][reg][inv_sector][fuel] += (trend_ones) * gsecpctem
                        continue   
                        
                    # if activity mapping does exist, find the activity trend 
                    # (the trend still might not exist for that region)
                    av = EMFact.loc[((EMFact.EMF30REGION == reg) & 
                                     (EMFact.SEC_UNIT == sec_map_GAINSGAINS[gsec]) &
                                     (EMFact.IDSCENARIOS == scen))]

                    # calculate EF of GAINS sector (emissions / activity) 
                    if len(ev) > 0 and len(av) > 0:
                        av = av.iloc[0,(i_yr_EMFact[i_scen]):] # subset the time series for simplicity
                        gsec_EF = (ev / av).fillna(0) # set an EF point to zero if its activity is zero

                        # calculate the CLE base year value (to normalize MFR)
                        if i_scen == 0:
                            cle_byr[spc_cat][reg][gsec] = gsec_EF[str(max_below_ref_yr)] + prop_to_ref_yr*(
                                gsec_EF[str(min_above_ref_yr)] - gsec_EF[str(max_below_ref_yr)])
                              
                    # based on Gomez Sanabria et al. (2021): assume CLE EF is const., MFR EF is MFRem/CLEem
                    if "Waste" in gsec and len(ev) > 0:
                        if i_scen == 0: # CLE
                            cle_byr[spc_cat][reg][gsec] = ev[str(max_below_ref_yr)] + prop_to_ref_yr*(
                                ev[str(min_above_ref_yr)] - ev[str(max_below_ref_yr)])
                            gsec_EF = trend_ones
                        else: # MFR
                            gsec_EF = ev # will be divided by CLE later
                    
                    # If we have EF data, add it to our EF dict and exports...
                    if len(gsec_EF) > 0:
                        # add the raw EF (unnormalized) for export 
                        # (using reversed EF vector beacuse MFR only has the latest years of 2050, 2030)
                        EF_EMF_raw.append([scen,spc_cat,reg,inv_sector,fuel,gsec,gsecpctem['2030']] 
                                     + list(reversed(list(gsec_EF))))

                        # nornmalize by cle_byr (if nonzero)
                        if cle_byr[spc_cat][reg][gsec] != 0:
                            gsec_EF = (gsec_EF / cle_byr[spc_cat][reg][gsec])
                        # or if cle_byr is zero (0.1% or less of each spc_cat's emissions), set EF constant
                        else: 
                            gsec_EF = trend_ones
                            
                        # then add to EF sum, with aggregation based on % of CEDS sector emissions in base year
                        EF_EMF[scen][spc_cat][reg][inv_sector][fuel] += gsec_EF * gsecpctem[int_scen_baseyr[i_scen]]
                
                    # If we don't have EF data but there are emissions, assume it's constant and add in proportion to emissions
                    elif len(gsec_EF) == 0 and len(ev) > 0:
                        EF_EMF[scen][spc_cat][reg][inv_sector][fuel] += (trend_ones) * gsecpctem[int_scen_baseyr[i_scen]]

                # if we still don't have data because inv_sector's emissions were zero, set constant at 1
                if (type(EF_EMF[scen][spc_cat][reg][inv_sector][fuel]) == float):
                    EF_EMF[scen][spc_cat][reg][inv_sector][fuel] = (trend_ones)

                # and renormalize after having aggregated from GAINS sector to inventory sector
                cle_byr_EF = (EF_EMF[int_scen[0]][spc_cat][reg][inv_sector][fuel][str(max_below_ref_yr)] + 
                              (EF_EMF[int_scen[0]][spc_cat][reg][inv_sector][fuel][str(min_above_ref_yr)]
                               - EF_EMF[int_scen[0]][spc_cat][reg][inv_sector][fuel][str(max_below_ref_yr)])
                              *prop_to_ref_yr)
                if cle_byr_EF != 0:
                    EF_EMF[scen][spc_cat][reg][inv_sector][fuel] /= cle_byr_EF
                else:
                    EF_EMF[scen][spc_cat][reg][inv_sector][fuel] = trend_ones  
                    
                # finally, if MFR 2030 is greater than CLE 2030, fix to CLE levels (shouldn't be >CLE):
                if scen == int_scen[1] and (EF_EMF[scen][spc_cat][reg][inv_sector][fuel]['2030'] > 
                                         EF_EMF[int_scen[0]][spc_cat][reg][inv_sector][fuel]['2030']):
                    EF_EMF[scen][spc_cat][reg][inv_sector][fuel][int_yrs[scen]] = EF_EMF[int_scen[0]][spc_cat][reg][inv_sector][fuel][int_yrs[scen]]
                
print('seconds taken: ' + str(time.time() - t0))


# %%


#~# now, consolidate GAINS regions (24) to EPPA regions (18)

EF_EMF_agg = {} # the main emissions factor data     
em_agg_reg = {} # the GAINS emissions in an EPPA region, for consolidation
for i_scen, scen in enumerate(int_scen):
    EF_EMF_agg[scen] = {}
    em_agg_reg[scen] = {}
    trend_ones = (EMFem.iloc[2,(i_yr_EMFem[i_scen]):]/
                  EMFem.iloc[2,(i_yr_EMFem[i_scen]):]) # trend of ones for filling constant EF trends
    # for each species category
    for s, spc_cat in enumerate(spc_GAINSEMF):
        EF_EMF_agg[scen][spc_cat] = {}
        em_agg_reg[scen][spc_cat] = {}
        # for each EPPA region in the groupings
        for reg,GAINSregs in reg_map_EPPAGAINSEMF.items(): 
            EF_EMF_agg[scen][spc_cat][reg] = {}
            em_agg_reg[scen][spc_cat][reg] = {}
            inv_sec_last = ''
            for (inv_sector,fuel),GAINS_sectors in sec_map_inv_GAINSfuels.items():
                if inv_sector != inv_sec_last:
                    em_agg_reg[scen][spc_cat][reg][inv_sector] = {}
                    EF_EMF_agg[scen][spc_cat][reg][inv_sector] = {}
                    inv_sec_last = inv_sector
                em_agg_reg[scen][spc_cat][reg][inv_sector][fuel] = 0.0
                EF_EMF_agg[scen][spc_cat][reg][inv_sector][fuel] = 0.0

                # default to the GAINS region if there's only one in the EPPA region
                if len(GAINSregs) < 2:
                    EF_EMF_agg[scen][spc_cat][reg][inv_sector][fuel] = EF_EMF[scen][spc_cat][GAINSregs[0]][inv_sector][fuel]
                    em_agg_reg[scen][spc_cat][reg][inv_sector][fuel] += em_GAINSreg[scen][spc_cat][GAINSreg][inv_sector][fuel]
                    
                # if more than one GAINS region, create a weighted average based on the GAINS region's sectoral emissions                    
                elif len(GAINSregs) > 1:
                    # first, sum the sectoral emissions across the GAINS regions in that EPPA region
                    for GAINSreg in GAINSregs:
                        if type(em_GAINSreg[scen][spc_cat][GAINSreg][inv_sector][fuel]) != float:
                            em_agg_reg[scen][spc_cat][reg][inv_sector][fuel] += em_GAINSreg[scen][spc_cat][GAINSreg][inv_sector][fuel]
                        
                    # then, find the proportion of EPPA-region emissions from that GAINS region, and build the weighted EF 
                    for GAINSreg in GAINSregs:
                        if sum(em_GAINSreg[scen][spc_cat][GAINSreg][inv_sector][fuel]) > 0:
                            # if there are emissions, calculate % from GAINS region in base year (if >0) or max year
                            if em_agg_reg[scen][spc_cat][reg][inv_sector][fuel][int_scen_baseyr[i_scen]] > 0:
                                pct_GAINS_in_EPPA = (em_GAINSreg[scen][spc_cat][GAINSreg][inv_sector][fuel][int_scen_baseyr[i_scen]] / 
                                                     em_agg_reg[scen][spc_cat][reg][inv_sector][fuel][int_scen_baseyr[i_scen]])
                            else: 
                                i_yr_max = pandas.to_numeric(em_agg_reg[scen][spc_cat][reg][inv_sector][fuel]).idxmax()
                                pct_GAINS_in_EPPA = (em_GAINSreg[scen][spc_cat][GAINSreg][inv_sector][fuel][i_yr_max] / 
                                                     em_agg_reg[scen][spc_cat][reg][inv_sector][fuel][i_yr_max])
                            # and then do the aggregation
                            EF_EMF_agg[scen][spc_cat][reg][inv_sector][fuel] += (EF_EMF[scen][spc_cat][GAINSreg][inv_sector][fuel]
                                                                               * pct_GAINS_in_EPPA)
                                
                # finally, if you haven't added a trend yet (due to lack of emissions), set as constant
                if type(EF_EMF_agg[scen][spc_cat][reg][inv_sector][fuel]) == float:
                    EF_EMF_agg[scen][spc_cat][reg][inv_sector][fuel] = trend_ones


# %%


#~# Now add the GAINS NH3 data (denoted by n; too many formatting differences to use a function structure)

t0 = time.time()
print('working on GAINS NH3')
EF_NH3 = {}           # scaled EF dictionary (main output)
em_GAINSreg_NH3 = {}  # emissions by GAINS region and CEDS sector (for EPPA-region weighted average in next code block)
cle_byr_NH3 = {}      # CLE 2030 value to scale
EF_NH3_raw = []       # one-time export: list of raw EF 

# year vector: fix later to be as simple as possible
nyrs = np.sort(n.index.unique())
nstartyrs = [int_yrs[scen][0] for scen in int_scen]
int_yrs_NH3 = {}

# for each emissions scenario (2: CLE, MFR)
for i_scen, scen in enumerate(int_scen):
    EF_NH3[scen] = {}
    em_GAINSreg_NH3[scen] = {} 
    int_yrs_NH3[scen] = nyrs[nyrs>=int(nstartyrs[i_scen])] 
    trend_ones = pandas.Series(np.ones(len(int_yrs_NH3[scen])),index=int_yrs_NH3[scen]) # trend of ones for filling constants
    
    # and set up the dataframe of total FAO production
    FAOi = (FAOtotal[FAOtotal['Domestic production (ktonne)'] == FAOscen[i_scen]].index[0])
    FAOtotals = FAOtotal.loc[(FAOi+1):(FAOi+len(reg_FAO))]
    
    # for each species (just NH3)
    for s, spc_cat in enumerate(spc_GAINSNH3):
        EF_NH3[scen][spc_cat] = {}
        em_GAINSreg_NH3[scen][spc_cat] = {}
        if i_scen == 0:
            cle_byr_NH3[spc_cat] = {} 
        
        # for each GAINS region (will be aggregated to EPPA regions in the next code block) 
        for r,reg in enumerate(reg_GAINSNH3):
            if r%2 == 0:
                print('working on ',i_scen,spc_cat,reg)

            if i_scen == 0:
                cle_byr_NH3[spc_cat][reg] = {} 
            EF_NH3[scen][spc_cat][reg] = {}
            em_GAINSreg_NH3[scen][spc_cat][reg] = {}
            
            # for FAO agricultural scaling
            FAOtotaltrend = FAOtotals[FAOtotals['Domestic production (ktonne)'] == 
                                      reg_map_NH3FAO[reg][0]].iloc[:,1:].squeeze().to_numpy()

            # for each CEDS sector
            inv_sec_last = ''
            for (inv_sector,fuel),GAINS_sectors in sec_map_CEDSGAINSNH3fuels.items():
                # separate treatment for waste and commercial in the next code block
                if inv_sector != 'Waste' and inv_sector != 'Commercial':
                    if inv_sector != inv_sec_last:
                        em_GAINSreg_NH3[scen][spc_cat][reg][inv_sector] = {}
                        EF_NH3[scen][spc_cat][reg][inv_sector] = {}
                        inv_sec_last = inv_sector
                        
                    # first, loop over GAINS sectors in each inventory sector to make a summed emissions trend 
                    em_GAINSreg_NH3[scen][spc_cat][reg][inv_sector][fuel] = 0.0
                    for gsec in GAINS_sectors:    
                        # find the emissions vector (ev)
                        ev = n.loc[((n.IDREGIONS == reg) & (n.ACT_SEC == gsec) & (n.IDSCENARIOS == scen))]['EMISS']
                        # if the emissions data exists,
                        if len(ev) > 0:
                            # fill in any missing values
                            if len(ev) < len(int_yrs_NH3[scen]):
                                ev = ev.reindex(int_yrs_NH3[scen])
                                ev = ev.interpolate(fill_value="extrapolate",limit_direction='both')
                            if sum(ev) > 0:
                                em_GAINSreg_NH3[scen][spc_cat][reg][inv_sector][fuel] += ev

                    # then, loop over GAINS sectors again to calculate EF trends
                    EF_NH3[scen][spc_cat][reg][inv_sector][fuel] = 0.0
                    for gsec in GAINS_sectors:
                        # initialize parameters
                        gsecpctem = trend_ones * 0.0       # % of inv. sector emissions from GAINS sector 
                        if i_scen == 0:
                            cle_byr_NH3[spc_cat][reg][gsec] = 1.0 # CLE EF at base year to normalize MFR
                            
                        # find the emissions and % of CEDS sector emissions in that GAINS region/scenario
                        v = n.loc[((n.IDREGIONS == reg) & (n.ACT_SEC == gsec) & (n.IDSCENARIOS == scen))][['IMPL_EF','EMISS']]
                        if len(v) > 0:
                            # fill any missing year values
                            if len(v) < len(int_yrs_NH3[scen]):
                                v = v.reindex(int_yrs_NH3[scen])
                                v = v.interpolate(fill_value="extrapolate",limit_direction='both')
                            gsecpctem = (v.EMISS / em_GAINSreg_NH3[scen][spc_cat][reg][inv_sector][fuel]).replace(np.inf,0).fillna(0)
                            
                            # grab the CLE base year value (to normalize MFR)
                            if i_scen == 0:
                                cle_byr_NH3[spc_cat][reg][gsec] = v.IMPL_EF[(max_below_ref_yr)] + prop_to_ref_yr*(
                                v.IMPL_EF[(min_above_ref_yr)] - v.IMPL_EF[(max_below_ref_yr)])

                        # If we have EF data, process it...
                        if len(v.IMPL_EF) > 0:
                            # for agriculture, scale by FAO activity trend (relative to total production)
                            if inv_sector == 'Agriculture':
                                FAOtrend = FAO.loc[((FAO.Region == reg_map_NH3FAO[reg][0]) & (FAO.Scenario == FAOscen[i_scen]) 
                                         & (FAO.Item == sec_map_NH3FAO[gsec][0][0]) & (FAO.Element == sec_map_NH3FAO[gsec][0][1])
                                         & (FAO.Units == sec_map_NH3FAO[gsec][0][2]) )]['Value']
                                if len(FAOtrend) > 0:
                                    # normalize (2030, 2050) values by the base year, which differs from FAO years
                                    FAOnorm = (FAOscale * (FAOtrend.iloc[1:] / FAOtrend.iloc[0]) / 
                                               (FAOtotaltrend[1:] / FAOtotaltrend[0]) )
                                    # for CLE, add pre-2030 values of one (since the rest are past years) to fill the trend
                                    if i_scen == 0:
                                        FAOnorm = pandas.concat([trend_ones[trend_ones.index<2030], FAOnorm])
                                    # and then scale the EF trend
                                    v.IMPL_EF *= FAOnorm
                                    
                            # add the raw EF for export (using reversed EF vector beacuse MFR only has the latest years of 2050, 2030)
                            EF_NH3_raw.append([scen,spc_cat,reg,inv_sector,fuel,gsec,gsecpctem[2050]] + list(reversed(list(v.IMPL_EF))))

                            # scale to CLE base year, or set constant if the base year value is zero
                            if cle_byr_NH3[spc_cat][reg][gsec] != 0:
                                v.IMPL_EF = (v.IMPL_EF / cle_byr_NH3[spc_cat][reg][gsec])
                            else:
                                v.IMPL_EF = trend_ones

                            # then add to EF sum based on % of inv_sector emissions in base year
                            EF_NH3[scen][spc_cat][reg][inv_sector][fuel] += v.IMPL_EF * gsecpctem[int(int_scen_baseyr[i_scen])] 
                 
                        # If we don't have EF data but there are emissions, assume it's constant and add in proportion to emissions
                        elif (len(v.IMPL_EF) == 0 or v.IMPL_EF[2030] == 0 or np.isnan(v.IMPL_EF[2030]) ) and len(v) > 0:
                            EF_NH3[scen][spc_cat][reg][inv_sector][fuel] += trend_ones * gsecpctem[int_scen_baseyr[i_scen]] 
                
                    # if we still have no data (because inv_sector's emissions were zero), set EF trend constant at 1
                    if (type(EF_NH3[scen][spc_cat][reg][inv_sector][fuel]) == float):
                        EF_NH3[scen][spc_cat][reg][inv_sector][fuel] = trend_ones

                    # and renormalize after having aggregated from GAINS sector to inventory sector
                    cle_byr_EF = (EF_NH3[int_scen[0]][spc_cat][reg][inv_sector][fuel][max_below_ref_yr] + 
                                  (EF_NH3[int_scen[0]][spc_cat][reg][inv_sector][fuel][min_above_ref_yr]
                                   - EF_NH3[int_scen[0]][spc_cat][reg][inv_sector][fuel][max_below_ref_yr])
                                  *prop_to_ref_yr)

                    if cle_byr_EF != 0:
                        EF_NH3[scen][spc_cat][reg][inv_sector][fuel] /= cle_byr_EF
                    else:
                        EF_NH3[scen][spc_cat][reg][inv_sector][fuel] = trend_ones 
                    
                    # and renormalize agriculture again because of FAO's different base year
                    if inv_sector == 'Agriculture': 
                        EF_NH3[scen][spc_cat][reg][inv_sector][fuel] /= cle_byr_EF 
                    
                    # and finally, if MFR 2030 is greater than CLE 2030, fix to CLE levels (shouldn't be >CLE):
                    if scen == int_scen[1] and (EF_NH3[scen][spc_cat][reg][inv_sector][fuel][2030] > 
                                             EF_NH3[int_scen[0]][spc_cat][reg][inv_sector][fuel][2030]):
                        EF_NH3[scen][spc_cat][reg][inv_sector][fuel] = EF_NH3[int_scen[0]][spc_cat][reg][inv_sector][fuel][int_yrs_NH3[scen]]

print('seconds taken: ' + str(time.time() - t0))


# %%


#~# and then move from GAINS G20 regions to EPPA regions

t0 = time.time()
EF_NH3_agg = {} # the main emissions factor data  
em_agg_reg_NH3 = {}   # the EPPA-sectoral GAINS emissions by EPPA region, for consolidation in the sectoral global figure 
for i_scen, scen in enumerate(int_scen):
    EF_NH3_agg[scen] = {}
    em_agg_reg_NH3[scen] = {}
    trend_ones = pandas.Series(np.ones(len(int_yrs[scen])),index=int_yrs[scen]) # vector of trend_ones for filling constants
    # for each species
    for s, spc_cat in enumerate(spc_GAINSNH3):
        EF_NH3_agg[scen][spc_cat] = {}
        em_agg_reg_NH3[scen][spc_cat] = {}
        # for each EPPA region in the groupings
        for reg,GAINSregs in reg_map_EPPAGAINSNH3.items(): 
            EF_NH3_agg[scen][spc_cat][reg] = {}
            em_agg_reg_NH3[scen][spc_cat][reg] = {}
            inv_sec_last = ''
            for (inv_sector,fuel),GAINS_sectors in sec_map_CEDSGAINSNH3fuels.items():
                if inv_sector != inv_sec_last:
                    em_agg_reg_NH3[scen][spc_cat][reg][inv_sector] = {}
                    EF_NH3_agg[scen][spc_cat][reg][inv_sector] = {}
                    inv_sec_last = inv_sector
                em_agg_reg_NH3[scen][spc_cat][reg][inv_sector][fuel] = 0.0
                EF_NH3_agg[scen][spc_cat][reg][inv_sector][fuel] = 0.0
                                
                # for the non-G20 countries, assume CLE is constant for non-ag, unless CLE of G20 is increasing
                if (reg not in reg_G20) and (i_scen == 0) and (inv_sector != 'Agriculture'):
                    if (EF_NH3[scen][spc_cat][GAINSreg][inv_sector][fuel].iloc[-1] <= 
                        EF_NH3[scen][spc_cat][GAINSreg][inv_sector][fuel].iloc[0] for GAINSreg in GAINSregs):
                        EF_NH3_agg[scen][spc_cat][reg][inv_sector][fuel] = trend_ones
                    continue

                # separate treatment for waste: mirror NOx from Gomez Sanabria et al. (2021) due to large NH3 data gaps
                # this will be constant for CLE and trending to zero for MFR (with differing trends by region)
                if inv_sector == 'Waste':
                    EF_NH3_agg[scen][spc_cat][reg][inv_sector][fuel] = EF_EMF_agg[scen]['NOX'][reg][inv_sector][fuel]
                    continue
                    
                # special treatment for commercial, whose sectors weren't distinct from residential in the GAINS NH3 input
                # mirror residential, but you'll have to do it below so residential can be read in first
                if inv_sector == 'Commercial':
                    continue 
                    
                else: 
                    # default to the GAINS region if there's only one in the EPPA region
                    if len(GAINSregs) < 2:
                        EF_NH3_agg[scen][spc_cat][reg][inv_sector][fuel] += EF_NH3[scen][spc_cat][GAINSregs[0]][inv_sector][fuel]
                        em_agg_reg_NH3[scen][spc_cat][reg][inv_sector][fuel] += em_GAINSreg_NH3[scen][spc_cat][GAINSregs[0]][inv_sector][fuel]

                    # if more than one GAINS region, use a weighted average based on the GAINS region's % of EPPA-region emissions                    
                    elif len(GAINSregs) > 1:
                        # first, sum the sectoral emissions across the GAINS regions in that EPPA region
                        for GAINSreg in GAINSregs:
                            if type(em_GAINSreg_NH3[scen][spc_cat][GAINSreg][inv_sector][fuel]) != float:
                                em_agg_reg_NH3[scen][spc_cat][reg][inv_sector][fuel] += em_GAINSreg_NH3[scen][spc_cat][GAINSreg][inv_sector][fuel]

                        # then, find the proportion of EPPA-region emissions from that GAINS region, and build the weighted EF 
                        EF_NH3_agg[scen][spc_cat][reg][inv_sector][fuel] = 0.0
                        for GAINSreg in GAINSregs:
                            # (if there are emissions at all, that is...otherwise don't add)
                            if type(em_GAINSreg_NH3[scen][spc_cat][GAINSreg][inv_sector][fuel]) != float:
                                if sum(em_agg_reg_NH3[scen][spc_cat][reg][inv_sector][fuel]) > 0:
                                    pct_GAINS_in_EPPA = (em_GAINSreg_NH3[scen][spc_cat][GAINSreg][inv_sector][fuel][int(int_scen_baseyr[i_scen])] / 
                                                         em_agg_reg_NH3[scen][spc_cat][reg][inv_sector][fuel][int(int_scen_baseyr[i_scen])]) 
                                    EF_NH3_agg[scen][spc_cat][reg][inv_sector][fuel] += (EF_NH3[scen][spc_cat][GAINSreg][inv_sector][fuel] 
                                                                                        * pct_GAINS_in_EPPA)
                    
                    # and if there were no emissions in the EPPA region, make the intensity vector constant
                    if type(EF_NH3_agg[scen][spc_cat][reg][inv_sector][fuel]) == float:
                        EF_NH3_agg[scen][spc_cat][reg][inv_sector][fuel] = trend_ones
                     
            # finally, add commercial after residential given the above note
            for (inv_sector,fuel),GAINS_sectors in sec_map_CEDSGAINSNH3fuels.items():
                if inv_sector == 'Commercial':
                    EF_NH3_agg[scen][spc_cat][reg][inv_sector][fuel] = EF_NH3_agg[scen][spc_cat][reg]['Residential'][fuel]


                    


# %%


#~# 2(c): Create intensity scenarios ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# in this implementation, we use full-century exponential fits to the 2000-2050 GAINS trends.
# different scenarios could be created for different research questions. 

# package documentation: https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.curve_fit.html
# from least-squares: https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.leastsq.html#scipy.optimize.leastsq

print('working on intensity scenarios')
t0 = time.time()

import scipy.optimize
def monoExp(x, m, gt, b):
    return m * np.exp(-gt * x) + b

#x values for exponential (where base year = 0). 
gbyr = int(int_yrs[int_scen[0]][0]) # GAINS initial year
# full vector for CLE (2000,05,10,20,30,50); past years + MFR years for MFR (2000,05,10,15,30,50)
xgcle = np.array(int_yrs[int_scen[0]]).astype(float)-gbyr
xgmfr = np.array(int_yrs[int_scen[0]][0:4]+int_yrs[int_scen[1]]).astype(float)-gbyr
# and append the 2014 = 1 value
xgcle = np.append(xgcle, (float(yr_list[0])-gbyr))
xgmfr = np.append(xgmfr, (float(yr_list[0])-gbyr))
xg = [xgcle,xgmfr]

efs_to_fit = 0.0
gainsex = {}  # dict for the exp fit parameters
gexfits = []  # export for the fits

# prepare the EMF vs NH3 inputs so you can assign the right dictionaries
EF_EMF_aggdict = {spc_cat: 'EF_EMF_agg' for spc_cat in spc_GAINSEMF}
EF_EMF_aggdict['NH3'] = 'EF_NH3_agg' 
fuelmapdict = {spc_cat: 'sec_map_inv_GAINSfuels' for spc_cat in spc_GAINSEMF}
fuelmapdict['NH3'] = 'sec_map_CEDSGAINSNH3fuels'

for i_scen, scen in enumerate(int_scen):
    gainsex[scen] = {}
    for s, spc_cat in enumerate(spc_GAINSEMF+spc_GAINSNH3):
        gainsex[scen][spc_cat] = {}
        print('working on ',spc_cat)
        for reg,GAINSregs in reg_map_EPPAGAINSEMF.items(): 
            gainsex[scen][spc_cat][reg] = {}
            inv_sec_last = ''
            for (inv_sector,fuel) in locals()[fuelmapdict[spc_cat]].keys():
                if inv_sector != inv_sec_last:
                    gainsex[scen][spc_cat][reg][inv_sector] = {}
                    inv_sec_last = inv_sector

                # read in points: full vector for CLE (2000,05,10,15,20,30,50)                    
                if i_scen == 0:
                    efs_to_fit = locals()[EF_EMF_aggdict[spc_cat]][scen][spc_cat][reg][inv_sector][fuel] 
                # modified for MFR (CLE 2000,05,10,15 plus MFR30,50)
                elif i_scen > 0:
                    efs_to_fit = locals()[EF_EMF_aggdict[spc_cat]][int_scen[0]][spc_cat][reg][inv_sector][fuel][0:4].append(
                        locals()[EF_EMF_aggdict[spc_cat]][scen][spc_cat][reg][inv_sector][fuel])

                # (for re-runs during testing, we had to start fresh by removing the ref_yr additions)
                if ref_yr in efs_to_fit.index:
                    efs_to_fit = efs_to_fit.drop(labels=[ref_yr])

                # and in both cases, append the ref_yr point = 1 
                efs_to_fit[str(ref_yr)] = 1
                efs_to_fit.index = efs_to_fit.index.astype(int) 
                # (with sigma weighting specified so the curve goes through it)
                sigma = np.ones(len(efs_to_fit))
                sigma[efs_to_fit.index.get_loc(ref_yr)] = 0.01
                # and add similar weighting for MFR waste at 2050 to follow Gomez Sanabria et al. (2021)
                if scen == int_scen[-1] and inv_sector == 'Waste':
                    sigma[efs_to_fit.index.get_loc(2050)] = 0.01

                # if it's a flat line at 1, don't fit an exponential 
                if (round(max(efs_to_fit),5) == 1.0 and round(min(efs_to_fit),5) == 1.0):
                    gainsex[scen][spc_cat][reg][inv_sector][fuel] = [1.0,1.0,1.0]

                else:
                    yg = efs_to_fit 
                    # perform the fit
                    params, cv = scipy.optimize.curve_fit(monoExp, xg[i_scen], yg, sigma=sigma,
                                                          bounds=([-np.inf,-np.inf,0], np.inf),maxfev=10000)
                    m, gt, b = params
                    gainsex[scen][spc_cat][reg][inv_sector][fuel] = params

                    # calculate and store r^2, alongside the CLE+MFR years (not 2020 which is only in CLE) of pre-fit data
                    squaredDiffs = np.square(yg - monoExp(xg[i_scen], m, gt, b))
                    squaredDiffsFromMean = np.square(yg - np.mean(yg))
                    rsquared = 1 - np.sum(squaredDiffs) / np.sum(squaredDiffsFromMean)
                    gexfits.append([scen,spc_cat,reg,inv_sector,fuel,rsquared,m,gt,b]+
                                  list(efs_to_fit.sort_index()[:5])+list(efs_to_fit.sort_index()[-2:]))   

print('seconds taken: ' + str(time.time() - t0))


# %%


#~# 3: Outputs

#~# 3(a): Export intensity outputs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# export the "raw EFs" by GAINS sector and region
def export_GAINS_raw_EFs(raw_EF_list,filename):
    column_headers = ['Scenario', 'Species','Region','Inventory Sector','CEDS Fuel','GAINS Sector','Sectoral Emissions Fraction']
    raw_EF_df = pandas.DataFrame(raw_EF_list, columns = column_headers+list(reversed(int_yrs[int_scen[0]]))) 
    # reorder to "un-reverse" the years
    raw_EF_df = raw_EF_df[column_headers+int_yrs[int_scen[0]]]
    # export to csv
    raw_EF_path = os.path.join(output_dir,filename)
    raw_EF_df.to_csv(raw_EF_path, index = False, header=True)
    return
# and complete the exports for EMF and NH3 
export_GAINS_raw_EFs(EF_EMF_raw,'GAINSEMFEFs.csv')
export_GAINS_raw_EFs(EF_NH3_raw,'GAINSNH3EFs.csv')

# export the GAINS EF trends by inventory sector and EPPA region
ef_eppa_reg = []
for i_scen, scen in enumerate(int_scen):
    for s, spc_cat in enumerate(spc_GAINSEMF):
        for reg,GAINSregs in reg_map_EPPAGAINSEMF.items(): 
            for (inv_sector,fuel),GAINS_sectors in sec_map_inv_GAINSfuels.items():
                ef_eppa_reg.append([scen,spc_cat,reg,inv_sector,fuel] + list(reversed(list(EF_EMF_agg[scen][spc_cat][reg][inv_sector][fuel]))))
for i_scen, scen in enumerate(int_scen):
    for s, spc_cat in enumerate(spc_GAINSNH3):
        for reg,GAINSregs in reg_map_EPPAGAINSNH3.items(): 
            for (inv_sector,fuel),GAINS_sectors in sec_map_CEDSGAINSNH3fuels.items():
                ef_eppa_reg.append([scen,spc_cat,reg,inv_sector,fuel] + list(reversed(list(EF_NH3_agg[scen][spc_cat][reg][inv_sector][fuel]))))
# make the dataframe to export to csv 
column_headers = ['Scenario', 'Species','Region','Inventory Sector','CEDS Fuel']
GAINS_then_ref_yr = list(sorted(int_yrs[int_scen[0]])+[str(ref_yr)]) 
ef_eppa_reg_df = pandas.DataFrame(ef_eppa_reg, columns = column_headers+list(reversed(GAINS_then_ref_yr))) 
# reorder to "un-reverse" the years and remove the ref_yr value of 1 that was added during the fits
ef_eppa_reg_df = ef_eppa_reg_df[column_headers+int_yrs[int_scen[0]]]
# export to csv
ef_eppa_reg_path = os.path.join(output_dir,'GAINS_EFs_by_EPPA_region.csv')
ef_eppa_reg_df.to_csv(ef_eppa_reg_path, index = False, header=True)

# export the fits, including the ref_yr since it was part of the fit
GAINS_and_ref_yr = list(sorted(int_yrs[int_scen[0]]+[str(ref_yr)]))
GAINS_and_ref_yr.remove('2020') # not included in fits since it wasn't in the MFR time series
gexfitspath = os.path.join(output_dir,'GAINSEF_Fits.csv')
gexfitsdf = pandas.DataFrame(gexfits, columns = column_headers+['R squared','m','gt','b']+GAINS_and_ref_yr) 
gexfitsdf.to_csv(gexfitspath, index = False, header=True)


# %%


#~# 3(b): Export scaling of CEDS inventory ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# years after initial GAINS year, used for implementing exponential fits
yrs_from_GAINS_start = np.array(yr_list) - int(int_yrs[int_scen[0]][0])     

t0 = time.time()
CEDS_scaling = {}        # the main scaling dict   
components_export = []   # to export components for evaluation: inventory, plus activity scaling and intensity scaling (2030, 2050, 2100) 
CEDS_scaling_export = [] # to export the main scaling trends
fpct_inGAINS = {}        # proportion (of sectoral emissions) covered by GAINS fuels

print('exporting CEDS scaling')

for i_scen, scen in enumerate(act_scen):
    CEDS_scaling[scen] = {}
    for s, (spc_cat,spcs_CEDS) in enumerate(spc_map_CEDS.items()):
        for spc_CEDS in spcs_CEDS:
            fpct_inGAINS[spc_CEDS] = {} 
            CEDS_scaling[scen][spc_CEDS] = {} 
            for reg in region_names:
                fpct_inGAINS[spc_CEDS][reg] = {}  
                CEDS_scaling[scen][spc_CEDS][reg] = {}

                # first, find the proportion of each sector's base year emissions covered by GAINS-scaled fuels
                # if this is less than 1, we will divide existing fuel trends by that to cover the gap
                inv_sec_last = ''
                for (inv_sector,fuel) in locals()[fuelmapdict[spc_map_GAINS[spc_cat]]].keys():
                    if inv_sector not in GFED_sec_dict.keys():
                        if inv_sector != inv_sec_last:
                            inv_sec_last = inv_sector
                            fpct_inGAINS[spc_CEDS][reg][inv_sector] = 0.0
                            sec = spc_CEDS + '_' + inv_sec_map[inv_sector][0]   
                        fpct_inGAINS[spc_CEDS][reg][inv_sector] += fpct[spc_CEDS][sec][reg][fuel]

                # then, loop through and scale GAINS intensity factors * EPPA activity factors
                for i_int_scen in int_scen:
                    CEDS_scaling[scen][spc_CEDS][reg][i_int_scen] = {}
                    inv_sec_last = ''
                    for (inv_sector,fuel) in locals()[fuelmapdict[spc_map_GAINS[spc_cat]]].keys():
                        if inv_sector not in GFED_sec_dict.keys():
                            # sum of emissions across fuels: to find inventory value for "components" export
                            if inv_sector != inv_sec_last:
                                fuelsum = 0.0
                                inv_sec_last = inv_sector
                                sec = spc_CEDS + '_' + inv_sec_map[inv_sector][0]
                                CEDS_scaling[scen][spc_CEDS][reg][i_int_scen][inv_sector] = {} 
                                for f in fuels:
                                    fuelsum += em_CEDS[spc_CEDS][f][sec][reg]                            

                            # extract the GAINS-based emissions intensity trends (int_scaling)
                            # use NH3 trends for agricultural sectors since GAINS only has ag. in NH3
                            if 'Agricult' in inv_sector: 
                                spc_here = 'NH3'
                            else: 
                                spc_here = spc_cat
                                
                            # calculate the intensity scaling trend
                            m,gt,b = (gainsex[i_int_scen][spc_map_GAINS[spc_here]][reg][inv_sector][fuel][0], 
                                      gainsex[i_int_scen][spc_map_GAINS[spc_here]][reg][inv_sector][fuel][1], 
                                      gainsex[i_int_scen][spc_map_GAINS[spc_here]][reg][inv_sector][fuel][2])
                            int_scaling = monoExp(yrs_from_GAINS_start, m,gt,b)

                            # if the fuels in GAINS cover this fuel, scale
                            if fpct_inGAINS[spc_CEDS][reg][inv_sector] > 0:
                                scaling = (int_scaling 
                                          * act_scaling[scen][inv_sector][reg][fuels.index(fuel)] 
                                          *  (1/fpct_inGAINS[spc_CEDS][reg][inv_sector]) # sector-fuel coverage adjustment
                                          )
                                CEDS_scaling[scen][spc_CEDS][reg][i_int_scen][inv_sector][fuel] = scaling
                                
                                # export components: inventory, plus activity scaling and intensity scaling (2030, 2050, 2100)
                                components_export.append([scen,spc_cat+'_'+spc_CEDS,reg,i_int_scen,inv_sector,fuel,
                                                fuelsum  * (1/fpct_inGAINS[spc_CEDS][reg][inv_sector]),
                                                act_scaling[scen][inv_sector][reg][fuels.index(fuel)][yr_list.index(2030)],
                                                act_scaling[scen][inv_sector][reg][fuels.index(fuel)][yr_list.index(2050)],
                                                act_scaling[scen][inv_sector][reg][fuels.index(fuel)][yr_list.index(2100)],
                                                int_scaling[yr_list.index(2030)], 
                                                    int_scaling[yr_list.index(2050)], 
                                                    int_scaling[yr_list.index(2100)]
                                               ])

                            # and if there's no scaling, add zero
                            else:
                                CEDS_scaling[scen][spc_CEDS][reg][i_int_scen][inv_sector][fuel] = np.zeros(len(int_scaling))
                                
                            # and finally, export the scaling
                            CEDS_scaling_export.append([scen,i_int_scen,spc_CEDS,reg,inv_sector,fuel]+
                                                       list(CEDS_scaling[scen][spc_CEDS][reg][i_int_scen][inv_sector][fuel]))

CEDS_scaling_path = os.path.join(output_dir,'CEDS_scaling.csv')
CEDS_scaling_df = pandas.DataFrame(CEDS_scaling_export, 
                                   columns = ['Activity_Scenario','Intensity_Scenario','Species',
                                              'Region','Sector','Fuel']+yr_list)
CEDS_scaling_df.to_csv(CEDS_scaling_path, index = False, header=True)

print('seconds taken: ' + str(time.time() - t0))


# %%


#~# 3(c): Export scaling of GFED inventory ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

GFED_scaling = {}        # dict for the scaling trends
GFED_scaling_export = [] # export for the scaling trends
fuel = 'process'         # all GFED sources are defined as "process"

print('exporting GFED scaling')

for i_scen, scen in enumerate(act_scen):
    GFED_scaling[scen] = {}
    for s, (spc_cat,spcs_GFED) in enumerate(spc_map_GFED.items()):
        for spc_GFED in spcs_GFED:
            GFED_scaling[scen][spc_GFED] = {}  
            for reg in region_names:
                GFED_scaling[scen][spc_GFED][reg] = {}
                for i_int_scen in int_scen:
                    GFED_scaling[scen][spc_GFED][reg][i_int_scen] = {}
                    for (inv_sector, fuel) in locals()[fuelmapdict[spc_map_GAINS[spc_cat]]].keys():
                        if inv_sector in GFED_sec_dict.keys():
                            GFED_scaling[scen][spc_GFED][reg][i_int_scen][inv_sector] = {} 
                            
                            # extract the GAINS-based intensity trends
                            # use GAINS NH3 trends for agricultural sectors since GAINS only has ag. in NH3
                            if 'Agricult' in inv_sector: 
                                spc_here = 'NH3'
                            m,gt,b = (gainsex[i_int_scen][spc_map_GAINS[spc_here]][reg][inv_sector][fuel][0], 
                                        gainsex[i_int_scen][spc_map_GAINS[spc_here]][reg][inv_sector][fuel][1], 
                                        gainsex[i_int_scen][spc_map_GAINS[spc_here]][reg][inv_sector][fuel][2])
                            int_scaling = monoExp(yrs_from_GAINS_start, m,gt,b)
                            scaling = (int_scaling *  act_scaling[scen][inv_sector][reg][fuels.index(fuel)])
                            GFED_scaling[scen][spc_GFED][reg][i_int_scen][inv_sector][fuel] = scaling
                            
                            # export components: inventory, plus activity scaling and intensity scaling (2030, 2050, 2100)
                            components_export.append([scen,spc_cat+'_'+spc_GFED,reg,i_int_scen,inv_sector,fuel,
                                                em_GFED[inv_sector][spc_GFED][reg],
                                                act_scaling[scen][inv_sector][reg][fuels.index(fuel)][yr_list.index(2030)],
                                                act_scaling[scen][inv_sector][reg][fuels.index(fuel)][yr_list.index(2050)],
                                                act_scaling[scen][inv_sector][reg][fuels.index(fuel)][yr_list.index(2100)],
                                                int_scaling[yr_list.index(2030)], 
                                                    int_scaling[yr_list.index(2050)], 
                                                    int_scaling[yr_list.index(2100)]
                                               ])
                            
                            # and export scaling
                            GFED_scaling_export.append([scen,i_int_scen,spc_GFED,reg,inv_sector,fuel]+
                                                       list(GFED_scaling[scen][spc_GFED][reg][i_int_scen][inv_sector][fuel]))

GFED_scaling_path = os.path.join(output_dir,'GFED_scaling.csv')
GFED_scaling_df = pandas.DataFrame(GFED_scaling_export, 
                                   columns = ['Activity_Scenario','Intensity_Scenario','Species',
                                              'Region','Sector','Fuel']+yr_list)
GFED_scaling_df.to_csv(GFED_scaling_path, index = False, header=True)

# and export the components
componentspath = os.path.join(output_dir,'Emissions_components.csv')
componentsdf = pandas.DataFrame(components_export, columns = [['Activity Scenario', 'Species','Region','Intensity Scenario',
                                                     'Inventory Sector','CEDS Fuel','Inventory_2014',
                                                     'Activity_2030','Activity_2050','Activity_2100',
                                                     'Intensity_2030','Intensity_2050','Intensity_2100']]) 
componentsdf.to_csv(componentspath, index = False, header=True)

print('exports complete')
# END
