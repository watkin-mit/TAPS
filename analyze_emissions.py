#!/usr/bin/env python
# coding: utf-8
# %%

# %%
'''
This script converts inventory and scaling inputs (from scale_emissions.py) to emissions trends.
It imports SSP emissions trends for comparison, and analyzes them in figure exports as submitted to GMD. 

Before running (for Fig. 3 and Fig. 4): 
Download stackedBarGraph.py into this directory from the following link: https://github.com/minillinim/stackedBarGraph

'''
# import packages 
import numpy as np
import os
import io
import matplotlib.pyplot as plt
import time
import pandas

input_dir = 'input_files'
output_dir = 'scaling_output'


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
    return nested_dict

print('reading in the scaling files')

# read in CEDS scaling file 
CEDS_scaling_file = pandas.read_csv(os.path.join(output_dir,'CEDS_scaling.csv'))
CEDS_scaling = df_to_dict(CEDS_scaling_file)

# read in GFED scaling file
GFED_scaling_file = pandas.read_csv(os.path.join(output_dir,'GFED_scaling.csv'))
GFED_scaling = df_to_dict(GFED_scaling_file)

# read in CEDS inventory data
em_CEDS_file = pandas.read_csv(os.path.join(output_dir,'em_CEDS.csv'))
em_CEDS = df_to_dict(em_CEDS_file)

# read in GFED inventory data
em_GFED_file = pandas.read_csv(os.path.join(output_dir,'em_GFED.csv'))
em_GFED = df_to_dict(em_GFED_file)


# %%


# define the parameters needed for analysis and figures

# year list: combination of inventory year and scaling years, avoiding duplicates if the reference year is in the scaling years
ref_yr = [int(yr_col) for yr_col in em_CEDS_file.columns if str.isdigit(yr_col)]
scaling_yr = [int(yr_col) for yr_col in CEDS_scaling_file.columns if str.isdigit(yr_col)]
yr_list = sorted(list(set(ref_yr) | set(scaling_yr)))

# species mapping
spc_map_csv = pandas.read_csv(os.path.join(input_dir,'spc_map.csv'))
spc_map_CEDS = spc_map_csv.groupby('EPPA')['CEDS'].apply(list).to_dict() 
spc_map_CEDS['VOC'] = set(spc_map_CEDS['VOC']) # remove duplicates from VOC speciation
spc_map_GFED = spc_map_csv.dropna().groupby('EPPA')['GFED'].apply(list).to_dict() 
spc_map_GAINS = dict(zip(spc_map_csv.dropna().EPPA,spc_map_csv.dropna().GAINS))

# sectoral mapping
inv_sec_map_file = pandas.read_csv(os.path.join(input_dir,'sectoral_mapping_EPPA7_inventories.csv'))

# map inventory sector names to codes
inv_sec_map = inv_sec_map_file.groupby('Inventory_Name')['Inventory_Code'].apply(list).to_dict()
for key, values in inv_sec_map.items():
    inv_sec_map[key] = list(set(inv_sec_map[key]))
    
# map inventory sector names to SSP sector names
sec_map_SSP = inv_sec_map_file.groupby('SSP')['Inventory_Name'].apply(list).to_dict()
for key, values in sec_map_SSP.items():
    sec_map_SSP[key] = list(set(sec_map_SSP[key]))
    
# map full SSP sector names to shortened names for figure
sec_names_SSP = dict(zip(inv_sec_map_file.SSP,inv_sec_map_file.Figure))

# map inventory sectors to GAINS sectors
sec_map_CEDSGAINSEMF_file = pandas.read_csv(os.path.join(input_dir,'CEDS_GAINSEMF_sectoral_mapping.csv'))
sec_map_CEDSGAINSEMFfuels = sec_map_CEDSGAINSEMF_file.groupby(['CEDS2020','FuelCEDS'])['GAINS_EMF'].apply(list).to_dict()
sec_map_CEDSGAINSNH3_file = pandas.read_csv(os.path.join(input_dir,'CEDS_GAINSNH3_sectoral_mapping.csv'))
sec_map_CEDSGAINSNH3fuels = sec_map_CEDSGAINSNH3_file.groupby(['CEDS','CEDSFuel'])['GAINS'].apply(list).to_dict() 
# point to the correct fuel mapping depending on the GAINS species (EMF or NH3)
spc_GAINSEMF = [spc for spc in spc_map_csv.dropna().GAINS.unique() if spc != 'NH3']
fuelmapdict = {spc_cat: 'sec_map_CEDSGAINSEMFfuels' for spc_cat in spc_GAINSEMF}
fuelmapdict['NH3'] = 'sec_map_CEDSGAINSNH3fuels'

# EPPA scenarios and SSP (CMIP6, IAM) connections
scen_map_EPPA_file = pandas.read_csv(os.path.join(input_dir,'scen_map_EPPA.csv'))
scen_map_EPPAfig = dict(zip(scen_map_EPPA_file.EPPA_energy,scen_map_EPPA_file.Figure))
scen_map_EPPA_CMIP6 = dict(zip(scen_map_EPPA_file.EPPA_energy,scen_map_EPPA_file.CMIP6))
scen_map_EPPA_IAM = dict(zip(scen_map_EPPA_file.EPPA_energy,scen_map_EPPA_file.IAM))

# general lists from the scaling files
list_act_scen = CEDS_scaling_file.Activity_Scenario.unique()
list_reg = CEDS_scaling_file.Region.unique()
list_int_scen = CEDS_scaling_file.Intensity_Scenario.unique()
list_sec_CEDS = CEDS_scaling_file.Sector.unique()
list_fuel = CEDS_scaling_file.Fuel.unique()
list_sec_GFED = GFED_scaling_file.Sector.unique() 
list_sec = sorted(list(set(list_sec_CEDS) | set(list_sec_GFED)))


# %%


# calculate emissions scenario trends in Tg yr^-1 (from inventory and scaling)

scaleNOx = 46/30 # from molecular weights of NO (CEDS data) vs NO2 (to follow McDuffie et al., 2020)
em_scen = {} 
em_scen_export = []

for i_act_scen, act_scen in enumerate(list_act_scen):
    em_scen[act_scen] = {}
    for int_scen in list_int_scen:
        em_scen[act_scen][int_scen] = {}
        # just the 7 species categores. VOC species are aggregated below for each inventory
        for spc in spc_map_CEDS.keys():
            em_scen[act_scen][int_scen][spc] = {}        
            for reg in list_reg:
                em_scen[act_scen][int_scen][spc][reg] = {}
                
                # first, add GFED 
                for sec_name in list_sec_GFED:
                    em_scen[act_scen][int_scen][spc][reg][sec_name] = 0.0
                    for spc_GFED in spc_map_GFED[spc]:
                        em_scen[act_scen][int_scen][spc][reg][sec_name] += (em_GFED[spc_GFED][reg][sec_name] * 
                                                                   GFED_scaling[act_scen][int_scen][spc_GFED][reg][sec_name]['process'])
                        
                # then, add CEDS
                for sec_name in list_sec_CEDS:
                    em_scen[act_scen][int_scen][spc][reg][sec_name] = 0.0
                    for spc_CEDS in spc_map_CEDS[spc]:
                        # full sector name in base-year data
                        sec = spc_CEDS + '_' + inv_sec_map[sec_name][0] # dict gives sector code
                        # add for each fuel, if it has a scaling from GAINS 
                        for fuel in list_fuel:
                            if (sec_name, fuel) in locals()[fuelmapdict[spc_map_GAINS[spc]]].keys():
                                em_scen[act_scen][int_scen][spc][reg][sec_name] += (em_CEDS[spc_CEDS][fuel][sec][reg] 
                                                                                    * CEDS_scaling[act_scen][int_scen][spc_CEDS][reg][sec_name][fuel])                        

                    # convert molecular weights from NO to NO2 
                    if spc == 'NOx':
                        em_scen[act_scen][int_scen][spc][reg][sec_name] *= scaleNOx

                    # and append the values to export into csv
                    em_scen_export.append([act_scen,int_scen,spc,reg,sec_name] + list(em_scen[act_scen][int_scen][spc][reg][sec_name]))
                    
em_scen_path = os.path.join(output_dir,'Emissions_scenarios.csv')
em_scendf = pandas.DataFrame(em_scen_export, columns = [['Activity_Scenario', 'Intensity_Scenario','Species','Region','Sector']+yr_list]) 
em_scendf.to_csv(em_scen_path, index = False, header=True)


# %%


# processing emissions time series for each figure (different aggregations)

# for Figure 4: aggregate the sectors, keeping regions
em_scen_byreg = {}
for i_act_scen, act_scen in enumerate(list_act_scen):
    em_scen_byreg[act_scen] = {}
    for int_scen in list_int_scen:            
        em_scen_byreg[act_scen][int_scen] = {}
        for s, spc in enumerate(spc_map_CEDS.keys()):
            em_scen_byreg[act_scen][int_scen][spc] = {}    
            for reg in list_reg:
                em_scen_byreg[act_scen][int_scen][spc][reg] = 0.0
                for sec_name in list_sec:
                    em_scen_byreg[act_scen][int_scen][spc][reg] += em_scen[act_scen][int_scen][spc][reg][sec_name]
                                        
# for Figure 2 and 3: aggregate the regions, keeping sectors 
em_scen_bysec = {} 
for i_act_scen, act_scen in enumerate(list_act_scen):
    em_scen_bysec[act_scen] = {}
    for int_scen in list_int_scen:            
        em_scen_bysec[act_scen][int_scen] = {}
        for s, spc in enumerate(spc_map_CEDS.keys()):
            em_scen_bysec[act_scen][int_scen][spc] = {}    
            for sec_name in list_sec: 
                em_scen_bysec[act_scen][int_scen][spc][sec_name] = 0.0
                for reg in list_reg:
                    em_scen_bysec[act_scen][int_scen][spc][sec_name] += em_scen[act_scen][int_scen][spc][reg][sec_name]

# For Figure 3, aggregate from 12 CEDS_GBD-MAPS and GFED to 9 SSP (CMIP) sectors 
em_scen_bySSPsec = {} 
for i_act_scen, act_scen in enumerate(list_act_scen):
    em_scen_bySSPsec[act_scen] = {}
    for int_scen in list_int_scen: 
        em_scen_bySSPsec[act_scen][int_scen] = {}
        for s, spc in enumerate(spc_map_CEDS.keys()):
            em_scen_bySSPsec[act_scen][int_scen][spc] = {}  
            for SSP_sector, inv_sector in sec_map_SSP.items():
                em_scen_bySSPsec[act_scen][int_scen][spc][SSP_sector] = 0.0
                # if the mapping is 1:1, keep it the same
                if len(inv_sector) < 2:
                    em_scen_bySSPsec[act_scen][int_scen][spc][SSP_sector] = em_scen_bysec[act_scen][int_scen][spc][inv_sector[0]]

                # if there's more than 1 inv_sector for the SSP sector, take the sum 
                if len(inv_sector) > 1:
                    for csec in inv_sector:
                        em_scen_bySSPsec[act_scen][int_scen][spc][SSP_sector] += em_scen_bysec[act_scen][int_scen][spc][csec]

# for Figure 2, aggregate the sectors to global totals 
em_scen_totals = {}
for i_act_scen, act_scen in enumerate(list_act_scen):
    em_scen_totals[act_scen] = {}  
    for int_scen in list_int_scen:
        em_scen_totals[act_scen][int_scen] = {}
        for s, spc in enumerate(spc_map_CEDS.keys()):
            em_scen_totals[act_scen][int_scen][spc] = 0.0  
            for SSP_sector, inv_sector in sec_map_SSP.items():
                em_scen_totals[act_scen][int_scen][spc] += em_scen_bySSPsec[act_scen][int_scen][spc][SSP_sector]


# %%


# prepare the SSP comparison data: CMIP6 scenarios for Figure 3-4 

# read the file
CMIP6_file = 'SSP_CMIP6_201811_solvents.csv' # new filename after shortening the "Solvents" sector name for processing
CMIP6 = pandas.read_csv(os.path.join(input_dir,CMIP6_file))
CMIP6w = CMIP6[CMIP6.REGION == 'World']

# species
spc_map_SSP = dict(zip(spc_map_csv.dropna().SSP,spc_map_csv.dropna().EPPA)) # note that EPPA names are also the SSP units
CMIP6spc = {}
for s,spc in enumerate(spc_map_SSP.keys()):
    CMIP6spc[spc] = 'CMIP6 Emissions|' + spc
    
# and extract the index of the SSP data columns where the years start 
CMIP6yrs = [int(yr_col) for yr_col in CMIP6.columns if (str.isdigit(yr_col))]
i_yr_CMIP6 = list(CMIP6.columns).index(str(CMIP6yrs[0]))


# populate the sectoral data for Figure 3
em_CMIP6_bysec = {} 
for scen in scen_map_EPPA_CMIP6.values(): 
    CMIP6w_scen = CMIP6w[CMIP6w.SCENARIO==scen] #pick the scenario from the sheet
    em_CMIP6_bysec[scen] = {}
    for spc_cat, spc_weight in spc_map_SSP.items(): # SSP species categories and their molecular weights
        em_CMIP6_bysec[scen][spc_cat] = {}    
        CMIP6w_scen_spc = CMIP6w_scen[CMIP6w_scen.UNIT == ('Mt ' + spc_weight + '/yr')] 
        
        for SSP_sec in sec_map_SSP.keys(): # 9 sectors that overlap with CEDS_GBD-MAPS 
            em_CMIP6_bysec[scen][spc_cat][SSP_sec] = {}    
            # if we have a time series for this scenario/species/sector, include it
            if (CMIP6spc[spc_cat] + '|' + SSP_sec) in str(CMIP6w_scen_spc.VARIABLE):
                em_CMIP6_bysec[scen][spc_cat][SSP_sec] = (CMIP6w_scen_spc[CMIP6w_scen_spc.VARIABLE ==
                                                                          (CMIP6spc[spc_cat] + '|' + SSP_sec)]).iloc[0,i_yr_CMIP6:]
        
# sum the sectors for Figure 4
em_CMIP6_totals = {}
for i_scen, (act_scen,SSP_scen) in enumerate(scen_map_EPPA_CMIP6.items()): 
    em_CMIP6_totals[SSP_scen] = {}
    for spc_SSP, spc_EPPA in spc_map_SSP.items():
        em_CMIP6_totals[SSP_scen][spc_EPPA] = 0.0
        for i_sec,SSP_sec in enumerate(sec_map_SSP.keys()): 
            # sum the sectors, if the trend has data
            if len(em_CMIP6_bysec[SSP_scen][spc_SSP][SSP_sec]) > 1: 
                em_CMIP6_totals[SSP_scen][spc_EPPA] += em_CMIP6_bysec[SSP_scen][spc_SSP][SSP_sec]


# %%


# import the SSP IAM totals for comparison in Figure 2

# first, read in the file and focus on the World regions, radiative forcing variable (to select scenarios)
IAM = pandas.read_csv(os.path.join(input_dir,'SSP_IAM_V2_201811.csv'))
IAMw = IAM[IAM.REGION=='World']
IAMwrf = IAMw[IAMw.VARIABLE=='Diagnostics|MAGICC6|Forcing']

IAMspc = {} 
for s,spc in enumerate(spc_map_SSP.keys()):
    IAMspc[spc] = 'Emissions|' + spc

# scenarios related to EPPA's 1.5 scenario: all the RF2.6 scenarios
IAM15 = IAMw[IAMw.SCENARIO.str.contains('26')]

# scenarios related to EPPA's 2 C scenario: all the RF3.4 scenarios 
IAM2 = IAMw[IAMw.SCENARIO.str.contains('34')]

# scenarios related to EPPA's Paris Forever scenario (within 0.5 of RF5.95 at 2100): could be RF6.0 or baseline scenarios
IAMpf = IAM15[IAM15.SCENARIO == 'notascenario'] 
# start with empty df of the same column format. 
# for each model, find the RF_2100 values in the range, and add that model-scenario combo to the empty df.
for model in IAMw.MODEL.unique():
    m = IAMwrf[IAMwrf.MODEL == model]
    for row in m.itertuples():
        if row[-1] > 5.45 and row[-1] < 6.45:
            IAMpf = IAMpf.append(IAMw[(IAMw.MODEL == row[1]) & (IAMw.SCENARIO == row[2])])

# store the 3 dataframes in the same place to iterate in the figure (using eval(df))
for act_scen,IAM_scen_df in scen_map_EPPA_IAM.items(): 
    # print the number of IAM scenarios in each, for Table 2
    print('number of IAM scenarios compared to the EPPA scenario of',act_scen,'=',
          len(eval(IAM_scen_df)[eval(IAM_scen_df).VARIABLE=='Diagnostics|MAGICC6|Forcing']))

# but first: subtract the SSP sectors that aren't scaled in TAPS
secsubtract = ['Aircraft','Forest Burning','Grassland Burning','Peat Burning']

# dict to map IAM classes to the most representative CMIP6 scenario (in the subset we compare with EPPA, i.e. RF 2.6 to 6.0)
CMIP6SSPdict = {'SSP1':'SSP1-26','SSP2':'SSP2-45','SSP3':'SSP3-LowNTCF',
                'SSP4-26':'SSP4-34',
                'SSP4-34':'SSP4-34',
                'SSP4-60':'SSP4-60',
                'SSP4-Baseline':'SSP4-60',
                'SSP5':'SSP5-34-OS'}

# define SSP year vector as the intersection of CMIP6 and IAM year vectors (2020, 2030...2100)
SSPyrs = [int(yr_col) for yr_col in IAM.columns if (yr_col in CMIP6.columns) and (str.isdigit(yr_col))]
# and extract the index of the SSP data columns where the years start 
i_yr_IAM = list(IAM.columns).index(str(SSPyrs[0]))
i_yr_CMIP6 = list(CMIP6.columns).index(str(SSPyrs[0]))
    
iam_adjusted = {}
for IAM_scen in scen_map_EPPA_IAM.values():
    iam_adjusted[IAM_scen] = {}
    for i_s,spc in enumerate(spc_map_SSP.keys()):
        iam_adjusted[IAM_scen][spc] = []
        # extract all the IAM scenarios that relate to each EPPA scenario 
        paths = eval(IAM_scen)[eval(IAM_scen).VARIABLE==IAMspc[spc]]
        # for every IAM scenario:
        for index, row in paths.iterrows():
            # sum the emissions of non-scaled sectors from the corollary CMIP6 scenario
            iam_ssp = [key for key in CMIP6SSPdict.keys() if key in row.SCENARIO][0]
            sumsecsubtract = 0.0
            for sec in secsubtract:
                if len(CMIP6w[CMIP6w.VARIABLE == CMIP6spc[spc]+'|'+sec])>0: 
                    # extract the scenarios 
                    sumsecsubtract += CMIP6w[(CMIP6w.SCENARIO.str.contains(CMIP6SSPdict[iam_ssp])) &
                                     (CMIP6w.VARIABLE == CMIP6spc[spc]+'|'+sec) ].iloc[:,i_yr_CMIP6:].reset_index().iloc[:,1:].squeeze()

            # find the total from the corollary CMIP6 scenario, and subtract the proportion from non-scaled sectors
            sectotal = CMIP6w[(CMIP6w.SCENARIO.str.contains(CMIP6SSPdict[iam_ssp])) &
                                     (CMIP6w.VARIABLE == CMIP6spc[spc]) ].iloc[:,i_yr_CMIP6:].reset_index().iloc[:,1:].squeeze()

            # and then store the result in ANOTHER dataframe, from which you'll take the max and min and plot them 
            iam_adjusted[IAM_scen][spc].append(row[i_yr_IAM:]*(1-(sumsecsubtract/sectotal)))


# %%


# Figure 2

spc_full = dict(zip(spc_map_csv.dropna().EPPA,spc_map_csv.dropna().Name)) 
spc_formatted = dict(zip(spc_map_csv.dropna().EPPA,spc_map_csv.dropna().Figure)) 
SSPcolor = 'paleturquoise' 
TAPScolor = 'darkviolet'

fig,ax = plt.subplots(len(list_act_scen),len(spc_map_CEDS.keys()),sharex=True,figsize=(24,10)) 

for i_act_scen, (act_scen,IAM_scen) in enumerate(scen_map_EPPA_IAM.items()): # 3 scenarios
    # make axis titles
    ax[i_act_scen,0].set_ylabel('RF~EPPA ' + str(scen_map_EPPAfig[act_scen]),fontsize=20)
    for i_s, (spc_SSP,spc) in enumerate(spc_map_SSP.items()): # 7 species  
        if i_act_scen == 0:
            ax[0,i_s].set_title(spc_formatted.get(spc) + ' (Tg)',fontsize=20) 
            #ax[0,i_s].set_title(spc_full[spc],fontsize=18) # alternate with full names

        # TAPS scenario range 
        em_scen_totalsmin = em_scen_totals[list_act_scen[i_act_scen]][list_int_scen[-1]][spc].astype(None)
        em_scen_totalsmax = em_scen_totals[list_act_scen[i_act_scen]][list_int_scen[0]][spc].astype(None) 
        ax[i_act_scen,i_s].fill_between(yr_list,em_scen_totalsmin,em_scen_totalsmax,alpha=0.2,color=TAPScolor)

        # IAM scenario range
        pathsndf = pandas.DataFrame(iam_adjusted[IAM_scen][spc_SSP])
        ax[i_act_scen,i_s].fill_between(SSPyrs,pathsndf.agg([min]).squeeze(),
                                    pathsndf.agg([max]).squeeze(),alpha=1,color=SSPcolor)

        # TAPS scenario range -- plot again for a better visual 
        ax[i_act_scen,i_s].fill_between(yr_list,em_scen_totalsmin,em_scen_totalsmax,alpha=0.2,color=TAPScolor)

        # after plotting, set y axis limit of zero
        ax[i_act_scen,i_s].set_ylim(ymin=0)

ax[1,3].legend(['TAPS','SSP IAMs'],fontsize=12)

fig_folder = './Figures/'
fig_name = 'f02'
plt.savefig(fig_folder+fig_name+'.png',dpi=300,bbox_inches='tight', format='png')
plt.savefig(fig_folder+fig_name+'.pdf',dpi=300,bbox_inches='tight', format='pdf')
plt.savefig(fig_folder+fig_name+'.svg',dpi=300,bbox_inches='tight', format='svg')
plt.savefig(fig_folder+fig_name+'.eps',dpi=300,bbox_inches='tight', format='eps')


# %%


# Figure 3

# before running: download stackedBarGraph.py into this directory from the following link: https://github.com/minillinim/stackedBarGraph
# as shown here: https://stackoverflow.com/questions/19060144/more-efficient-matplotlib-stacked-bar-chart-how-to-calculate-bottom-values
from matplotlib.ticker import FormatStrFormatter
from stackedBarGraph import StackedBarGrapher
SBG = StackedBarGrapher()

d_labels = ['Base Year','','TAPS CLE','TAPS MFR','SSP Analog']
label_text = ['\\n'.join(d_labels)]
d_colors = ['#d62728', '#2ca02c', '#1f77b4', '#ff7f0e', '#9467bd', '#8c564b', '#e377c2', '#bcbd22', '#7f7f7f']
d_widths = [1.,0.5,1.,1.,1.]

fig,ax = plt.subplots(len(list_act_scen),len(spc_map_CEDS.keys()),sharex=True,figsize=(24,10)) 

for i_act_scen, (act_scen,SSP_scen) in enumerate(scen_map_EPPA_CMIP6.items()): # 3 scenarios    
    # make axis titles
    ax[i_act_scen,0].set_ylabel('RF~EPPA ' + str(scen_map_EPPAfig[act_scen]),fontsize=20)
    for i_s, (spc_SSP,spc_cat) in enumerate(spc_map_SSP.items()): # 7 species
        if i_act_scen == 0:
            ax[0,i_s].set_title(spc_formatted.get(spc_cat) + ' (Tg)',fontsize=20) 
            #ax[0,i_s].set_title(spc_full[spc],fontsize=18) # alternate with full names
            
        # Now, set up a lists of lists for the stacked bars  
        bars = []
        for i_sec,SSP_sec in enumerate(sec_map_SSP.keys()): # 9 sectors that overlap with CEDS_GBD_MAPS
            # for every sector, add that color to each of the "5" bar stacks in d_labels
            bar = []
            bar.append(em_scen_bySSPsec[act_scen][list_int_scen[0]][spc_cat][SSP_sec][yr_list.index(2014)])
            bar.append(0.0)
            bar.append(em_scen_bySSPsec[act_scen][list_int_scen[0]][spc_cat][SSP_sec][yr_list.index(2050)])
            bar.append(em_scen_bySSPsec[act_scen][list_int_scen[-1]][spc_cat][SSP_sec][yr_list.index(2050)])
            if len(em_CMIP6_bysec[SSP_scen][spc_SSP][SSP_sec]) > 1:
                bar.append(em_CMIP6_bysec[SSP_scen][spc_SSP][SSP_sec]['2050'])
            else:
                bar.append(0.0)

            # then, get the bars ready and plot!
            bars.append(bar)
        d = np.array(bars).T 
        SBG.stackedBarPlot(ax[i_act_scen,i_s],
               d,
               d_colors,
               widths=d_widths, gap = 0.1
              )

        # correct y-ticks from having too many decimals
        ax[-1,-1].yaxis.set_major_formatter(FormatStrFormatter('%.f'))

        # and set xlabels on the bottom row
        if i_act_scen == len(list_act_scen)-1:
            ax[-1,i_s].set_xticks(list(np.cumsum(d_widths)-1))
            ax[-1,i_s].set_xticklabels(d_labels,rotation=90,fontsize=20)

# add a legend
import matplotlib.patches as mpatches
legends = []      
for i,(sec) in enumerate(sec_map_SSP.keys()):
    legends.append(mpatches.Patch(color=d_colors[i], label=sec_names_SSP[sec]))
ax[0,0].legend(bbox_to_anchor=(-0.35,1.2), ncol=len(sec_map_SSP.keys()),handles=legends,loc='lower left',fontsize=16.5) 

# save 
fig_folder = './Figures/'
fig_name = 'f03'
plt.savefig(fig_folder+fig_name+'.png',dpi=300,bbox_inches='tight', format='png')
plt.savefig(fig_folder+fig_name+'.pdf',dpi=300,bbox_inches='tight', format='pdf')
plt.savefig(fig_folder+fig_name+'.svg',dpi=300,bbox_inches='tight', format='svg')
plt.savefig(fig_folder+fig_name+'.eps',dpi=300,bbox_inches='tight', format='eps')


# %%


# Figure 4 

regbarnames = ['IND','REA','MES','CHN','ASI','JPN','IDZ','KOR','ANZ','AFR','BRA','LAM','USA','CAN','MEX','EUR','ROE','RUS']
regbarcolors = ['tab:green','lime','darkgreen',
               'deepskyblue','blue','navy','cadetblue','darkturquoise',
               'saddlebrown','yellow',
               'magenta','darkviolet',
               'darkorange','red','darkgoldenrod',
               'silver','gray','black',
                'oldlace'] # last one is for the SSP comparison 
# W/S Asia greens, E Asia blues, ANZ brown, AFR yellow, S. Am purple/pink, N. Am red/orange/brown, Europe grays, SSP "oldlace" 

d_labels = ['Base Year','','TAPS CLE','TAPS MFR','SSP Analog']
label_text = ['\\n'.join(d_labels)]
d_widths = [1.,0.5,1.,1.,1.]
edge_colors = ([None] * len(regbarnames)) + ['black']

fig,ax = plt.subplots(len(list_act_scen),len(spc_map_CEDS.keys()),sharex=True,figsize=(24,10)) 

for i_act_scen, (act_scen,SSP_scen) in enumerate(scen_map_EPPA_CMIP6.items()): # 3 scenarios
    # set axis labels
    ax[i_act_scen,0].set_ylabel('RF~EPPA ' + str(scen_map_EPPAfig[act_scen]),fontsize=20)
    for i_s, (spc_SSP,spc_cat) in enumerate(spc_map_SSP.items()): # 7 species  
        if i_act_scen == 0:
            ax[0,i_s].set_title(spc_formatted.get(spc_cat) + ' (Tg)',fontsize=20) 
            #ax[0,i_s].set_title(spc_full[spc],fontsize=18) # alternate with full names
            
        # Now, set up a lists of lists for the stacked bars  
        bars = []
        for r,reg in enumerate(regbarnames): # EPP7 regions
            # for every region, add that color to each of the "5" bar stacks in d_labels
            bar = []
            bar.append(em_scen_byreg[act_scen][list_int_scen[0]][spc_cat][reg][yr_list.index(2014)])
            bar.append(0)
            bar.append(em_scen_byreg[act_scen][list_int_scen[0]][spc_cat][reg][yr_list.index(2050)])
            bar.append(em_scen_byreg[act_scen][list_int_scen[-1]][spc_cat][reg][yr_list.index(2050)])
            bar.append(0)
            bars.append(bar)

        # and then add the global SSP bars separately 
        bar = []
        bar.append(0)
        bar.append(0)
        bar.append(0)
        bar.append(0)
        bar.append(em_CMIP6_totals[SSP_scen][spc_cat]['2050']) 
        bars.append(bar)

        # then, get the bars ready and plot!
        d = np.array(bars).T
        SBG.stackedBarPlot(ax[i_act_scen,i_s],
                           d,
                           regbarcolors,
                           edgeCols=edge_colors,
                           widths=d_widths, gap = 0.1,
          )

        # correct y-ticks from having too many decimals
        ax[-1,-1].yaxis.set_major_formatter(FormatStrFormatter('%.f'))

        # and set xlabels on the bottom row
        if i_act_scen == len(list_act_scen)-1:
            ax[-1,i_s].set_xticks(list(np.cumsum(d_widths)-1))
            ax[-1,i_s].set_xticklabels(d_labels,rotation=90,fontsize=20)

# add a legend (with EPPA regions and "SSP" marker)
import matplotlib.patches as mpatches
legends = []      
for i,reg_name in enumerate(regbarnames):
    legends.append(mpatches.Patch(color=regbarcolors[i], label=reg_name))
legends.append(mpatches.Patch(color=regbarcolors[-1], label='SSP'))
ax[0,0].legend(bbox_to_anchor=(-0.35,1.2), ncol=len(regbarnames)+1,handles=legends,loc='lower left',fontsize=11) 

# save
fig_folder = './Figures/'
fig_name = 'f04'
plt.savefig(fig_folder+fig_name+'.png',dpi=300,bbox_inches='tight', format='png')
plt.savefig(fig_folder+fig_name+'.pdf',dpi=300,bbox_inches='tight', format='pdf')
plt.savefig(fig_folder+fig_name+'.svg',dpi=300,bbox_inches='tight', format='svg')
plt.savefig(fig_folder+fig_name+'.eps',dpi=300,bbox_inches='tight', format='eps')


# %%
print('figures complete')



