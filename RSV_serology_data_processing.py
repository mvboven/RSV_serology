#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  6 12:59:33 2020

@author: Stijn Andeweg
"""

# %% import packages
import pandas as pd 
import numpy as np
from sas7bdat import SAS7BDAT
import matplotlib.pyplot as plt
import matplotlib as mpl
import datetime
import matplotlib.dates as mdates
from matplotlib.lines import Line2D
from mpl_toolkits.axes_grid1 import make_axes_locatable

# %% set paths
# Set read location serology files
path_raw_data = '/data/BioGrid/projecten/pienter/1_raw_data/' 	#Path raw PIENTER data
PIENTER2 = 'p2voorrsvanalyses.sas7bdat'    			#File name PIENTER2
PIENTER3 = 'p3voorrsvanalyses.sas7bdat'				#File name PIENTER3

# Set figure save path
path_figures = '/data/BioGrid/andewegs/1_figures/python_main_figures/'

# %% import PIENTER2 and PIENTER 3
'''
Open PIENTER2 & 3 SAS files and transform into DataFrame
'''
with SAS7BDAT(path_raw_data + PIENTER2, skip_header=False) as sasfile:
    df_PIENTER2 = sasfile.to_data_frame()

with SAS7BDAT(path_raw_data + PIENTER3, skip_header=False) as sasfile:
    df_PIENTER3 = sasfile.to_data_frame()

'''
Remove two incomplete samples ID = 2120 and 67425
'''

#remove this sample because incomplete MIA IgG, PostF and N absent
df_PIENTER2 = df_PIENTER2[df_PIENTER2.ParticipantId != 2120]   

# patient 67425 IgA incomplete MIA for GB, so I remove all IgA data for this patient
df_PIENTER3.at[df_PIENTER3['id_patient'] == 67425, 'RSV_IgA_GA'] = np.NaN
df_PIENTER3.at[df_PIENTER3['id_patient'] == 67425, 'RSV_IgA_postF'] = np.NaN
df_PIENTER3.at[df_PIENTER3['id_patient'] == 67425, 'RSV_IgA_PreF_NIH'] = np.NaN
df_PIENTER3.at[df_PIENTER3['id_patient'] == 67425, 'RSV_IgA_N'] = np.NaN

#print(list(df_PIENTER2))   #print all variables
#print(list(df_PIENTER3))

'''
Create dictionaries containing for IgG and IgA for the 5 RSV proteins 
'''
IgG = {}      					#make empty dictionary
for index, row in df_PIENTER2.iterrows():		#loop over PIENTER2
    if row['RSV_IgG_GA'] == row['RSV_IgG_GA']:		#check if IgG_GA data is a number and not NaN
        ID = 'P2_' + str(round(row['ParticipantId']))	#use patientID with P2_ in front of it
        IgG[ID] = {}			#Use patientID as key and IgG_GA as value
        IgG[ID]['Ga'] = row['RSV_IgG_GA']
        IgG[ID]['Gb'] = row['RSV_IgG_GB']
        IgG[ID]['PreF'] = row['RSV_IgG_PreF_NIH']
        IgG[ID]['PostF'] = row['RSV_IgG_postF']
        IgG[ID]['N'] = row['RSV_IgG_N']

for index, row in df_PIENTER3.iterrows():		#loop over PIENTER3
    if row['RSV_IgG_GA'] == row['RSV_IgG_GA']:		#check if IgG_GA data is a number and not NaN
        ID = 'P3_' + str(round(row['id_patient']))	#use patientID with P3_ in front of it
        IgG[ID] = {}			#Use patientID as key and IgG_GA as value
        IgG[ID]['Ga'] = row['RSV_IgG_GA']
        IgG[ID]['Gb'] = row['RSV_IgG_GB']
        IgG[ID]['PreF'] = row['RSV_IgG_PreF_NIH']
        IgG[ID]['PostF'] = row['RSV_IgG_postF']
        IgG[ID]['N'] = row['RSV_IgG_N']

IgA = {}      					#make empty dictionary
for index, row in df_PIENTER2.iterrows():		#loop over PIENTER2
    #check if IgG_GA data is a number and not NaN
    if row['RSV_IgA_GA'] == row['RSV_IgA_GA']:
        #use patientID with P2_ in front of it
        ID = 'P2_' + str(round(row['ParticipantId']))	
        IgA[ID] = {}			#Use patientID as key and IgG_GA as value
        IgA[ID]['Ga'] = row['RSV_IgA_GA']
        IgA[ID]['Gb'] = row['RSV_IgA_GB']
        IgA[ID]['PreF'] = row['RSV_IgA_PreF_NIH']
        IgA[ID]['PostF'] = row['RSV_IgA_postF']
        IgA[ID]['N'] = row['RSV_IgA_N']

for index, row in df_PIENTER3.iterrows():		#loop over PIENTER3
    #check if IgG_GA data is a number and not NaN
    if row['RSV_IgA_GA'] == row['RSV_IgA_GA']:
        #use patientID with P3_ in front of it
        ID = 'P3_' + str(round(row['id_patient']))	
        IgA[ID] = {}			#Use patientID as key and IgG_GA as value
        IgA[ID]['Ga'] = row['RSV_IgA_GA']
        IgA[ID]['Gb'] = row['RSV_IgA_GB']
        IgA[ID]['PreF'] = row['RSV_IgA_PreF_NIH']
        IgA[ID]['PostF'] = row['RSV_IgA_postF']
        IgA[ID]['N'] = row['RSV_IgA_N']

# Print strange output and change to lowest possible value. 
# Not a problem because I do not use Ga and Gb.
for k, v in IgA.items():
    for k2, v2 in v.items():
        if v2 == 0:
            print(k,k2,v2)
            IgA[k][k2] = 0.01   #set 0 to lowest possible value
            
# %% make dictionary with time related data
'''
 Make a dictionary with all time/age/date related variables.
'''
onnodig_nan = []
#dictionary age in months: lftinmnd[P2/3_IDpatient] = age_in_months
time_dict = 	{}
for index, row in df_PIENTER2.iterrows():		     #loop over PIENTER2
    ID = 'P2_' + str(round(row['ParticipantId'])) 	 #set ID with P2_ infront
    time_dict[ID] = {}
    birthday = row['gebdat2']				         #extract birthday
    cons_day = row['consultdatum']			         #extract consult day
    delta = cons_day - birthday				         #calculate diference between days
    time_dict[ID]['age_days'] = delta.days		 	 #set age in days at consult day
    age_months_c = delta.days/30.44                  #average month has 30.44 days
    time_dict[ID]['age_months']  = age_months_c      #continues months 
    time_dict[ID]['age_months_d'] = row['lftinmnd2'] #assign age in months to patientID
    time_dict[ID]['age_years'] = row['LFTB']  		 #assign age in years to patientID
    consult_day = row['consultdatum']			     #extract consult day
    time_dict[ID]['consultdate'] = consult_day	     #set age in days at consult day
    birthday = row['gebdat2']			             #extract consult day
    time_dict[ID]['birthday'] = birthday			     #set age in days at consult day

for index, row in df_PIENTER3.iterrows():		     #loop over PIENTER3
    ID = 'P3_' + str(round(row['id_patient']))		 #set ID with P3_ infront
    time_dict[ID] = {}
    birthday = row['dateofbirth_true']			     #extract birthday
    cons_day = row['consultdate']			         #extract consultday
    
    if cons_day is not None:				             #only calculate difference if there is a consultdate
        delta = cons_day - birthday			         #calculate diffference between days
        time_dict[ID]['age_days'] = delta.days		 #set age in days at consult day
        age_months_c = delta.days/30.44              #average month has 30.44 days
        time_dict[ID]['age_months']  = age_months_c  #continues months 
        
    time_dict[ID]['age_months_d']  = row['age_months']#assign age in months to patientID, months naar beneden afgerond.
    time_dict[ID]['age_years'] = row['age_year']		 #assign age in years to patientID        
    consult_day = row['consultdate']			         #extract consultday
    time_dict[ID]['consultdate'] = consult_day		 #set age in days at consult day
    birthday = row['dateofbirth_true']			     #extract consultday
    time_dict[ID]['birthday'] = birthday			     #set age in days at consult day

# %% make dictionary with additional variables

'''
'pershuis5', 'huisgenoot0tot17', 'deelnemer0tot4', 'huisgenoot0tot4', 
'huisgenoot1', 'huisgenoot2', 'huisgenoot3', 'n_huisgenoot0tot4', 'kind0tot4', 
'crechhuis_pos', 'crechkind_pos (1=ja, 0=nee)', 'crech', 'crechdeel', 'tot_kind0tot4', 'n_dagdeel_per_kind'
'''
# list of variables related to the age of household
householdage = ['householdage1','householdage2','householdage3','householdage4'\
        ,'householdage5','householdage6','householdage7','householdage8',
        'householdage9','householdage10']
# create a dictionary with additional variables
additional_variables = {}
for index, row in df_PIENTER3.iterrows():		#loop over PIENTER3
    ID = 'P3_' + str(round(row['id_patient']))	#use patientID with P3_ in front of it
    additional_variables[ID] = {}
    
    geslacht = row['geslacht']
    additional_variables[ID]['sex']  = geslacht
    if row['pregnancytime'] == row['pregnancytime']: 
        if row['pregnancytime'] > 3:            #remove the strange pregnancytime of 3 weeks for 1 sample. 
            additional_variables[ID]['pregnancytime'] = row['pregnancytime'] # in weeks
    
    if row['breastfeeding'] == row['breastfeeding']: #1=Yes|2=No|3=I don't know
        additional_variables[ID]['breastfeeding'] = row['breastfeeding']
    
    if row['howlongbreastfed'] == row['howlongbreastfed']: # months
        additional_variables[ID]['howlongbreastfed'] = row['howlongbreastfed']  
    
    if row['household'] == row['household']:
        additional_variables[ID]['n_household'] = int(row['household'])
    
    if row['visitnursery_child'] == row['visitnursery_child']: #1=yes, 2=no
        if row['visitnursery_child'] == 2.0:
            additional_variables[ID]['visitnursery_child'] = 0.0 # 0=no, 1=yes
        else:
            additional_variables[ID]['visitnursery_child'] = row['visitnursery_child']
    
    if row['visitnursery'] == row['visitnursery']: 
        if row['visitnursery'] == 2.0:
            additional_variables[ID]['visitnursery_house'] = 0.0 # 0=no, 1=yes
        else:
            additional_variables[ID]['visitnursery_house'] = row['visitnursery']
    
    if row['otherchildren'] == row['otherchildren']:
        if row['otherchildren'] == 2.0:
            additional_variables[ID]['visitnursery_house'] = 0.0 # 0=no, 1=yes
        else:
            additional_variables[ID]['visitnursery_house'] = row['otherchildren']
    
    ages_housemates = []
    for item in householdage:
        if row[item] == row[item]:
            ages_housemates.append(row[item])
            additional_variables[ID]['ages_housemates'] = ages_housemates
    
    if row['household'] - len(ages_housemates) == 1 and row['age_year'] not in ages_housemates:
        ages_housemates.append(row['age_year'])
        additional_variables[ID]['ages_housemates'] = ages_housemates
    household02 = 0
    household03 = 0
    household04 = 0
    household59 = 0
    for i in ages_housemates:
        if i in range(0,3):
            household02 += 1
        if i in range(0,4):
            household03 += 1
        if i in range(0,5):
            household04 += 1
        if i in range(5,10):
            household59 += 1
    additional_variables[ID]['household02'] = household02
    additional_variables[ID]['household03'] = household03
    additional_variables[ID]['household04'] = household04
    additional_variables[ID]['household59'] = household59
    if row['contacttotal'] == row['contacttotal']:
        additional_variables[ID]['contacttotal'] = row['contacttotal']
        additional_variables[ID]['contactday'] = row['contactday'] # 1=Monday|2=Tuesday|3=Wednesday|4=Thursday|5=Friday|6=Saturday|7=Sunday
    if row['contactday'] == row['contactday']:
        contact04 = 0
        for i in ['contact04ymen1', 'contact04ywom1']:
            if row[i] == row[i]:
                contact04 += row[i]
        additional_variables[ID]['contact04'] = contact04
        contact59 = 0
        for i in ['contact59ymen1', 'contact59ywom1']:
            if row[i] == row[i]:
                contact59 += row[i]
        additional_variables[ID]['contact59'] = contact59
    
    if not bool(additional_variables[ID]):
        del additional_variables[ID]
# list of variables related to the age of household
householdage = ['Lftpershuis1', 'Lftpershuis2', 'Lftpershuis3', 'Lftpershuis4', 'Lftpershuis5',\
        'Lftpershuis6', 'Lftpershuis7', 'Lftpershuis8', 'Lftpershuis9','Lftpershuis10']
for index, row in df_PIENTER2.iterrows():		#loop over PIENTER3
    ID = 'P2_' + str(round(row['ParticipantId']))	#use patientID with P3_ in front of it
    additional_variables[ID] = {}
    
    geslacht = row['gesl2']
    additional_variables[ID]['sex']  = geslacht
    if row['pershuis'] == row['pershuis']:
        if row['pershuis'].isnumeric():
            additional_variables[ID]['n_household'] = int(row['pershuis'])
    if row['crechkind_pos'] == row['crechkind_pos']:
        additional_variables[ID]['visitnursery_child'] = row['crechkind_pos']
    if row['crechhuis_pos'] == row['crechhuis_pos']:
        additional_variables[ID]['visitnursery_house'] = row['crechhuis_pos']
    ages_housemates = []
    for item in householdage:
        if row[item].isdigit():
            ages_housemates.append(int(row[item]))
            additional_variables[ID]['ages_housemates'] = ages_housemates

    if row['pershuis'] == row['pershuis'] and row['pershuis'].isnumeric():
        if int(row['pershuis'])-len(ages_housemates) == 1 and row['LFTB'] not in ages_housemates:
            ages_housemates.append(round(row['LFTB']))
            additional_variables[ID]['ages_housemates'] = ages_housemates
    household02 = 0
    household03 = 0
    household04 = 0
    household59 = 0
    for i in ages_housemates:
        if i in range(0,3):
            household02 += 1
        if i in range(0,4):
            household03 += 1
        if i in range(0,5):
            household04 += 1
        if i in range(5,10):
            household59 += 1
    additional_variables[ID]['household02'] = household02
    additional_variables[ID]['household03'] = household03
    additional_variables[ID]['household04'] = household04
    additional_variables[ID]['household59'] = household59
    
    if row['contdag'] == row['contdag']:
        additional_variables[ID]['contactday'] = row['contdag']
        if row['contpers'].isdigit():
            additional_variables[ID]['contacttotal'] = row['contpers']  #'contpers_nieuw_v2', 'contpers_v2' zijn deze beter?
        else:
            additional_variables[ID]['contacttotal'] = 0
        if row['contpers0tot4'].isdigit():
            additional_variables[ID]['contact04'] = row['contpers0tot4']
        else:
            additional_variables[ID]['contact04'] = 0
        if row['contpers5tot9'].isdigit():
            additional_variables[ID]['contact59'] = row['contpers5tot9']
        else:
            additional_variables[ID]['contact59'] = 0
    if not bool(additional_variables[ID]):
        del additional_variables[ID]


# %% Figure histogram of age in months distribution all data

max_months = 61

#histogram distribution of age
list_age_year_IgG = []			#create empty list for IgG
for k in IgG.keys():			    #loop over IgG_GA keys
    if time_dict[k]['age_months'] < max_months:
        #append age in months to list IgG
        list_age_year_IgG.append(time_dict[k]['age_months'])
list_age_year_IgA = []			#create empty list for IgA
for k in IgA.keys():			    #loop over IgA_GA keys
    if time_dict[k]['age_months'] < max_months:
        #append age in months to list IgA
        list_age_year_IgA.append(time_dict[k]['age_months'])

# Set bins
bins = [-0.000000001,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,\
        21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,
        41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61
        ]

# Plot Supplement Figure 1
plt.figure(dpi=300, figsize=(8,4))

plt.xticks(np.arange(min(bins), max(bins), 5))
plt.yticks(np.arange(0, 91, 10))
plt.hist(list_age_year_IgG, bins, facecolor='g', label='IgG',\
         rwidth = 0.9)	        # plot histogram IgG
plt.hist(list_age_year_IgA, bins, facecolor='orange', label='IgG & IgA',\
         rwidth = 0.9)	        # plot histogram IgG
plt.legend()				        # make legend

plt.xlabel('Age bins (months)')	# make x-axis label
plt.ylabel('Counts')				# make y-axis label

# vertical line at age = 500 days
plt.vlines((500/30.44),0,100,linestyles='dotted',colors='r')
plt.ylim(0,95)                  # set y limits

plt.savefig(path_figures + 'histogram_age_months_distr.svg', \
            format='svg', bbox_inches='tight')
plt.show()

# %% Produce dataframe IgA and IgG samples
max_months = 61   #cut-off age in months

# create empty dataframe
index = range(len([(k) for k in IgA.keys() if time_dict[k]['age_months'] < max_months]))
df_IgAandG = pd.DataFrame(columns=['ID','age_months', 'age_days', 'IgG_PreF', 'IgA_PreF',\
            'IgG_PostF', 'IgA_PostF', 'IgG_Ga', 'IgA_Ga', 'IgG_Gb', 'IgA_Gb',\
                'IgG_N', 'IgA_N', 'ratio_PreF', 'ratio_PostF', 'ratio_N',\
                    'ratio_Ga', 'ratio_Gb'], index=index)
# fill the dataframe		
index_count = 0
for k, v in IgA.items():
    if k in IgG.keys():
        if time_dict[k]['age_months'] < max_months:
        #print(k)
            df_IgAandG.loc[index_count].ID = k
            df_IgAandG.loc[index_count].age_months = time_dict[k]['age_months']
            df_IgAandG.loc[index_count].age_days = time_dict[k]['age_days']
            df_IgAandG.loc[index_count].IgG_PreF = IgG[k]['PreF']
            df_IgAandG.loc[index_count].IgA_PreF = IgA[k]['PreF']
            df_IgAandG.loc[index_count].IgG_PostF = IgG[k]['PostF']
            df_IgAandG.loc[index_count].IgA_PostF = IgA[k]['PostF']
            df_IgAandG.loc[index_count].IgG_Ga = IgG[k]['Ga']
            df_IgAandG.loc[index_count].IgA_Ga = IgA[k]['Ga']
            df_IgAandG.loc[index_count].IgG_Gb = IgG[k]['Gb']
            df_IgAandG.loc[index_count].IgA_Gb = IgA[k]['Gb']
            df_IgAandG.loc[index_count].IgG_N = IgG[k]['N']
            df_IgAandG.loc[index_count].IgA_N = IgA[k]['N']
            df_IgAandG.loc[index_count].ratio_PreF = IgA[k]['PreF'] / IgG[k]['PreF']
            df_IgAandG.loc[index_count].ratio_PostF = IgA[k]['PostF'] / IgG[k]['PostF']
            df_IgAandG.loc[index_count].ratio_N = IgA[k]['N'] / IgG[k]['N']
            df_IgAandG.loc[index_count].ratio_Ga = IgA[k]['Ga'] / IgG[k]['Ga']
            df_IgAandG.loc[index_count].ratio_Gb = IgA[k]['Gb'] / IgG[k]['Gb']
            index_count += 1
  
index = range(len([(k) for k in IgG.keys() if time_dict[k]['age_months'] < max_months]))
df_IgG = pd.DataFrame(columns=['ID','age_months', 'age_days', 'IgG_PreF', \
            'IgG_PostF', 'IgG_Ga','IgG_Gb', 'IgG_N'], index=index)

index_count = 0
for k, v in IgG.items():
    if time_dict[k]['age_months'] < max_months:
        #print(k)
        df_IgG.loc[index_count].ID = k
        df_IgG.loc[index_count].age_months = time_dict[k]['age_months']
        df_IgG.loc[index_count].age_days = time_dict[k]['age_days']
        df_IgG.loc[index_count].IgG_PreF = IgG[k]['PreF']
        df_IgG.loc[index_count].IgG_PostF = IgG[k]['PostF']
        df_IgG.loc[index_count].IgG_Ga = IgG[k]['Ga']
        df_IgG.loc[index_count].IgG_Gb = IgG[k]['Gb']
        df_IgG.loc[index_count].IgG_N = IgG[k]['N']
        index_count += 1
        
# %% Figure scatter plot antibody concentration IgG (y) vs age (days) (x), colored by IgA
# Figure did not end up in the manuscript
plt.figure(dpi=300, figsize=(20,10))

RSV_proteins = ['PreF', 'PostF', 'N', 'Ga', 'Gb']
position_list = [1,4,2,3,6]
for inx, i in enumerate(RSV_proteins):
    ax = plt.subplot(2, 3, position_list[inx])
    plt.title(i)
    plt.scatter(np.array(df_IgAandG['age_days']), np.array(df_IgAandG['IgG_' + i]), \
                c=np.array(df_IgAandG['IgA_' + i]), marker="1", norm=mpl.colors.LogNorm())
    plt.yscale('log')   
    plt.ylabel('AU/ml (IgG)')
    cbar = plt.colorbar()
    cbar.set_label('IgA_'+ i)
    if inx == 0 or inx == 3:
        plt.setp(ax.get_xticklabels(), visible=False)
    else:
        plt.xlabel('Age (days)')
plt.tight_layout()
plt.savefig(path_figures + 'scatter_RSV_proteins_IgG_days.svg')
plt.show()   


# %% Figure scatter plot antibody concentration IgG (y) vs age (days) (x), colored by IgA
# Figure did not end up in the manuscript
#fig2 = plt.figure(dpi=300, figsize=(16,8))

fig, axes = plt.subplots(nrows=1, ncols=3, sharey=True, figsize=(18,6))

RSV_proteins = ['PreF', 'PostF', 'N']
position_list = [0,1,2]
for inx, i in enumerate(RSV_proteins):
   #ax = plt.subplot(2, 3, position_list[inx])
    axes[inx].set_title(i)
    im = axes[inx].scatter(np.array(df_IgAandG['age_days']), np.array(df_IgAandG['IgG_' + i]), \
                c=np.array(df_IgAandG['IgA_' + i]), marker="1", norm=mpl.colors.LogNorm())
    axes[inx].set_yscale('log')   
    axes[inx].set_xlabel('Age (days)')
    axes[inx].set_ylim(10**-2, 10**4)
axes[0].set_ylabel('AU/ml (IgG)')

fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.35, -0.05, 0.35, 0.05])

cbar = plt.colorbar(im, orientation="horizontal", cax=cbar_ax)
cbar.ax.set_xlabel('IgA AU/ml')

plt.tight_layout()
plt.savefig(path_figures + 'scatter_pre_post_N_IgG_days.svg', bbox_inches='tight')
plt.show()   
# %% Figure scatter plot antibody concentration IgG (y) vs age (days) (x), colored by IgA for PostF

plt.figure(dpi=300, figsize=(10,6))
plt.title('PostF')
plt.scatter(np.array(df_IgAandG['age_days']), np.array(df_IgAandG['IgG_PostF']), \
                c=np.array(df_IgAandG['IgA_' + 'PostF']), marker="1", norm=mpl.colors.LogNorm())
plt.yscale('log')   
plt.ylabel('AU/ml (IgG)')
cbar = plt.colorbar()
cbar.set_label('IgA_'+ 'PostF')
plt.xlabel('Age (days)')
plt.savefig(path_figures + 'scatter_PostF_IgG_days_all_ages.svg')
plt.show()   

# %% 2x3 Scatter plot only the proteins we use for the infected/non-infected classification

titles_RSV = ['PreF','PostF','N']
plt.figure(dpi=300, figsize=(16,8))
for index, i in enumerate(titles_RSV):
    ax1 = plt.subplot(2,3,index+1)
    ax1.set_yscale('log')
    
    #Histogram on the side
    divider = make_axes_locatable(ax1)
    ax1Histy = divider.append_axes("right", size=0.6, pad=0.1, sharey=ax1)
    
    logbins = np.logspace(-2, 4, num=31)
    ax1Histy.hist(df_IgG['IgG_' + i].loc[df_IgG['IgG_' + i] <= 1.0]\
                 .loc[df_IgG['age_days'] > 500]\
                 , bins = logbins, orientation='horizontal', color = 'C0')
    ax1Histy.hist(df_IgG['IgG_' + i].loc[df_IgG['IgG_' + i] > 1.0]\
                 .loc[df_IgG['age_days'] > 500]\
                 , bins = logbins, orientation='horizontal', color = 'C1')
    ax1Histy.tick_params(axis="y", labelleft=False)

    ax1Histy.set_xticks(np.arange(0,51, step = 25))
    for ind in range(len(df_IgG['IgG_' + i])):
        if df_IgG['age_days'][ind] <= 500:
            ax1.scatter(df_IgG['age_days'][ind], df_IgG['IgG_' + i][ind],\
                        color = 'C7' ,marker="1")
        elif df_IgG['IgG_' + i][ind] > 1.0:
            ax1.scatter(df_IgG['age_days'][ind], df_IgG['IgG_' + i][ind],\
                        color = 'C1',marker="1")
        else:
            ax1.scatter(df_IgG['age_days'][ind], df_IgG['IgG_' + i][ind],\
                        color = 'C0',marker="1")
    
    ax1.set_title(i)
    if index < 3:
        ax1.vlines(500,10**-2,10**4,linestyles='dotted',colors='r')
        ax1.hlines(10**0,500,max(df_IgG['age_days']),linestyles='dashed',colors='black') #IgG infected cut off
        ax1.text(1500,10**0.1,'1.0')
        ax1.set_ylim(10**-2,10**4)
    
    ax2 = plt.subplot(2,3,index+4)
    #make y hist
    divider = make_axes_locatable(ax2)
    ax2Histy = divider.append_axes("right", size=0.6, pad=0.1, sharey=ax2)
    
    logbins2 = np.logspace(-3, 3, num=31)
    ax2Histy.hist(df_IgAandG['IgA_' + i].loc[df_IgAandG['IgA_' + i] <= 0.2]\
                 .loc[df_IgAandG['age_days'] <= 500]\
                 , bins = logbins2, orientation='horizontal', color = 'C0')
    ax2Histy.hist(df_IgAandG['IgA_' + i].loc[df_IgAandG['IgA_' + i] > 0.2]\
                 .loc[df_IgAandG['age_days'] <= 500]\
                 , bins = logbins2, orientation='horizontal', color = 'C1')
    ax2Histy.tick_params(axis="y", labelleft=False)
    ax2Histy.set_xlabel('Counts')
    
    
    for ind in range(len(df_IgAandG['IgA_' + i])):
        if df_IgAandG['age_days'][ind] > 500:
            ax2.scatter(df_IgAandG['age_days'][ind], df_IgAandG['IgA_' + i][ind],\
                        marker="1", c='C7')
        elif df_IgAandG['IgA_' + i][ind] > 0.2:
            ax2.scatter(df_IgAandG['age_days'][ind], df_IgAandG['IgA_' + i][ind],\
                        color = 'C1',marker="1")
        else:
            ax2.scatter(df_IgAandG['age_days'][ind], df_IgAandG['IgA_' + i][ind],\
                        color = 'C0',marker="1")
    
    if index < 3:
        ax2.hlines(0.2,min(df_IgAandG['age_days']),500,linestyles='dashed',colors='black') #IgG infected cut off
        ax2.text(400,.25,'0.2')
        ax2.vlines(500,10**-3,10**3,linestyles='dotted',colors='r')
        ax2.set_ylim(10**-3,10**3)
    else:
        ax2.set_ylim(10**-3,10**2)
    
    ax1.grid(True, which='major')
    ax2.grid(True, which='major')
    
    ax1Histy.grid(True, which='major')
    ax2Histy.grid(True, which='major')
    
    ax2.set_yscale('log')
    if index in [1,2]:
        ax1.yaxis.set_ticklabels([])
        ax2.yaxis.set_ticklabels([])
        ax2Histy.set_xticks(np.arange(0,151, step = 75))
        
    ax1.xaxis.set_ticklabels([])
    ax2.set_xlabel('Age (days)')
    
    if index == 0:
        ax1.set_ylabel('AU/ml (IgG)')
        ax2.set_ylabel('AU/ml (IgA)')
        ax2Histy.set_xticks(np.arange(0,251, step = 125))
    
    

plt.tight_layout(pad=0.1)    
plt.savefig(path_figures + 'scatter_RSV_for_classification_withcolors_days2x3.svg')
plt.show()

# %% produce dataframe with infection classification
max_months = 61	   #we are intrested in indivduals below 5 years
cut_off_days_IgG = 500

index = range(len([(k) for k in IgG.keys() if time_dict[k]['age_months'] \
                   < max_months\
                  and time_dict[k]['age_days'] > cut_off_days_IgG]) \
                  + len([(k) for k in IgA.keys() if time_dict[k]['age_days'] \
                         <= cut_off_days_IgG]))
df_infection = pd.DataFrame(columns=['ID','age_months' ,'age_days','birthday',\
            'consultdate', 'infection', 'pregnancytime', 'IgG_PreF', 'IgA_PreF',\
            'IgG_PostF', 'IgA_PostF', 'IgG_Ga', 'IgA_Ga', 'IgG_Gb', 'IgA_Gb',\
                'IgG_N', 'IgA_N'], index=index)


#new test
index_count = 0
#define infected
RSV_selection = ['PreF','PostF','N']
for k, v in IgG.items():
    if all( [time_dict[k]['age_months'] < max_months, time_dict[k]['age_days'] \
             > cut_off_days_IgG]):
        df_infection.loc[index_count].ID = k
        df_infection.loc[index_count].age_months = time_dict[k]['age_months']
        df_infection.loc[index_count].age_days = time_dict[k]['age_days']
        df_infection.loc[index_count].consultdate = time_dict[k]['consultdate']
        df_infection.loc[index_count].birthday = time_dict[k]['birthday']
        df_infection.loc[index_count].IgG_PreF = IgG[k]['PreF']
        df_infection.loc[index_count].IgG_PostF = IgG[k]['PostF']
        df_infection.loc[index_count].IgG_Ga = IgG[k]['Ga']
        df_infection.loc[index_count].IgG_Gb = IgG[k]['Gb']
        df_infection.loc[index_count].IgG_N = IgG[k]['N']
        if k in additional_variables.keys():
            if 'pregnancytime' in additional_variables[k].keys():
                df_infection.loc[index_count].pregnancytime = \
                    additional_variables[k]['pregnancytime']
        
        true_count = sum(1 for item in RSV_selection if v[item] > 1.0)
        if true_count > 1: # 2 or 3 above cutoff
            df_infection.loc[index_count].infection = 1
        else:
            df_infection.loc[index_count].infection = 0
        if k in IgA.keys():
            df_infection.loc[index_count].IgA_PreF = IgA[k]['PreF']
            df_infection.loc[index_count].IgA_PostF = IgA[k]['PostF']
            df_infection.loc[index_count].IgA_Ga = IgA[k]['Ga']
            df_infection.loc[index_count].IgA_Gb = IgA[k]['Gb']
            df_infection.loc[index_count].IgA_N = IgA[k]['N']
        else:
            df_infection.loc[index_count].IgA_PreF = np.NaN
            df_infection.loc[index_count].IgA_PostF = np.NaN
            df_infection.loc[index_count].IgA_Ga = np.NaN
            df_infection.loc[index_count].IgA_Gb = np.NaN
            df_infection.loc[index_count].IgA_N = np.NaN
        index_count += 1

for k, v in IgA.items():
    if time_dict[k]['age_days'] <= cut_off_days_IgG:
        df_infection.loc[index_count].ID = k
        df_infection.loc[index_count].age_months = time_dict[k]['age_months']
        df_infection.loc[index_count].age_days = time_dict[k]['age_days']
        df_infection.loc[index_count].consultdate = time_dict[k]['consultdate']
        df_infection.loc[index_count].birthday = time_dict[k]['birthday']
        df_infection.loc[index_count].IgA_PreF = IgA[k]['PreF']
        df_infection.loc[index_count].IgA_PostF = IgA[k]['PostF']
        df_infection.loc[index_count].IgA_Ga = IgA[k]['Ga']
        df_infection.loc[index_count].IgA_Gb = IgA[k]['Gb']
        df_infection.loc[index_count].IgA_N = IgA[k]['N']
        df_infection.loc[index_count].IgG_PreF = IgG[k]['PreF']
        df_infection.loc[index_count].IgG_PostF = IgG[k]['PostF']
        df_infection.loc[index_count].IgG_Ga = IgG[k]['Ga']
        df_infection.loc[index_count].IgG_Gb = IgG[k]['Gb']
        df_infection.loc[index_count].IgG_N = IgG[k]['N']
        if k in additional_variables.keys():
            if 'pregnancytime' in additional_variables[k].keys():
                df_infection.loc[index_count].pregnancytime = \
                    additional_variables[k]['pregnancytime']
                
        true_count = sum(1 for item in RSV_selection if IgA[k][item] > 0.2)#0.06683439)
        if true_count > 1: # 2 or 3 above cutoff
            df_infection.loc[index_count].infection = 1
        else:
            df_infection.loc[index_count].infection = 0
        index_count += 1

# %% 2x3 Scatter plot only the proteins we use for the infected/non-infected classification


IgG_colouring = []
for index, row in df_IgG.iterrows():
    #if row['age_days'] > 500:
    if row['ID'] in list(df_infection['ID']):
        if int(df_infection['infection'].loc[df_infection['ID'] == row['ID']]) == 1:
            IgG_colouring.append('C1')
        if int(df_infection['infection'].loc[df_infection['ID'] == row['ID']]) == 0:
            IgG_colouring.append('C0')   
    else:
        IgG_colouring.append('darkgrey')
df_IgG['IgG_colouring'] = IgG_colouring 
  
IgA_colouring = []
for index, row in df_IgAandG.iterrows():
    #print(row['ID'])
    #if row['age_days'] <= 500:
        
    if int(df_infection['infection'].loc[df_infection['ID'] == row['ID']]) == 1:
        IgA_colouring.append('C1')
    if int(df_infection['infection'].loc[df_infection['ID'] == row['ID']]) == 0:
        IgA_colouring.append('C0')   
    #else:
     #   IgA_colouring.append('C7')
df_IgAandG['IgA_colouring'] = IgA_colouring 


    
titles_RSV = ['PreF','PostF','N']
plt.figure(dpi=300, figsize=(16,8))
for index, i in enumerate(titles_RSV):
    ax1 = plt.subplot(2,3,index+1)
    ax1.set_yscale('log')
    
    #Histogram on the side
    divider = make_axes_locatable(ax1)
    ax1Histy = divider.append_axes("right", size=0.6, pad=0.1, sharey=ax1)
    
    logbins = np.logspace(-2, 4, num=31)
    ax1Histy.hist(df_IgG['IgG_' + i].loc[df_IgG['IgG_colouring'] == 'C0']\
                 .loc[df_IgG['age_days'] > 500], stacked=True\
                 , bins = logbins, orientation='horizontal', color = 'C0')
    ax1Histy.hist(df_IgG['IgG_' + i].loc[df_IgG['IgG_colouring'] == 'C1']\
                 .loc[df_IgG['age_days'] > 500], stacked=True\
                 , bins = logbins, orientation='horizontal', color = 'C1')
    ax1Histy.tick_params(axis="y", labelleft=False)

    ax1Histy.set_xticks(np.arange(0,51, step = 25))
    
    
    ax1.scatter(df_IgG['age_days'], df_IgG['IgG_' + i],\
                        color = df_IgG['IgG_colouring'] ,marker="1", alpha = 0.8)

    
    ax1.set_title(i)
    if index < 3:
        ax1.vlines(500,10**-2,10**4,linestyles='dotted',colors='r')
        ax1.hlines(10**0,500,max(df_IgG['age_days']),linestyles='dashed',colors='black') #IgG infected cut off
        ax1.text(1500,10**0.1,'1.0')
        ax1.set_ylim(10**-2,10**4)
    
    ax2 = plt.subplot(2,3,index+4)
    #make y hist
    divider = make_axes_locatable(ax2)
    ax2Histy = divider.append_axes("right", size=0.6, pad=0.1, sharey=ax2)
    
    logbins2 = np.logspace(-3, 3, num=31)
    ax2Histy.hist(df_IgAandG['IgA_' + i].loc[df_IgAandG['IgA_colouring'] == 'C0']\
                 .loc[df_IgAandG['age_days'] <= 500], stacked=True\
                 , bins = logbins2, orientation='horizontal', color = 'C0')
    ax2Histy.hist(df_IgAandG['IgA_' + i].loc[df_IgAandG['IgA_colouring'] == 'C1']\
                 .loc[df_IgAandG['age_days'] <= 500], stacked=True\
                 , bins = logbins2, orientation='horizontal', color = 'C1')
    ax2Histy.tick_params(axis="y", labelleft=False)
    ax2Histy.set_xlabel('Counts')
    
    # the scatter plot
    ax2.scatter(df_IgAandG['age_days'], df_IgAandG['IgA_' + i],\
                        color = df_IgAandG['IgA_colouring'] ,marker="1", alpha = 0.8)
    
    if index < 3:
        ax2.hlines(0.2,min(df_IgAandG['age_days']),500,linestyles='dashed',colors='black') #IgG infected cut off
        ax2.text(400,.25,'0.2')
        ax2.vlines(500,10**-3,10**3,linestyles='dotted',colors='r')
        ax2.set_ylim(10**-3,10**3)
    else:
        ax2.set_ylim(10**-3,10**2)
    
    ax1.grid(True, which='major')
    ax2.grid(True, which='major')
    
    ax1Histy.grid(True, which='major')
    ax2Histy.grid(True, which='major')
    
    ax2.set_yscale('log')
    if index in [1,2]:
        ax1.yaxis.set_ticklabels([])
        ax2.yaxis.set_ticklabels([])
        ax2Histy.set_xticks(np.arange(0,151, step = 75))
        
    ax1.xaxis.set_ticklabels([])
    ax2.set_xlabel('Age (days)')
    
    if index == 0:
        ax1.set_ylabel('AU/ml (IgG)')
        ax2.set_ylabel('AU/ml (IgA)')
        ax2Histy.set_xticks(np.arange(0,251, step = 125))
    
    

plt.tight_layout(pad=0.1)    
plt.savefig(path_figures + 'scatter_RSV_for_classification_withcolors_days2x3.svg')
plt.show()


# %% added additional variables to infection dataframe (maybe want to undo this)
df_infection['age_group'] = pd.Series(['none' for x in \
                                       range(len(df_infection['age_months']))]\
                                      ,index=df_infection.index)
for item in ['n_household', 'visitnursery_child', 'visitnursery_house',\
             'contacttotal', 'contactday', 'contact04', 'contact59','household04',\
             'household03','household02','household59', 'sex']:
    df_infection[item] = pd.Series([np.nan for x in \
                                    range(len(df_infection['age_months']))]\
                                  ,index=df_infection.index)

    
for index, row in df_infection.iterrows():
    
    if row['ID'] in additional_variables.keys():
        df_infection.at[index, 'sex'] = additional_variables[row['ID']]['sex']
        if 'n_household' in additional_variables[row['ID']].keys():
            if additional_variables[row['ID']]['n_household'] < 20: #remove the strange big households
                df_infection.at[index, 'n_household'] = \
                    additional_variables[row['ID']]['n_household']
        if 'visitnursery_child' in additional_variables[row['ID']].keys():
            df_infection.at[index, 'visitnursery_child'] = \
                additional_variables[row['ID']]['visitnursery_child']
        if 'visitnursery_house' in additional_variables[row['ID']].keys():
            df_infection.at[index, 'visitnursery_house'] = \
                additional_variables[row['ID']]['visitnursery_house']
        if 'contacttotal' in additional_variables[row['ID']].keys():
            df_infection.at[index, 'contacttotal'] = \
                additional_variables[row['ID']]['contacttotal']
        if 'contactday' in additional_variables[row['ID']].keys():
            df_infection.at[index, 'contactday'] = \
                additional_variables[row['ID']]['contactday']
        if 'contact04' in additional_variables[row['ID']].keys():
            df_infection.at[index, 'contact04'] = \
                additional_variables[row['ID']]['contact04']
        if 'contact59' in additional_variables[row['ID']].keys():
            df_infection.at[index, 'contact59'] = \
                additional_variables[row['ID']]['contact59']
        if 'household04' in additional_variables[row['ID']].keys():
            if additional_variables[row['ID']]['household04'] >= 1:
                df_infection.at[index, 'household04'] = \
                    additional_variables[row['ID']]['household04']
                df_infection.at[index, 'household03'] = \
                    additional_variables[row['ID']]['household03']
                df_infection.at[index, 'household02'] = \
                    additional_variables[row['ID']]['household02']
        if 'household59' in additional_variables[row['ID']].keys():
            df_infection.at[index, 'household59'] = \
                additional_variables[row['ID']]['household59']
        if 'pregnancytime' in additional_variables[row['ID']].keys():
            df_infection.at[index, 'pregnancytime'] = \
                additional_variables[row['ID']]['pregnancytime']


#risk factors to csv
header = ['ID', 'age_days', 'birthday', 'consultdate', 'infection',\
        'n_household', 'household04','household03','household02', 'household59'\
            ,'visitnursery_child',\
        'visitnursery_house','pregnancytime', 'contacttotal', 'contact04', 'contact59', 'sex']
df_infection.to_csv('/home/andewegs/1_RSV_scripts/tensor_splines/data/Riskfactors2_csv.txt'\
                    , columns=header,index=False)

# %% Figure barplot fraction infected over age in months
fraction_infectd_list, number_ind_list, fraction_infectd_P2, number_ind_P2, \
    fraction_infectd_P3, number_ind_P3 = [], [], [], [], [], []

for i in range(60):
    i += 1  # start at month one
    infected_ind = 0
    total_ind = 0
    for index, row in df_infection.iterrows():
        if round(row['age_months']) == i:
            infected_ind += row['infection']
            total_ind += 1
    if infected_ind == 0 and total_ind > 0:
        fraction_infectd_list.append(0)
    elif infected_ind == 0 and total_ind == 0:
        fraction_infectd_list.append(1.0)
    else:
        fraction_infectd_list.append((infected_ind/total_ind))
    number_ind_list.append(total_ind)
    #Pienter2
    infected_ind = 0
    total_ind = 0
    for index, row in df_infection.loc[df_infection["ID"].str.startswith('P2')].iterrows():
        if round(row['age_months']) == i:
            infected_ind += row['infection']
            total_ind += 1
    if infected_ind == 0 and total_ind > 0:
        fraction_infectd_P2.append(0)
    elif infected_ind == 0 and total_ind == 0:
        fraction_infectd_P2.append(1.0)
    else:
        fraction_infectd_P2.append((infected_ind/total_ind))
    number_ind_P2.append(total_ind)
    #Pienter3
    infected_ind = 0
    total_ind = 0
    for index, row in df_infection.loc[df_infection["ID"].str.startswith('P3')].iterrows():
        if round(row['age_months']) == i:
            infected_ind += row['infection']
            total_ind += 1
    if infected_ind == 0 and total_ind > 0:
        fraction_infectd_P3.append(0)
    elif infected_ind == 0 and total_ind == 0:
        fraction_infectd_P3.append(1.0)
    else:
        fraction_infectd_P3.append((infected_ind/total_ind))
    number_ind_P3.append(total_ind)

number_ind = number_ind_list
fraction_infected = fraction_infectd_list

left, width = 0.1, 1.2
bottom, height = 0.1, 0.6
spacing = 0.005
rect_bar = [left, bottom, width, height]
rect_histx = [left, bottom + height + spacing, width, 0.2]

plt.figure(dpi=300)

ind = np.arange(1, len(fraction_infected)+1)
width_bars = 0.9
ax_bar = plt.axes(rect_bar)  
ax_histx = plt.axes(rect_histx)  

ax_bar.bar(ind, np.ones(len(fraction_infected)), width_bars, color='C0', label='seronegative')
ax_bar.bar(ind, fraction_infected, width_bars, color='C1', label='seropositive')
ax_bar.set_xticks(np.arange(0, len(fraction_infected)+2, 5))
ax_bar.legend(loc='upper right', fontsize=7)
ax_bar.set_xlabel('Age (months)')
ax_bar.set_ylabel('Proportion')

ax_histx.bar(ind, number_ind, width_bars, color='green')
ax_histx.set_xticklabels([])
ax_histx.set_ylabel('Counts')
ax_histx.set_yticks(np.arange(0,60,10))

plt.savefig(path_figures + 'barplot_infected.svg', bbox_inches='tight')
plt.show()

# %% Virologische weekstaten

#Load vir weekstaten
Vir_weekstaten = pd.read_excel('/home/andewegs/3_Documents/Vir_Weekstaten_2020.02.10rsv.xlsx')

Vir_weekstaten['year'] = Vir_weekstaten.apply(lambda row: \
                            int(str(row.JaarWeek)[0:4]), axis = 1)
Vir_weekstaten['week'] = Vir_weekstaten.apply(lambda row: \
                            int(str(row.JaarWeek)[4:6]), axis = 1)

#translate week number to dates each week is a monday    
Vir_weekstaten['date'] = Vir_weekstaten.apply(lambda row: \
                        datetime.datetime.strptime(str(row.JaarWeek) \
                        + '-1', "%Y%W-%w"), axis = 1)

# Plot figure, not in paper
years = Vir_weekstaten['year'].unique()
fig = plt.figure(dpi=300, figsize=(20,12))
for index, item in enumerate(years):
    data = Vir_weekstaten[Vir_weekstaten['year'] == item]
    ax = plt.subplot(4,5,index+1)
    plt.bar(list(data['week']), list(data['AantalRSV']))
    plt.ylim(0,350)
    plt.xlim(0.5,len(data)+0.5)
    plt.title(item, pad = 0)
    if index < 15:
        plt.xticks([])
    else:
        plt.xticks(np.arange(min(data['week'])+9, max(data['week'])+1, 10.0))
    if index not in [0,5,10,15]:
        plt.yticks([])
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)
    else:
        plt.yticks(np.arange(0, 301, 100.0))
        ax.spines['right'].set_visible(False)
    if index in (4,9,14,19):
        ax.spines['right'].set_visible(True)
    
fig.text(-0.006, 0.5, 'number of RSV positive samples', va='center', \
         rotation='vertical', fontsize=20)
fig.text(0.5, -0.01, 'week number', ha='center', fontsize = 20)
plt.tight_layout(h_pad = 0.01, w_pad = 0.0000001)

plt.savefig(path_figures + 'RSV_seasons_histogram.svg', bbox_inches = 'tight')
plt.show() 

# %% Figure lexis diagram infected/uninfected individuals below 1 year
#https://matplotlib.org/gallery/api/two_scales.html

# set sizes figure
left, width = 0, 1
bottom, height = 0, 0.5
spacing = 0.2
rect_lexisP3 = [left, bottom, width, height]
rect_lexisP2 = [left, bottom+0.5+spacing, width, height]

# set cut-off left y-axis
month_y_axis = 12

# plot Figure 2
fig = plt.figure(dpi=300)
ax_lexisP2 = plt.axes(rect_lexisP2)
ax_lexisP3 = plt.axes(rect_lexisP3)

P2_colors = df_infection['infection'].loc[df_infection['age_months'] <= month_y_axis]\
                .loc[df_infection['ID'].str.startswith('P2')].tolist()
color = ['C0', 'C1']
dot_colors = []
for i in P2_colors:
    dot_colors.append('C' +str(i))
counter = 0

for row in df_infection.loc[df_infection['age_months'] <= month_y_axis].itertuples():
    if row.ID.startswith('P2'):
        ax_lexisP2.plot([row.birthday, row.consultdate], [0,row.age_months],\
                        color=color[P2_colors[counter]], linewidth=0.7)
        counter +=1

ax_lexisP2.scatter(df_infection['consultdate'].loc[df_infection['age_months'] \
                                                   <= month_y_axis]\
                .loc[df_infection['ID'].str.startswith('P2')],\
              df_infection['age_months'].loc[df_infection['age_months'] \
                                             <= month_y_axis]\
                .loc[df_infection['ID'].str.startswith('P2')],\
                  marker='.', c=dot_colors, label='consultdate')

    
ax_lexisP2.plot(df_infection['birthday'].loc[df_infection['age_months'] \
                                             <= month_y_axis]\
                .loc[df_infection['ID'].str.startswith('P2')],\
              [0 for i in range(len(df_infection['birthday']\
                                    .loc[df_infection['age_months'] <= month_y_axis]\
                .loc[df_infection['ID'].str.startswith('P2')]))],\
                  'ro', marker='.', c='C2', label='birthday')

ax2_lexisP2 = ax_lexisP2.twinx()
ax2_lexisP2.set_ylabel('RSV counts', color = 'grey')
ax2_lexisP2.set_ylim(0,260)

small_date = pd.Timestamp('2005-02-01') #xmin ax_lexisP2
large_date = pd.Timestamp('2007-07-01') #xmax ax_lexisP2
Vir_data_P2 = Vir_weekstaten.loc[Vir_weekstaten['date'] >= small_date]\
                        .loc[Vir_weekstaten['date'] <= large_date]

ax2_lexisP2.bar(list(Vir_data_P2['date']), list(Vir_data_P2['AantalRSV']), 
                width = 0.8*7, color = 'grey', alpha = 0.7)

years = mdates.YearLocator()
months = mdates.MonthLocator(bymonth = range(2,14,2))
monthsFmt = mdates.DateFormatter('%m') 
yearsFmt = mdates.DateFormatter('\n\n%Y')  # add some space for the year label
ax_lexisP2.xaxis.set_minor_locator(months)
ax_lexisP2.xaxis.set_minor_formatter(monthsFmt)
plt.setp(ax_lexisP2.xaxis.get_minorticklabels(), rotation=90)

ax_lexisP2.set_xlim(np.datetime64(small_date), np.datetime64(large_date))

ax_lexisP2.xaxis.set_major_locator(years)
ax_lexisP2.xaxis.set_major_formatter(yearsFmt)
start, end = ax_lexisP2.get_xlim()
ax_lexisP2.yaxis.set_ticks(np.arange(0, round(max(df_infection["age_months"]\
    .loc[df_infection['ID'].str.startswith('P2')].loc[df_infection['age_months'] <= month_y_axis]))+1, 2))

    
ax_lexisP2.set_ylabel('Age (months)')
ax_lexisP2.set_xlabel('Time (months)')

ax_lexisP2.xaxis.grid(True, which='minor')
ax_lexisP2.yaxis.grid(True, which='major')

P3_colors = df_infection['infection'].loc[df_infection['age_months'] <= month_y_axis]\
                .loc[df_infection['ID'].str.startswith('P3')].tolist()

counter = 0
for row in df_infection.loc[df_infection['age_months'] <= month_y_axis].itertuples():
    if row.ID.startswith('P3'):
        ax_lexisP3.plot([row.birthday, row.consultdate], [0,row.age_months],\
                        color=color[P3_colors[counter]], linewidth=0.7)
        counter +=1
dot_colors=[]
for i in P3_colors:
    dot_colors.append('C' +str(i))
counter = 0

ax_lexisP3.scatter(df_infection['consultdate'].loc[df_infection['age_months'] <= month_y_axis]\
                .loc[df_infection['ID'].str.startswith('P3')],\
              df_infection['age_months'].loc[df_infection['age_months'] <= month_y_axis]\
                .loc[df_infection['ID'].str.startswith('P3')],\
                  marker='.', c=dot_colors, label='consultdate')

    
ax_lexisP3.plot(df_infection['birthday'].loc[df_infection['age_months'] <= month_y_axis]\
                .loc[df_infection['ID'].str.startswith('P3')],\
              [0 for i in range(len(df_infection['birthday'].loc[df_infection['age_months'] <= month_y_axis]\
                .loc[df_infection['ID'].str.startswith('P3')]))],\
                  'ro', marker='.', c='C2', label='birthday')

ax2_lexisP3 = ax_lexisP3.twinx()
ax2_lexisP3.set_ylabel('RSV counts', color = 'grey')
ax2_lexisP3.set_ylim(0,260)

small_date = pd.Timestamp('2015-02-01') #xmin ax_lexisP2
large_date = pd.Timestamp('2017-07-01') #xmax ax_lexisP2
Vir_data_P3 = Vir_weekstaten.loc[Vir_weekstaten['date'] >= small_date]\
                        .loc[Vir_weekstaten['date'] <= large_date]

ax2_lexisP3.bar(list(Vir_data_P3['date']), list(Vir_data_P3['AantalRSV']), 
                width = 0.8*7, color = 'grey', alpha = 0.7)    
    
    
years = mdates.YearLocator()
months = mdates.MonthLocator(bymonth = range(2,14,2))
monthsFmt = mdates.DateFormatter('%m') 
yearsFmt = mdates.DateFormatter('\n\n%Y')  # add some space for the year label
ax_lexisP3.xaxis.set_minor_locator(months)
ax_lexisP3.xaxis.set_minor_formatter(monthsFmt)
plt.setp(ax_lexisP3.xaxis.get_minorticklabels(), rotation=90)

ax_lexisP3.set_xlim(np.datetime64(small_date), np.datetime64(large_date))

ax_lexisP3.xaxis.set_major_locator(years)
ax_lexisP3.xaxis.set_major_formatter(yearsFmt)
start, end = ax_lexisP3.get_xlim()
ax_lexisP3.yaxis.set_ticks(np.arange(0, round(max(df_infection["age_months"]\
    .loc[df_infection['ID'].str.startswith('P3')].loc[df_infection['age_months'] <= month_y_axis]))+1, 2))


ax_lexisP3.set_ylabel('Age (months)')
ax_lexisP3.set_xlabel('Time (months)')

ax_lexisP3.xaxis.grid(True, which='minor')
ax_lexisP3.yaxis.grid(True, which='major')

custom_lines = [Line2D([0], [0], color='C0', lw=0.7),
                Line2D([0], [0], color='C1', lw=0.7)]
ax_lexisP2.legend(custom_lines, ['seronegative', 'seropositive'])

#ax_lexisP2.set_title(RSV_protein)
fig.savefig(path_figures + 'Lexis_diagram_infected_1years'\
            '.svg', format='svg', bbox_inches='tight')
plt.show()
# %% Figure lexis diagram infected/ uninfected individuals below 5 year

# set sizes figure
left, width = 0, 1
bottom, height = 0, 0.5
spacing = 0.2
rect_lexisP3 = [left, bottom, width, height]
rect_lexisP2 = [left, bottom+0.5+spacing, width, height]

month_y_axis = 60

# Plot Supplement Figure
fig = plt.figure(dpi=300)
ax_lexisP2 = plt.axes(rect_lexisP2)
ax_lexisP3 = plt.axes(rect_lexisP3)

P2_colors = df_infection['infection'].loc[df_infection['age_months'] <= month_y_axis]\
                .loc[df_infection['ID'].str.startswith('P2')].tolist()
color = ['C0', 'C1']
dot_colors = []
for i in P2_colors:
    dot_colors.append('C' +str(i))
counter = 0

for row in df_infection.loc[df_infection['age_months'] <= month_y_axis].itertuples():
    if row.ID.startswith('P2'):
        ax_lexisP2.plot([row.birthday, row.consultdate], [0,row.age_months],\
                        color=color[P2_colors[counter]], linewidth=0.7)
        counter +=1

ax_lexisP2.scatter(df_infection['consultdate'].loc[df_infection['age_months'] <= month_y_axis]\
                .loc[df_infection['ID'].str.startswith('P2')],\
              df_infection['age_months'].loc[df_infection['age_months'] <= month_y_axis]\
                .loc[df_infection['ID'].str.startswith('P2')],\
                  marker='.', c=dot_colors, label='consultdate')

    
ax_lexisP2.plot(df_infection['birthday'].loc[df_infection['age_months'] <= month_y_axis]\
                .loc[df_infection['ID'].str.startswith('P2')],\
              [0 for i in range(len(df_infection['birthday'].loc[df_infection['age_months'] <= month_y_axis]\
                .loc[df_infection['ID'].str.startswith('P2')]))],\
                  'ro', marker='.', c='C2', label='birthday')

ax2_lexisP2 = ax_lexisP2.twinx()
ax2_lexisP2.set_ylabel('RSV counts', color = 'grey')
ax2_lexisP2.set_ylim(0,260)

small_date = pd.Timestamp('2001-04-01')

large_date = pd.Timestamp('2008-01-01')
    
Vir_data_P2 = Vir_weekstaten.loc[Vir_weekstaten['date'] >= small_date]\
                        .loc[Vir_weekstaten['date'] <= large_date]

ax2_lexisP2.bar(list(Vir_data_P2['date']), list(Vir_data_P2['AantalRSV']), 
                width = 0.8*7, color = 'grey', alpha = 0.7)    
    
years = mdates.YearLocator()
months = mdates.MonthLocator(bymonth = range(2,14,2))
monthsFmt = mdates.DateFormatter('%m') 
yearsFmt = mdates.DateFormatter('\n\n%Y')  # add some space for the year label
ax_lexisP2.xaxis.set_minor_locator(months)
ax_lexisP2.xaxis.set_minor_formatter(monthsFmt)
plt.setp(ax_lexisP2.xaxis.get_minorticklabels(), rotation=90)

ax_lexisP2.set_xlim(np.datetime64(small_date), np.datetime64(large_date))

ax_lexisP2.xaxis.set_major_locator(years)
ax_lexisP2.xaxis.set_major_formatter(yearsFmt)
start, end = ax_lexisP2.get_xlim()
ax_lexisP2.yaxis.set_ticks(np.arange(0, round(max(df_infection["age_months"]\
    .loc[df_infection['ID'].str.startswith('P2')].loc[df_infection['age_months'] \
                                                      <= month_y_axis]))+1, 5))

ax_lexisP2.set_ylabel('Age (months)')
ax_lexisP2.set_xlabel('Time (months)')

ax_lexisP2.xaxis.grid(True, which='minor')
ax_lexisP2.yaxis.grid(True, which='major')

P3_colors = df_infection['infection'].loc[df_infection['age_months'] <= month_y_axis]\
                .loc[df_infection['ID'].str.startswith('P3')].tolist()

counter = 0
for row in df_infection.loc[df_infection['age_months'] <= month_y_axis].itertuples():
    if row.ID.startswith('P3'):
        ax_lexisP3.plot([row.birthday, row.consultdate], [0,row.age_months],\
                        color=color[P3_colors[counter]], linewidth=0.7)
        counter +=1

dot_colors=[]
for i in P3_colors:
    dot_colors.append('C' +str(i))
counter = 0

ax_lexisP3.scatter(df_infection['consultdate'].loc[df_infection['age_months'] \
                                                   <= month_y_axis]\
                .loc[df_infection['ID'].str.startswith('P3')],\
              df_infection['age_months'].loc[df_infection['age_months'] \
                                             <= month_y_axis]\
                .loc[df_infection['ID'].str.startswith('P3')],\
                  marker='.', c=dot_colors, label='consultdate')

    
ax_lexisP3.plot(df_infection['birthday'].loc[df_infection['age_months'] \
                                             <= month_y_axis]\
                .loc[df_infection['ID'].str.startswith('P3')],\
              [0 for i in range(len(df_infection['birthday']\
                                    .loc[df_infection['age_months'] <= month_y_axis]\
                .loc[df_infection['ID'].str.startswith('P3')]))],\
                  'ro', marker='.', c='C2', label='birthday')

ax2_lexisP3 = ax_lexisP3.twinx()
ax2_lexisP3.set_ylabel('RSV counts', color = 'grey')
ax2_lexisP3.set_ylim(0,260)

small_date = pd.Timestamp('2011-04-01')

large_date = pd.Timestamp('2018-01-01')

Vir_data_P3 = Vir_weekstaten.loc[Vir_weekstaten['date'] >= small_date]\
                        .loc[Vir_weekstaten['date'] <= large_date]

ax2_lexisP3.bar(list(Vir_data_P3['date']), list(Vir_data_P3['AantalRSV']), 
                width = 0.8*7, color = 'grey', alpha = 0.7)        
    
years = mdates.YearLocator()
months = mdates.MonthLocator(bymonth = range(2,14,2))
monthsFmt = mdates.DateFormatter('%m') 
yearsFmt = mdates.DateFormatter('\n\n%Y')  # add some space for the year label
ax_lexisP3.xaxis.set_minor_locator(months)
ax_lexisP3.xaxis.set_minor_formatter(monthsFmt)
plt.setp(ax_lexisP3.xaxis.get_minorticklabels(), rotation=90)

ax_lexisP3.xaxis.set_major_locator(years)
ax_lexisP3.xaxis.set_major_formatter(yearsFmt)
ax_lexisP3.yaxis.set_ticks(np.arange(0, round(max(df_infection["age_months"]\
    .loc[df_infection['ID'].str.startswith('P3')].loc[df_infection['age_months']\
                                                      <= month_y_axis]))+1, 5))

ax_lexisP3.set_xlim(np.datetime64(small_date), np.datetime64(large_date))

ax_lexisP3.set_ylabel('Age (months)')
ax_lexisP3.set_xlabel('Time (months)')

ax_lexisP3.set_yticks(np.arange(0,61,5))

ax_lexisP3.xaxis.grid(True, which='minor')
ax_lexisP3.yaxis.grid(True, which='major')

custom_lines = [Line2D([0], [0], color='C0', lw=0.7),
                Line2D([0], [0], color='C1', lw=0.7)]
ax_lexisP2.legend(custom_lines, ['seronegative', 'seropositive'])


fig.savefig(path_figures + 'Lexis_diagram_infected_5years'\
            '.svg', format='svg', bbox_inches='tight')
plt.show()
