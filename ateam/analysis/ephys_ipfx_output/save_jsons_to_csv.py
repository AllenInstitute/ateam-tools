
# coding: utf-8

# # import libraries

# In[92]:


import utils_ipfx_linear as utils
import os
import csv
import numpy as np
import pandas as pd


# replace all non-numerical values with NaN's

def isnumber(x):
    try:
        float(x)
        return True
    except:
        return False



# # get all directory names

# In[93]:


dir_list= [name for name in os.listdir(".") if os.path.isdir(name)]

# remove the first directory with '.ipynb_checkpoints'
dir_list=dir_list[1:]

print dir_list


# # extract features cell by cell

# In[94]:


# dictionary to store cells dictionary
cells_dict={}

for i in np.arange(len(dir_list)):
    # extract the current dictionary
    current_dict=utils.parse_json_ipfx(dir_list[i])
    
    # save dictionary only if it is not empty
    if current_dict:
        cells_dict[dir_list[i]]=current_dict

# print the cell dictionary
cells_dict


# # save all dictionaries to csv

# In[95]:


# artificially add keys based on successful cell

feature_names=cells_dict['539537044'].keys()
cell_ids=cells_dict.keys()


# In[108]:


# define the total number of cells and total number of features

#keys=cells_dict.keys()

# determine the length based on one dictionary
# number of features is hard-coded
#n_features=len(my_keys)
#n_cells=len(cells_dict)


#print 'Total number of features'
#print n_features
#print
#print 'Total number of cells'
#print n_cells
#print




# # Convert data to dataframe

# In[104]:


# test to create the dataframe from dict

my_data_frame=pd.DataFrame.from_dict(cells_dict,orient='index',)

# add index column
my_data_frame.index.name='spec_id_label'

# show the dataframe
my_data_frame


# In[105]:




my_data_frame=my_data_frame[my_data_frame.applymap(isnumber)]

# replace empty values with NaN's
my_data_frame=my_data_frame.replace(r'\s+', np.nan, regex=True)

print
my_data_frame


# In[107]:


# save the final data to csv

my_data_frame.to_csv('MET_ephys_scaling.csv',sep=',')

