#%% This code focuses on analyzing the radiomic features considering the different variations

# This code performs the radiogenomic analysis

# The radiomic feature position
# Shape Feature: 0 - 13
# First Order Feature: 14 - 31
# GLCM: 32 - 55
# GLDM: 56 - 69
# GLRLM: 70 - 85
# GLSZM: 86 - 101
# NGTDM: 102 - 106

# import numpy as np
import numpy as np
import matplotlib.pyplot as plt

import os
import pandas as pd
from sklearn.preprocessing import StandardScaler, MinMaxScaler
from scipy.stats import linregress

#%% Load existing files ##################################################################################################################

# Load the clinical dataset
file_path = '/Users/hafeng/Documents/Research_Data/UCSF_Dataset/UCSF_PDGM/'
sheet_name = 'preop'

filename = '/preop_glioma_9-20-15--4-1-21-update-surv.xlsx'

# Load the data from CW to ABA columns
gene = pd.read_excel(file_path+filename, sheet_name=sheet_name)
# Remove the last 7 rows using iloc
# gene = gene.iloc[:-10]

# Extract the numerical values into a NumPy array (matrix)
# matrix = gene.to_numpy()
# matrix = matrix[:-7,:]

os_time    = pd.read_excel(file_path+filename, sheet_name=sheet_name, usecols='AD')
os_status  = pd.read_excel(file_path+filename, sheet_name=sheet_name, usecols='AC')
idh_status = pd.read_excel(file_path+filename, sheet_name=sheet_name, usecols='Q')

#%% Load the patient shape features

data_loc = '/Users/hafeng/Documents/Postdoc_Research_Meetings/TCGA/UCSF_PDGM/radiogenomic_whole/'
file_numbers = np.loadtxt(data_loc + 'patient_list_shape.txt')
ce1_shape    = np.load(data_loc + 'patient_shape_ce1.npy')
ce2_shape    = np.load(data_loc + 'patient_shape_ce2.npy')
ce4_shape    = np.load(data_loc + 'patient_shape_ce4.npy')
fl1_shape    = np.load(data_loc + 'patient_shape_fl1.npy')
fl2_shape    = np.load(data_loc + 'patient_shape_fl2.npy')
fl4_shape    = np.load(data_loc + 'patient_shape_fl4.npy')
adc1_shape   = np.load(data_loc + 'patient_shape_adc1.npy')
adc2_shape   = np.load(data_loc + 'patient_shape_adc2.npy')
adc4_shape   = np.load(data_loc + 'patient_shape_adc4.npy')





#%% Check for the value of the existing files ##############################################################################################

# Specify the folder where your files are located
# radiomic_path = '/Users/hafeng/Documents/Postdoc_Research_Meetings/TCGA/UCSF_PDGM/radiogenomic_tumor/'
radiomic_path = '/Users/hafeng/Documents/Postdoc_Research_Meetings/TCGA/UCSF_PDGM/radiogenomic_whole/'

# Initialize an empty list to store the integers
file_numbers1 = []
file_numbers2 = []
file_numbers4 = []
file_numbers0 = []

# Loop through all the files in the folder
for filename in os.listdir(radiomic_path):
    
    # if filename.endswith('_adc_combined_iso_256.xlsx') and filename[:4].isdigit():
    if filename.endswith('_label_1_shell_19_CE_features_multi8_norm_update_l1.xlsx') and filename[:4].isdigit():
        # Extract the first 4 characters and convert to an integer
        file_number = int(filename[:4])
        # Append the integer to the list
        file_numbers1.append(file_number)

    if filename.endswith('_label_2_shell_19_CE_features_multi8_norm_update_l1.xlsx') and filename[:4].isdigit():
        # Extract the first 4 characters and convert to an integer
        file_number = int(filename[:4])
        # Append the integer to the list
        file_numbers2.append(file_number)

    if filename.endswith('_label_4_shell_19_CE_features_multi8_norm_update_l1.xlsx') and filename[:4].isdigit():
        # Extract the first 4 characters and convert to an integer
        file_number = int(filename[:4])
        # Append the integer to the list
        file_numbers4.append(file_number)
        
    if filename.endswith('_label_0_shell_4_CE_features_multi8_norm_update_l1.xlsx') and filename[:4].isdigit():
        # Extract the first 4 characters and convert to an integer
        file_number = int(filename[:4])
        # Append the integer to the list
        file_numbers0.append(file_number)

# Display the resulting list of integers
file_numbers_unique = list(set(file_numbers1) & set(file_numbers2) & set(file_numbers4) & set(file_numbers0))
file_numbers_unique.sort()

print(file_numbers_unique)
print(len(file_numbers_unique))

common_elements = file_numbers_unique

#%% Check for the value of the existing files ##############################################################################################

# # Specify the folder where your files are located
# radiomic_path1 = '/Users/hafeng/Documents/Postdoc_Research_Meetings/TCGA/UCSF_PDGM/radiogenomic_shell_whole/'

# # Initialize an empty list to store the integers
# file_numbers1 = []
# file_numbers2 = []
# file_numbers4 = []
# file_numbers0 = []

# # Loop through all the files in the folder
# for filename in os.listdir(radiomic_path1):
    
#     # if filename.endswith('_adc_combined_iso_256.xlsx') and filename[:4].isdigit():
#     if filename.endswith('_label_1_shell_19_CE_features_multi8_norm_update.xlsx') and filename[:4].isdigit():
#         # Extract the first 4 characters and convert to an integer
#         file_number = int(filename[:4])
#         # Append the integer to the list
#         file_numbers1.append(file_number)

#     if filename.endswith('_label_2_shell_19_CE_features_multi8_norm_update.xlsx') and filename[:4].isdigit():
#         # Extract the first 4 characters and convert to an integer
#         file_number = int(filename[:4])
#         # Append the integer to the list
#         file_numbers2.append(file_number)

#     if filename.endswith('_label_4_shell_19_CE_features_multi8_norm_update.xlsx') and filename[:4].isdigit():
#         # Extract the first 4 characters and convert to an integer
#         file_number = int(filename[:4])
#         # Append the integer to the list
#         file_numbers4.append(file_number)
        
#     # if filename.endswith('_label_0_shell_4_CE_features_multi8_norm_update.xlsx') and filename[:4].isdigit():
#     #     # Extract the first 4 characters and convert to an integer
#     #     file_number = int(filename[:4])
#     #     # Append the integer to the list
#     #     file_numbers0.append(file_number)

# # Display the resulting list of integers
# file_numbers_unique1 = list(set(file_numbers1) & set(file_numbers2) & set(file_numbers4))
# file_numbers_unique1.sort()

# print(file_numbers_unique1)
# print(len(file_numbers_unique1))

# # Find the common list
# common_number    = list(set(file_numbers_unique) & set(file_numbers_unique1))
# common_elements  = [num for num in common_number]
# print(common_elements)  # Output: [4, 5]
# print(len(common_elements))

#%% Extract the genomic matrix with the given unique values ######################################################################################

gene_filter = gene[gene['Assigned ID'].isin(file_numbers_unique)]
gene_matrix = gene_filter.to_numpy()

# Add Clinical Information
sex = gene_filter['Sex']
Age = gene_filter['Age at MRI'] / 100

sex_status    = sex.map(lambda x: -1 if pd.isna(x) else (0 if str(x).lower() == 'f' else 1)).to_numpy()
sex_status_df = sex.map(lambda x: -1 if pd.isna(x) else (0 if str(x).lower() == 'f' else 1))
print(sex_status_df)

age = Age.to_numpy() / 100
print(age)
# os_time   = gene_fitler['']
# os_status = pd.read_excel(file_path+filename, sheet_name=sheet_name, usecols='AC')



#%%
skip_row  = 23
data_dir  = radiomic_path
stacked_data_ce0  = []
stacked_data_fl0  = []
stacked_data_ce1  = []
stacked_data_fl1  = []
stacked_data_ce2  = []
stacked_data_fl2  = []
stacked_data_ce4  = []
stacked_data_fl4  = []
stacked_data_adc0 = []
stacked_data_adc1 = []
stacked_data_adc2 = []
stacked_data_adc4 = []

import re

def extract_shell_num(filename):
    # Extract the number after "shell_"
    match = re.search(r'shell_(\d+)', filename)
    return int(match.group(1)) if match else -1

for i in range(len(common_elements)):
    idx = file_numbers_unique[i]
    # if idx > 520:
        # idx = 443
    folder_name = f'{idx:04d}'
    print(f'Working for Patient {folder_name}')
 
    pattern1 = f"{idx:04d}_label_1_shell_"
    pattern2 = f"{idx:04d}_label_2_shell_"
    pattern4 = f"{idx:04d}_label_4_shell_"
    pattern0 = f"{idx:04d}_label_0_shell_"
    
    files_ce1  = [f for f in os.listdir(data_dir) if f.startswith(pattern1) and f.endswith('_CE_features_multi8_norm_update_l1.xlsx')]
    files_ce2  = [f for f in os.listdir(data_dir) if f.startswith(pattern2) and f.endswith('_CE_features_multi8_norm_update_l1.xlsx')]
    files_ce4  = [f for f in os.listdir(data_dir) if f.startswith(pattern4) and f.endswith('_CE_features_multi8_norm_update_l1.xlsx')]
    files_ce0  = [f for f in os.listdir(data_dir) if f.startswith(pattern0) and f.endswith('_CE_features_multi8_norm_update_l1.xlsx')]
    files_fl1  = [f for f in os.listdir(data_dir) if f.startswith(pattern1) and f.endswith('_FL_features_multi8_norm_update_l1.xlsx')]
    files_fl2  = [f for f in os.listdir(data_dir) if f.startswith(pattern2) and f.endswith('_FL_features_multi8_norm_update_l1.xlsx')]
    files_fl4  = [f for f in os.listdir(data_dir) if f.startswith(pattern4) and f.endswith('_FL_features_multi8_norm_update_l1.xlsx')]
    files_fl0  = [f for f in os.listdir(data_dir) if f.startswith(pattern0) and f.endswith('_FL_features_multi8_norm_update_l1.xlsx')]
    files_adc1 = [f for f in os.listdir(data_dir) if f.startswith(pattern1) and f.endswith('_ADC_features_multi8_norm_update_l1.xlsx')]
    files_adc2 = [f for f in os.listdir(data_dir) if f.startswith(pattern2) and f.endswith('_ADC_features_multi8_norm_update_l1.xlsx')]
    files_adc4 = [f for f in os.listdir(data_dir) if f.startswith(pattern4) and f.endswith('_ADC_features_multi8_norm_update_l1.xlsx')]
    files_adc0 = [f for f in os.listdir(data_dir) if f.startswith(pattern0) and f.endswith('_ADC_features_multi8_norm_update_l1.xlsx')]

    files_ce0.sort(key=extract_shell_num)
    files_ce1.sort(key=extract_shell_num)
    files_ce2.sort(key=extract_shell_num)
    files_ce4.sort(key=extract_shell_num)
    files_fl0.sort(key=extract_shell_num)
    files_fl1.sort(key=extract_shell_num)
    files_fl2.sort(key=extract_shell_num)
    files_fl4.sort(key=extract_shell_num)
    files_adc0.sort(key=extract_shell_num)
    files_adc1.sort(key=extract_shell_num)
    files_adc2.sort(key=extract_shell_num)
    files_adc4.sort(key=extract_shell_num)

    data_list_ce1  = []
    data_list_fl1  = []
    data_list_ce2  = []
    data_list_fl2  = []
    data_list_ce4  = []
    data_list_fl4  = []
    data_list_ce0  = []
    data_list_fl0  = []
    data_list_adc1  = []
    data_list_adc2  = []
    data_list_adc4  = []
    data_list_adc0  = []
        
    for file in files_ce0:
        file_path = os.path.join(radiomic_path, file)
        df = pd.read_excel(file_path, sheet_name='Sheet1', skiprows=skip_row, usecols='B', header=None)  # Adjust for your Excel file format if needed
        # Convert the dataframe to a 1D array (or whatever structure you need)
        data_array = df.values.flatten()  # Flatten into a 1D array
        data_list_ce0.append(data_array)
    # Stack all data arrays into a matrix
    stacked_matrix_ce0 = np.vstack(data_list_ce0)
    stacked_data_ce0.append(stacked_matrix_ce0)
    
    for file in files_fl0:
        file_path = os.path.join(radiomic_path, file)
        df = pd.read_excel(file_path, sheet_name='Sheet1', skiprows=skip_row, usecols='B', header=None)  # Adjust for your Excel file format if needed
        # Convert the dataframe to a 1D array (or whatever structure you need)
        data_array = df.values.flatten()  # Flatten into a 1D array
        data_list_fl0.append(data_array)
    # Stack all data arrays into a matrix
    stacked_matrix_fl0 = np.vstack(data_list_fl0)
    stacked_data_fl0.append(stacked_matrix_fl0)
    
    for file in files_ce1:
        file_path = os.path.join(radiomic_path, file)
        df = pd.read_excel(file_path, sheet_name='Sheet1', skiprows=skip_row, usecols='B', header=None)  # Adjust for your Excel file format if needed
        # Convert the dataframe to a 1D array (or whatever structure you need)
        data_array = df.values.flatten()  # Flatten into a 1D array
        data_list_ce1.append(data_array)
    # Stack all data arrays into a matrix
    stacked_matrix_ce1 = np.vstack(data_list_ce1)
    stacked_data_ce1.append(stacked_matrix_ce1)

    for file in files_ce2:
        file_path = os.path.join(radiomic_path, file)
        df = pd.read_excel(file_path, sheet_name='Sheet1', skiprows=skip_row, usecols='B', header=None)  # Adjust for your Excel file format if needed
        # Convert the dataframe to a 1D array (or whatever structure you need)
        data_array = df.values.flatten()  # Flatten into a 1D array
        data_list_ce2.append(data_array)
    # Stack all data arrays into a matrix
    stacked_matrix_ce2 = np.vstack(data_list_ce2)
    stacked_data_ce2.append(stacked_matrix_ce2)

    for file in files_ce4:
        file_path = os.path.join(radiomic_path, file)
        df = pd.read_excel(file_path, sheet_name='Sheet1', skiprows=skip_row, usecols='B', header=None)  # Adjust for your Excel file format if needed
        # Convert the dataframe to a 1D array (or whatever structure you need)
        data_array = df.values.flatten()  # Flatten into a 1D array
        data_list_ce4.append(data_array)
    # Stack all data arrays into a matrix
    stacked_matrix_ce4 = np.vstack(data_list_ce4)
    stacked_data_ce4.append(stacked_matrix_ce4)
        
    for file in files_fl1:
        file_path = os.path.join(radiomic_path, file)
        df = pd.read_excel(file_path, sheet_name='Sheet1', skiprows=skip_row, usecols='B', header=None)  # Adjust for your Excel file format if needed
        # Convert the dataframe to a 1D array (or whatever structure you need)
        data_array = df.values.flatten()  # Flatten into a 1D array
        data_list_fl1.append(data_array)
    # Stack all data arrays into a matrix
    stacked_matrix_fl1 = np.vstack(data_list_fl1)
    stacked_data_fl1.append(stacked_matrix_fl1)

    for file in files_fl2:
        file_path = os.path.join(radiomic_path, file)
        df = pd.read_excel(file_path, sheet_name='Sheet1', skiprows=skip_row, usecols='B', header=None)  # Adjust for your Excel file format if needed
        # Convert the dataframe to a 1D array (or whatever structure you need)
        data_array = df.values.flatten()  # Flatten into a 1D array
        data_list_fl2.append(data_array)
    # Stack all data arrays into a matrix
    stacked_matrix_fl2 = np.vstack(data_list_fl2)
    stacked_data_fl2.append(stacked_matrix_fl2)

    for file in files_fl4:
        file_path = os.path.join(radiomic_path, file)
        df = pd.read_excel(file_path, sheet_name='Sheet1', skiprows=skip_row, usecols='B', header=None)  # Adjust for your Excel file format if needed
        # Convert the dataframe to a 1D array (or whatever structure you need)
        data_array = df.values.flatten()  # Flatten into a 1D array
        data_list_fl4.append(data_array)
    # Stack all data arrays into a matrix
    stacked_matrix_fl4 = np.vstack(data_list_fl4)
    stacked_data_fl4.append(stacked_matrix_fl4)
    
    for file in files_adc0:
        file_path = os.path.join(data_dir, file)
        df = pd.read_excel(file_path, sheet_name='Sheet1', skiprows=skip_row, usecols='B', header=None)  # Adjust for your Excel file format if needed
        # Convert the dataframe to a 1D array (or whatever structure you need)
        data_array = df.values.flatten()  # Flatten into a 1D array
        data_list_adc0.append(data_array)
    # Stack all data arrays into a matrix
    stacked_matrix_adc0 = np.vstack(data_list_adc0)
    stacked_data_adc0.append(stacked_matrix_adc0)
    
    for file in files_adc1:
        file_path = os.path.join(data_dir, file)
        df = pd.read_excel(file_path, sheet_name='Sheet1', skiprows=skip_row, usecols='B', header=None)  # Adjust for your Excel file format if needed
        # Convert the dataframe to a 1D array (or whatever structure you need)
        data_array = df.values.flatten()  # Flatten into a 1D array
        data_list_adc1.append(data_array)
    # Stack all data arrays into a matrix
    stacked_matrix_adc1 = np.vstack(data_list_adc1)
    stacked_data_adc1.append(stacked_matrix_adc1)

    for file in files_adc2:
        file_path = os.path.join(data_dir, file)
        df = pd.read_excel(file_path, sheet_name='Sheet1', skiprows=skip_row, usecols='B', header=None)  # Adjust for your Excel file format if needed
        # Convert the dataframe to a 1D array (or whatever structure you need)
        data_array = df.values.flatten()  # Flatten into a 1D array
        data_list_adc2.append(data_array)
    # Stack all data arrays into a matrix
    stacked_matrix_adc2 = np.vstack(data_list_adc2)
    stacked_data_adc2.append(stacked_matrix_adc2)

    for file in files_adc4:
        file_path = os.path.join(data_dir, file)
        df = pd.read_excel(file_path, sheet_name='Sheet1', skiprows=skip_row, usecols='B', header=None)  # Adjust for your Excel file format if needed
        # Convert the dataframe to a 1D array (or whatever structure you need)
        data_array = df.values.flatten()  # Flatten into a 1D array
        data_list_adc4.append(data_array)
    # Stack all data arrays into a matrix
    stacked_matrix_adc4 = np.vstack(data_list_adc4)
    stacked_data_adc4.append(stacked_matrix_adc4)
        
#%% Uniformly reshape the features #############################################################
from collections import Counter

# Example: If stacked_data_ce1t is a list of matrices
valid_indices0 = [i for i, matrix in enumerate(stacked_data_ce0) if matrix.shape[0] == 20]
valid_indices1 = [i for i, matrix in enumerate(stacked_data_ce1) if matrix.shape[0] == 20]
valid_indices2 = [i for i, matrix in enumerate(stacked_data_ce2) if matrix.shape[0] == 20]
valid_indices4 = [i for i, matrix in enumerate(stacked_data_ce4) if matrix.shape[0] == 20]
valid_indices  = set(valid_indices1) & set(valid_indices2) & set(valid_indices4) & set(valid_indices0)

filtered_matrices = [stacked_data_ce1[i] for i in valid_indices]
# Convert to numpy array if needed
filtered_data = np.array([matrix.flatten() for matrix in filtered_matrices])

print("Indices of matrices with 20 rows:", valid_indices)
print("Filtered data shape:", filtered_data.shape)

common_indices = [common_elements[i] for i in list(valid_indices)]
#common_indices = common_elements
print("Real Indices of amtrices with 20 rows:", common_indices)

#%% Obtain the shape features from common_indices

common_final = [x for x in common_indices if x in file_numbers]

select_fea   = [i for i, pid in enumerate(file_numbers) if pid in common_final]
ce1_shape_final  = ce1_shape[:,select_fea]
ce2_shape_final  = ce2_shape[:,select_fea]
ce4_shape_final  = ce4_shape[:,select_fea]
fl1_shape_final  = fl1_shape[:,select_fea]
fl2_shape_final  = fl2_shape[:,select_fea]
fl4_shape_final  = fl4_shape[:,select_fea]
adc1_shape_final = adc1_shape[:,select_fea]
adc2_shape_final = adc2_shape[:,select_fea]
adc4_shape_final = adc4_shape[:,select_fea]

ce_shape_final   = np.concatenate((ce1_shape_final,  ce2_shape_final, ce4_shape_final), axis=0)
fl_shape_final   = np.concatenate((fl1_shape_final,  fl2_shape_final, fl4_shape_final), axis=0)
adc_shape_final  = np.concatenate((adc1_shape_final, adc2_shape_final, adc4_shape_final), axis=0)

#%%
from collections import Counter

# select_indices = [i for i, pid in enumerate(common_indices) if pid in common_final]
select_indices = [i for i in valid_indices if common_elements[i] in common_final]
# select_indices = valid_indices

stacked_data_ce00  = [stacked_data_ce0[i]  for i in select_indices]
stacked_data_ce11  = [stacked_data_ce1[i]  for i in select_indices]
stacked_data_ce22  = [stacked_data_ce2[i]  for i in select_indices]
stacked_data_ce44  = [stacked_data_ce4[i]  for i in select_indices]
stacked_data_fl00  = [stacked_data_fl0[i]  for i in select_indices]
stacked_data_fl11  = [stacked_data_fl1[i]  for i in select_indices]
stacked_data_fl22  = [stacked_data_fl2[i]  for i in select_indices]
stacked_data_fl44  = [stacked_data_fl4[i]  for i in select_indices]
stacked_data_adc00 = [stacked_data_adc0[i] for i in select_indices]
stacked_data_adc11 = [stacked_data_adc1[i] for i in select_indices]
stacked_data_adc22 = [stacked_data_adc2[i] for i in select_indices]
stacked_data_adc44 = [stacked_data_adc4[i] for i in select_indices]

length = []
data_test = stacked_data_ce4

for i in range(len(data_test)):
    length.append(data_test[i].shape[1])

length_new = np.sort(length)
count_dict = Counter(length_new)
print(count_dict)

#.reshape(-1)
stacked_data_ce0  = np.array([matrix.reshape(-1) for matrix in stacked_data_ce00])
stacked_data_fl0  = np.array([matrix.reshape(-1) for matrix in stacked_data_fl00])
stacked_data_ce1  = np.array([matrix.reshape(-1) for matrix in stacked_data_ce11])
stacked_data_fl1  = np.array([matrix.reshape(-1) for matrix in stacked_data_fl11])
stacked_data_ce2  = np.array([matrix.reshape(-1) for matrix in stacked_data_ce22])
stacked_data_fl2  = np.array([matrix.reshape(-1) for matrix in stacked_data_fl22])
stacked_data_ce4  = np.array([matrix.reshape(-1) for matrix in stacked_data_ce44])
stacked_data_fl4  = np.array([matrix.reshape(-1) for matrix in stacked_data_fl44])
stacked_data_adc0 = np.array([matrix.reshape(-1) for matrix in stacked_data_adc00])
stacked_data_adc1 = np.array([matrix.reshape(-1) for matrix in stacked_data_adc11])
stacked_data_adc2 = np.array([matrix.reshape(-1) for matrix in stacked_data_adc22])
stacked_data_adc4 = np.array([matrix.reshape(-1) for matrix in stacked_data_adc44])

# stacked_data_ce0  = np.array([matrix.flatten() for matrix in stacked_data_ce00])
# stacked_data_fl0  = np.array([matrix.flatten() for matrix in stacked_data_fl00])
# stacked_data_ce1  = np.array([matrix.flatten() for matrix in stacked_data_ce11])
# stacked_data_fl1  = np.array([matrix.flatten() for matrix in stacked_data_fl11])
# stacked_data_ce2  = np.array([matrix.flatten() for matrix in stacked_data_ce22])
# stacked_data_fl2  = np.array([matrix.flatten() for matrix in stacked_data_fl22])
# stacked_data_ce4  = np.array([matrix.flatten() for matrix in stacked_data_ce44])
# stacked_data_fl4  = np.array([matrix.flatten() for matrix in stacked_data_fl44])
# stacked_data_adc0 = np.array([matrix.flatten() for matrix in stacked_data_adc00])
# stacked_data_adc1 = np.array([matrix.flatten() for matrix in stacked_data_adc11])
# stacked_data_adc2 = np.array([matrix.flatten() for matrix in stacked_data_adc22])
# stacked_data_adc4 = np.array([matrix.flatten() for matrix in stacked_data_adc44])



#%% New stacks
stacked_data_ce0_new = stacked_data_ce0
stacked_data_ce1_new = stacked_data_ce1
stacked_data_ce2_new = stacked_data_ce2
stacked_data_ce4_new = stacked_data_ce4

stacked_data_fl0_new = stacked_data_fl0
stacked_data_fl1_new = stacked_data_fl1
stacked_data_fl2_new = stacked_data_fl2
stacked_data_fl4_new = stacked_data_fl4

stacked_data_adc0_new = stacked_data_adc0
stacked_data_adc1_new = stacked_data_adc1
stacked_data_adc2_new = stacked_data_adc2
stacked_data_adc4_new = stacked_data_adc4

# stacked_data_ce0_new = np.concatenate((stacked_data_ce0, ce_shape_final.T), axis=1)
# stacked_data_ce1_new = np.concatenate((stacked_data_ce1, ce_shape_final.T), axis=1)
# stacked_data_ce2_new = np.concatenate((stacked_data_ce2, ce_shape_final.T), axis=1)
# stacked_data_ce4_new = np.concatenate((stacked_data_ce4, ce_shape_final.T), axis=1)

# stacked_data_fl0_new = np.concatenate((stacked_data_fl0, fl_shape_final.T), axis=1)
# stacked_data_fl1_new = np.concatenate((stacked_data_fl1, fl_shape_final.T), axis=1)
# stacked_data_fl2_new = np.concatenate((stacked_data_fl2, fl_shape_final.T), axis=1)
# stacked_data_fl4_new = np.concatenate((stacked_data_fl4, fl_shape_final.T), axis=1)

# stacked_data_adc0_new = np.concatenate((stacked_data_adc0, adc_shape_final.T), axis=1)
# stacked_data_adc1_new = np.concatenate((stacked_data_adc1, adc_shape_final.T), axis=1)
# stacked_data_adc2_new = np.concatenate((stacked_data_adc2, adc_shape_final.T), axis=1)
# stacked_data_adc4_new = np.concatenate((stacked_data_adc4, adc_shape_final.T), axis=1)

#%% Reform the matrix ########################################################################
diff1 = list(set(valid_indices) - set(select_indices))
print("In list1 but not in list2:", diff1)

ce_df    = pd.read_excel(file_path, sheet_name='Sheet1', skiprows=skip_row, usecols='A', header=None)
ce_list  = ce_df.iloc[:,0].tolist()

#%% Normalize the radiomic data ##############################################################

# Apply the Min-Max normalization
def min_max_norm(matrix):
    # row_min = np.min(matrix, axis=1)  # Minimum of each row
    # row_max = np.max(matrix, axis=1)  # Maximum of each row
    
    row_min = np.min(matrix, axis=0, keepdims=True)  # Minimum of each row
    row_max = np.max(matrix, axis=0, keepdims=True)  # Maximum of each row
    
    normalized_matrix = (matrix - row_min) / (row_max - row_min)
    
    return normalized_matrix

def min_max_norm_df(df):
    
    # Function for min-max normalization
    row_min = df.min(axis=1)  # Minimum of each row
    row_max = df.max(axis=1)  # Maximum of each row
    
    # Perform min-max normalization for each row
    normalized_df = (df.subtract(row_min, axis=0)).div(row_max - row_min, axis=0)
    
    return normalized_df
# 
min_max_scaler = MinMaxScaler()
# min_max_scaler = StandardScaler()
# 
# ce0_mat_norm  = min_max_scaler.fit_transform(stacked_data_ce0_new)
# fl0_mat_norm  = min_max_scaler.fit_transform(stacked_data_fl0_new)
# ce1_mat_norm  = min_max_scaler.fit_transform(stacked_data_ce1_new)
# fl1_mat_norm  = min_max_scaler.fit_transform(stacked_data_fl1_new)
# ce2_mat_norm  = min_max_scaler.fit_transform(stacked_data_ce2_new)
# fl2_mat_norm  = min_max_scaler.fit_transform(stacked_data_fl2_new)
# ce4_mat_norm  = min_max_scaler.fit_transform(stacked_data_ce4_new)
# fl4_mat_norm  = min_max_scaler.fit_transform(stacked_data_fl4_new)
# adc0_mat_norm = min_max_scaler.fit_transform(stacked_data_adc0_new)
# adc1_mat_norm = min_max_scaler.fit_transform(stacked_data_adc1_new)
# adc2_mat_norm = min_max_scaler.fit_transform(stacked_data_adc2_new)
# adc4_mat_norm = min_max_scaler.fit_transform(stacked_data_adc4_new)

ce0_mat_norm  = stacked_data_ce0_new
fl0_mat_norm  = stacked_data_fl0_new
ce1_mat_norm  = stacked_data_ce1_new
fl1_mat_norm  = stacked_data_fl1_new
ce2_mat_norm  = stacked_data_ce2_new
fl2_mat_norm  = stacked_data_fl2_new
ce4_mat_norm  = stacked_data_ce4_new
fl4_mat_norm  = stacked_data_fl4_new
adc0_mat_norm = stacked_data_adc0_new
adc1_mat_norm = stacked_data_adc1_new
adc2_mat_norm = stacked_data_adc2_new
adc4_mat_norm = stacked_data_adc4_new

ce_shape_final_norm  = min_max_scaler.fit_transform(ce_shape_final.T)
fl_shape_final_norm  = min_max_scaler.fit_transform(fl_shape_final.T)
adc_shape_final_norm = min_max_scaler.fit_transform(adc_shape_final.T)

#%%
shape_names = [f"Shape_Feature{i}" for i in range(13)]
# radiomic_feature_online  = ce_list * 20
radiomic_feature_online1 = ce_list * 20
shape_feature_online     = shape_names * 3

radiomic_feature_online = []
for suffix_num in range(20):  # S0 to S19
    suffixed_features = [f"{feature}_S{suffix_num}" for feature in ce_list]
    radiomic_feature_online.extend(suffixed_features)

ce0_radiomic  = pd.DataFrame(ce0_mat_norm,   columns=radiomic_feature_online)
fl0_radiomic  = pd.DataFrame(fl0_mat_norm,   columns=radiomic_feature_online)
ce1_radiomic  = pd.DataFrame(ce1_mat_norm,   columns=radiomic_feature_online)
fl1_radiomic  = pd.DataFrame(fl1_mat_norm,   columns=radiomic_feature_online)
ce2_radiomic  = pd.DataFrame(ce2_mat_norm,   columns=radiomic_feature_online)
fl2_radiomic  = pd.DataFrame(fl2_mat_norm,   columns=radiomic_feature_online)
ce4_radiomic  = pd.DataFrame(ce4_mat_norm,   columns=radiomic_feature_online)
fl4_radiomic  = pd.DataFrame(fl4_mat_norm,   columns=radiomic_feature_online)
adc0_radiomic = pd.DataFrame(adc0_mat_norm,  columns=radiomic_feature_online)
adc1_radiomic = pd.DataFrame(adc1_mat_norm,  columns=radiomic_feature_online)
adc2_radiomic = pd.DataFrame(adc2_mat_norm,  columns=radiomic_feature_online)
adc4_radiomic = pd.DataFrame(adc4_mat_norm,  columns=radiomic_feature_online)

ce_shape_radiomic  = pd.DataFrame(ce_shape_final_norm,  columns=shape_feature_online)
fl_shape_radiomic  = pd.DataFrame(fl_shape_final_norm,  columns=shape_feature_online)
adc_shape_radiomic = pd.DataFrame(adc_shape_final_norm, columns=shape_feature_online)

#%% Define the function to remove constant in the dataframe
from sklearn.feature_selection import VarianceThreshold
# from statsmodels.stats.outliers_influence import variance_inflation_factor

def const_remove(df):
    # Find columns with only one unique value
    constant_columns = df.columns[df.nunique() <= 1]
    print(constant_columns)
    df_cleaned = df.drop(columns=constant_columns)
    return df_cleaned

def var_remove(df, threshold):
    low_variance_columns = df.columns[df.var() < threshold]
    print(low_variance_columns)
    
    if len(low_variance_columns)>0:
        df.drop(columns=low_variance_columns, inplace=True)
    
    return df

#%% Remove constant and collineary parameters
ce1_radiomic_clean   = const_remove(ce1_radiomic)
fl1_radiomic_clean   = const_remove(fl1_radiomic)
ce2_radiomic_clean   = const_remove(ce2_radiomic)
fl2_radiomic_clean   = const_remove(fl2_radiomic)
ce4_radiomic_clean   = const_remove(ce4_radiomic)
fl4_radiomic_clean   = const_remove(fl4_radiomic)
ce0_radiomic_clean   = const_remove(ce0_radiomic)
fl0_radiomic_clean   = const_remove(fl0_radiomic)
adc0_radiomic_clean  = const_remove(adc0_radiomic)
adc1_radiomic_clean  = const_remove(adc1_radiomic)
adc2_radiomic_clean  = const_remove(adc2_radiomic)
adc4_radiomic_clean  = const_remove(adc4_radiomic)

ce_shape_clean   = const_remove(ce_shape_radiomic)
fl_shape_clean   = const_remove(fl_shape_radiomic)
adc_shape_clean  = const_remove(adc_shape_radiomic)

#%% MGMT status check
gene_filter = gene[gene['Assigned ID'].isin(file_numbers_unique)]
gene_matrix = gene_filter.to_numpy()

# gene_matrix = gene_filter.to_numpy()
MGMT = gene_filter['MGMT']
gene_filter_1 = gene_filter
gene_filter_1['MGMT'] = gene_filter_1['MGMT'].map({'positive': 1, 'negative': 0})
gene_filter_1['EGFR'] = gene_filter_1['EGFR'].map(lambda x: 0 if x == 'negative' else 1)
gene_filter_1['PTEN'] = gene_filter_1['PTEN'].map(lambda x: 0 if x == 'negative' else 1)

threshold = 0.05
gene_filter_1['TP53'] = gene_filter_1['TP53'].map(lambda x: 0 if x == 'negative' else 1 if x in ['positive', 'rare'] else 1 if isinstance(x, (int,float)) and x >= threshold else 0 if isinstance(x, (int,float)) and x < threshold else None)

# Apply the Min-Max normalization
def min_max_final(matrix):
    # row_min = np.min(matrix, axis=1)  # Minimum of each row
    # row_max = np.max(matrix, axis=1)  # Maximum of each row
    
    row_min = np.min(matrix, axis=0, keepdims=True)  # Minimum of each row
    row_max = np.max(matrix, axis=0, keepdims=True)  # Maximum of each row
    normalized_matrix = (matrix - row_min) / (row_max - row_min)
    
    return normalized_matrix

#%% Extract for Certain radiomics

# First Order Feature: 0 - 17
# GLCM: 18 - 41
# GLDM: 42 - 55
# GLRLM: 56 - 71
# GLSZM: 72 - 87
# NGTDM: 88 - 92

n_features_per_region = 93
n_keep_first = 18
n_keep_glcm  = 24
n_keep_gldm  = 14
n_keep_glrlm = 16
n_keep_glszm = 16
n_keep_ngtdm = 5

start_first = 0
start_glcm  = 18
start_gldm  = 42
start_glrlm = 56
start_glszm = 72
start_ngtdm = 88

keep_indices_first = []
keep_indices_glcm  = []
keep_indices_gldm  = []
keep_indices_glrlm = []
keep_indices_glszm = []
keep_indices_ngtdm = []

for region_idx in range(20):
    s_first = start_first + region_idx * n_features_per_region
    keep_indices_first.extend(range(s_first, s_first + n_keep_first))
    
    s_glcm = start_glcm + region_idx * n_features_per_region
    keep_indices_glcm.extend(range(s_glcm, s_glcm + n_keep_glcm))
    
    s_gldm = start_gldm + region_idx * n_features_per_region
    keep_indices_gldm.extend(range(s_gldm, s_gldm + n_keep_gldm))
    
    s_glrlm = start_glrlm + region_idx * n_features_per_region
    keep_indices_glrlm.extend(range(s_glrlm, s_glrlm + n_keep_glrlm))
    
    s_glszm = start_glszm + region_idx * n_features_per_region
    keep_indices_glszm.extend(range(s_glszm, s_glszm + n_keep_glszm))
    
    s_ngtdm = start_ngtdm + region_idx * n_features_per_region
    keep_indices_ngtdm.extend(range(s_ngtdm, s_ngtdm + n_keep_ngtdm))

#%% 
ce1_radiomic_clean_first = ce1_radiomic_clean.iloc[:, keep_indices_first]
ce2_radiomic_clean_first = ce2_radiomic_clean.iloc[:, keep_indices_first]
ce4_radiomic_clean_first = ce4_radiomic_clean.iloc[:, keep_indices_first]
ce1_radiomic_clean_glcm  = ce1_radiomic_clean.iloc[:, keep_indices_glcm]
ce2_radiomic_clean_glcm  = ce2_radiomic_clean.iloc[:, keep_indices_glcm]
ce4_radiomic_clean_glcm  = ce4_radiomic_clean.iloc[:, keep_indices_glcm]
ce1_radiomic_clean_gldm  = ce1_radiomic_clean.iloc[:, keep_indices_gldm]
ce2_radiomic_clean_gldm  = ce2_radiomic_clean.iloc[:, keep_indices_gldm]
ce4_radiomic_clean_gldm  = ce4_radiomic_clean.iloc[:, keep_indices_gldm]
ce1_radiomic_clean_glrlm = ce1_radiomic_clean.iloc[:, keep_indices_glrlm]
ce2_radiomic_clean_glrlm = ce2_radiomic_clean.iloc[:, keep_indices_glrlm]
ce4_radiomic_clean_glrlm = ce4_radiomic_clean.iloc[:, keep_indices_glrlm]
ce1_radiomic_clean_glszm = ce1_radiomic_clean.iloc[:, keep_indices_glszm]
ce2_radiomic_clean_glszm = ce2_radiomic_clean.iloc[:, keep_indices_glszm]
ce4_radiomic_clean_glszm = ce4_radiomic_clean.iloc[:, keep_indices_glszm]
ce1_radiomic_clean_ngtdm = ce1_radiomic_clean.iloc[:, keep_indices_ngtdm]
ce2_radiomic_clean_ngtdm = ce2_radiomic_clean.iloc[:, keep_indices_ngtdm]
ce4_radiomic_clean_ngtdm = ce4_radiomic_clean.iloc[:, keep_indices_ngtdm]

fl1_radiomic_clean_first = fl1_radiomic_clean.iloc[:, keep_indices_first]
fl2_radiomic_clean_first = fl2_radiomic_clean.iloc[:, keep_indices_first]
fl4_radiomic_clean_first = fl4_radiomic_clean.iloc[:, keep_indices_first]
fl1_radiomic_clean_glcm  = fl1_radiomic_clean.iloc[:, keep_indices_glcm]
fl2_radiomic_clean_glcm  = fl2_radiomic_clean.iloc[:, keep_indices_glcm]
fl4_radiomic_clean_glcm  = fl4_radiomic_clean.iloc[:, keep_indices_glcm]
fl1_radiomic_clean_gldm  = fl1_radiomic_clean.iloc[:, keep_indices_gldm]
fl2_radiomic_clean_gldm  = fl2_radiomic_clean.iloc[:, keep_indices_gldm]
fl4_radiomic_clean_gldm  = fl4_radiomic_clean.iloc[:, keep_indices_gldm]
fl1_radiomic_clean_glrlm = fl1_radiomic_clean.iloc[:, keep_indices_glrlm]
fl2_radiomic_clean_glrlm = fl2_radiomic_clean.iloc[:, keep_indices_glrlm]
fl4_radiomic_clean_glrlm = fl4_radiomic_clean.iloc[:, keep_indices_glrlm]
fl1_radiomic_clean_glszm = fl1_radiomic_clean.iloc[:, keep_indices_glszm]
fl2_radiomic_clean_glszm = fl2_radiomic_clean.iloc[:, keep_indices_glszm]
fl4_radiomic_clean_glszm = fl4_radiomic_clean.iloc[:, keep_indices_glszm]
fl1_radiomic_clean_ngtdm = fl1_radiomic_clean.iloc[:, keep_indices_ngtdm]
fl2_radiomic_clean_ngtdm = fl2_radiomic_clean.iloc[:, keep_indices_ngtdm]
fl4_radiomic_clean_ngtdm = fl4_radiomic_clean.iloc[:, keep_indices_ngtdm]

adc1_radiomic_clean_first = adc1_radiomic_clean.iloc[:, keep_indices_first]
adc2_radiomic_clean_first = adc2_radiomic_clean.iloc[:, keep_indices_first]
adc4_radiomic_clean_first = adc4_radiomic_clean.iloc[:, keep_indices_first]
adc1_radiomic_clean_glcm  = adc1_radiomic_clean.iloc[:, keep_indices_glcm]
adc2_radiomic_clean_glcm  = adc2_radiomic_clean.iloc[:, keep_indices_glcm]
adc4_radiomic_clean_glcm  = adc4_radiomic_clean.iloc[:, keep_indices_glcm]
adc1_radiomic_clean_gldm  = adc1_radiomic_clean.iloc[:, keep_indices_gldm]
adc2_radiomic_clean_gldm  = adc2_radiomic_clean.iloc[:, keep_indices_gldm]
adc4_radiomic_clean_gldm  = adc4_radiomic_clean.iloc[:, keep_indices_gldm]
adc1_radiomic_clean_glrlm = adc1_radiomic_clean.iloc[:, keep_indices_glrlm]
adc2_radiomic_clean_glrlm = adc2_radiomic_clean.iloc[:, keep_indices_glrlm]
adc4_radiomic_clean_glrlm = adc4_radiomic_clean.iloc[:, keep_indices_glrlm]
adc1_radiomic_clean_glszm = adc1_radiomic_clean.iloc[:, keep_indices_glszm]
adc2_radiomic_clean_glszm = adc2_radiomic_clean.iloc[:, keep_indices_glszm]
adc4_radiomic_clean_glszm = adc4_radiomic_clean.iloc[:, keep_indices_glszm]
adc1_radiomic_clean_ngtdm = adc1_radiomic_clean.iloc[:, keep_indices_ngtdm]
adc2_radiomic_clean_ngtdm = adc2_radiomic_clean.iloc[:, keep_indices_ngtdm]
adc4_radiomic_clean_ngtdm = adc4_radiomic_clean.iloc[:, keep_indices_ngtdm]

#%%

ce1_radiomic_clean_final  = ce1_radiomic_clean
ce2_radiomic_clean_final  = ce2_radiomic_clean
ce4_radiomic_clean_final  = ce4_radiomic_clean
adc1_radiomic_clean_final = adc1_radiomic_clean
adc2_radiomic_clean_final = adc2_radiomic_clean
adc4_radiomic_clean_final = adc4_radiomic_clean
fl1_radiomic_clean_final  = fl1_radiomic_clean
fl2_radiomic_clean_final  = fl2_radiomic_clean
fl4_radiomic_clean_final  = fl4_radiomic_clean

#%% Conduct CDKN2A/B Test
from scipy.stats import pointbiserialr
import seaborn as sns
# Predefine the dataframe needed
# OS (Days from diagnosis surgery)

gene_filter_1 = gene_filter_1[gene_filter_1['Assigned ID'].isin(common_indices)]
os_name       = gene_filter_1.columns[29]
status_name   = gene_filter_1.columns[28]

# Single
temp_df = ce1_radiomic_clean_final.add_suffix('_ce0')
temp_df = pd.merge(temp_df, ce2_radiomic_clean_final.add_suffix('_ce1'), left_index=True, right_index=True)
temp_df = pd.merge(temp_df, ce4_radiomic_clean_final.add_suffix('_ce2'), left_index=True, right_index=True)
# temp_df = pd.merge(temp_df, ce0_radiomic_clean, left_index=True, right_index=True, suffixes=('','_ce4'))
# temp_df = pd.merge(temp_df, ce_shape_clean.add_suffix('_ce'), left_index=True, right_index=True)
# temp_df = fl1_radiomic_clean_final.add_suffix('_fl0')
temp_df = pd.merge(temp_df, fl1_radiomic_clean_final.add_suffix('_fl0'), left_index=True, right_index=True)
temp_df = pd.merge(temp_df, fl2_radiomic_clean_final.add_suffix('_fl1'), left_index=True, right_index=True)
temp_df = pd.merge(temp_df, fl4_radiomic_clean_final.add_suffix('_fl2'), left_index=True, right_index=True)
# temp_df = pd.merge(temp_df, fl0_radiomic_clean, left_index=True, right_index=True, suffixes=('', '_fl4'))
# temp_df = pd.merge(temp_df, fl_shape_clean.add_suffix('_fl'), left_index=True, right_index=True)
# temp_df = adc1_radiomic_clean_final.add_suffix('_adc0')
# temp_df = pd.merge(temp_df, adc1_radiomic_clean_final.add_suffix('_adc0'), left_index=True, right_index=True)
# temp_df = pd.merge(temp_df, adc2_radiomic_clean_final.add_suffix('_adc1'), left_index=True, right_index=True)
# temp_df = pd.merge(temp_df, adc4_radiomic_clean_final.add_suffix('_adc2'), left_index=True, right_index=True)
# temp_df = pd.merge(temp_df, adc0_radiomic_clean, left_index=True, right_index=True, suffixes=('', '_adc4'))
# temp_df = pd.merge(temp_df, adc_shape_clean.add_suffix('_adc'), left_index=True, right_index=True)
# 
# temp_df = ce_shape_clean.add_suffix('_ce')
# temp_df = pd.merge(temp_df, fl_shape_clean.add_suffix('_ce'), left_index=True, right_index=True)
# temp_df = pd.merge(temp_df, adc_shape_clean.add_suffix('_adc'), left_index=True, right_index=True)

# Summary
final_df  = temp_df
merged_df = pd.merge(final_df.reset_index(drop=True), gene_filter_1[['MGMT']].reset_index(drop=True), left_index=True, right_index=True)
merged_df = pd.merge(merged_df.reset_index(drop=True), gene_filter_1[['EGFR']].reset_index(drop=True), left_index=True, right_index=True)
merged_df = pd.merge(merged_df, gene_filter_1['IDH'].reset_index(drop=True), left_index=True, right_index=True)
merged_df = pd.merge(merged_df, gene_filter_1[os_name].reset_index(drop=True), left_index=True, right_index=True)
merged_df = pd.merge(merged_df, gene_filter_1[status_name].reset_index(drop=True), left_index=True, right_index=True)
merged_df = merged_df[merged_df['IDH'] == 'negative']

merged_df = pd.merge(merged_df, gene_filter_1['Sex'].reset_index(drop=True), left_index=True, right_index=True)
merged_df = pd.merge(merged_df, gene_filter_1['Age at MRI'].reset_index(drop=True), left_index=True, right_index=True)
# merged_df = pd.merge(merged_df, gene_filter_1['Survival Status (1-dead,0-alive)'].reset_index(drop=True), left_index=True, right_index=True)
merged_df = pd.merge(merged_df, gene_filter_1['PTEN'].reset_index(drop=True), left_index=True, right_index=True)

merged_df = merged_df.dropna(subset=['MGMT'])
merged_df = merged_df.dropna(subset=['EGFR'])
file_numbers_idh = merged_df.index.tolist()

file_numbers_common = [common_indices[i] for i in file_numbers_idh]

#%% Change the gene name here
new_sex = merged_df['Sex']
new_age = merged_df['Age at MRI']

X_data = merged_df.drop(['MGMT','EGFR','IDH','Sex','Age at MRI','PTEN',os_name,status_name], axis=1)
# y    = merged_df[status_name]
y      = merged_df['MGMT']
# y      = merged_df['EGFR']
# y      = merged_df['PTEN']
# y      = merged_df['Survival Status (1-dead,0-alive)']

Labels = y.to_numpy()

# Load external data for external validation

data_path = '/Users/hafeng/Documents/Postdoc_Research_Meetings/TCGA/BraTS2021/'

X_upenn      = np.load(data_path + 'upenn_feature_origin.npy')
# X_upenn      = min_max_scaler.fit_transform(X_upenn)
upenn_MGMT   = np.load(data_path + 'upenn_mgmt_label.npy')
# UCSF_EGFR    = np.load(data_path + 'ucsf_egfr_label.npy')
y_upenn       = upenn_MGMT

#%%
import random
from scipy import stats
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import classification_report, roc_auc_score, confusion_matrix, accuracy_score
from sklearn.model_selection import GridSearchCV
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.model_selection import cross_val_score, StratifiedKFold, KFold
from imblearn.over_sampling import SMOTE
# from sklearn.model_selection import KFold
# from sklearn.metrics import accuracy_score
from sklearn.feature_selection import RFE
from sklearn.preprocessing import StandardScaler, MinMaxScaler
import torch
from torch import nn

X_numeric = X_data.to_numpy()
y_numeric = y.to_numpy()

# Convert your object array to float
X = np.array(X_numeric, dtype=np.float32)
y = np.array(y_numeric).astype(np.float32)

# min_max_scaler = StandardScaler()
min_max_scaler = MinMaxScaler()
X       = min_max_scaler.fit_transform(X)
X_upenn = min_max_scaler.transform(X_upenn)

# X_combined = np.concatenate((X, X_upenn), axis=0)
# X_combined = min_max_scaler.fit_transform(X_combined)

# X = X_combined[:X.shape[0],:]
# X_upenn = X_combined[X.shape[0]:,:]

# Normalize features
# scaler = StandardScaler()
# scaler = MinMaxScaler()
# X = scaler.fit_transform(X)
# X_ucsf = scaler.fit_transform(X_ucsf) * 2

# selector = SelectKBest(score_func=f_classif, k=350)  # e.g., top 100
# X        = selector.fit_transform(X, y)

from sklearn.decomposition import PCA
# pca = PCA(n_components=60)  # or 20–50 depending on your variance explained
# X_scaled = pca.fit_transform(X_scaled)  # use scaled features
from sklearn.feature_selection import SelectKBest, f_classif

selector = SelectKBest(score_func=f_classif, k=100)  # e.g., top 100
# X_combined = selector.fit_transform(np.concatenate((X, X_upenn),axis=0), np.concatenate((y, y_upenn),axis=0))
XX       = selector.fit_transform(X, y)
internal_indices = selector.get_support(indices=True)
 
# print(selector.get_support(indices=True))
XU       = selector.fit_transform(X_upenn, y_upenn)
external_indices = selector.get_support(indices=True)

# X_upenn  = selector.fit_transform(X_upenn)
# print(selector.get_support(indices=True))
# X        = selector.transform(X)

# Combine features
# common_features = np.union1d(internal_indices, external_indices)
# common_features = external_indices
common_features = internal_indices

# --- Apply to original datasets ---
# X       = X[:, common_features]
# X_upenn = X_upenn[:, common_features]

X       = X[:, common_features]
X_upenn = X_upenn[:, common_features]

#%% Add Robust Scaler

from sklearn.preprocessing import RobustScaler
from MLstatkit.stats import Delong_test
# scaler_R = RobustScaler().fit(X)
# X      = scaler_R.transform(X)
# X_ucsf = scaler_R.transform(X_ucsf)

#%% Logistic Regression
from sklearn.model_selection import StratifiedKFold
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
kfold  = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

auc_clf = []
auc_ucsf   = []
best_auc = -1
best_model = None
best_fold = -1
preds_clf = []
true_label = []

# kfold  = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

for fold, (train_idx, test_idx) in enumerate(kfold.split(X, y)):
    print(f"\nFold {fold + 1}")

    X_train, X_test = X[train_idx], X[test_idx]
    y_train, y_test = y[train_idx], y[test_idx]
    
    # selector1 = SelectKBest(score_func=f_classif, k=50)  # e.g., top 100
    # X_train1   = selector1.fit_transform(X_train, y_train)
    # X_test1    = selector1.transform(X_test)
    # X_upenn1   = selector1.transform(X_upenn)

    # Train logistic regression
    # clf = LogisticRegression(penalty='l2', solver='liblinear', max_iter=6000)
    clf = LogisticRegression()
    clf.fit(X_train, y_train)

    # Predict probabilities for class 1
    preds = clf.predict_proba(X_test)[:, 1]

    auc = roc_auc_score(y_test, preds)
    auc_clf.append(auc)
    print(f"AUC for fold {fold + 1}: {auc:.4f}")
    preds_clf.append(preds)
    true_label.append(y_test)
    
    preds  = clf.predict_proba(X_upenn)[:, 1]
    auc1   = roc_auc_score(y_upenn, preds)
    auc_ucsf.append(auc1)
    print(f"UCSF AUC for fold {fold + 1}: {auc1:.4f}")

    if auc > best_auc:
        best_auc = auc
        best_model_clf = clf
        best_fold      = fold + 1

# After loop
mean_auc = np.mean(auc_clf)
std_auc  = np.std(auc_clf, ddof=1)
ci95     = stats.t.interval(0.95, len(auc_clf) - 1, loc=mean_auc, scale=std_auc / np.sqrt(len(auc_clf)))
print(f"\nBest AUC: {best_auc:.4f} from Fold {best_fold}")
print(f"Average AUC: {mean_auc:.4f}")
print(f"95% CI: ({ci95[0]:.4f}, {ci95[1]:.4f})")

mean_auc = np.mean(auc_ucsf)
std_auc  = np.std(auc_ucsf, ddof=1)
ci95     = stats.t.interval(0.95, len(auc_ucsf) - 1, loc=mean_auc, scale=std_auc / np.sqrt(len(auc_ucsf)))
print(f"\nAverage AUC: {mean_auc:.4f}")
print(f"95% CI: ({ci95[0]:.4f}, {ci95[1]:.4f})")

#%% Logistic Regression in Reverse

from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score

auc_clf = []
auc_ucsf   = []
best_auc = -1
best_model = None
best_fold = -1
preds_clf = []
true_label = []

# kfold  = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

for fold, (train_idx, test_idx) in enumerate(kfold.split(X_upenn, y_upenn)):
    print(f"\nFold {fold + 1}")

    X_train, X_test = X_upenn[train_idx], X_upenn[test_idx]
    y_train, y_test = y_upenn[train_idx], y_upenn[test_idx]

    # Train logistic regression
    # clf = LogisticRegression(penalty='l2', solver='liblinear', max_iter=6000)
    clf = LogisticRegression()
    clf.fit(X_train, y_train)

    # Predict probabilities for class 1
    preds = clf.predict_proba(X_test)[:, 1]

    auc = roc_auc_score(y_test, preds)
    auc_clf.append(auc)
    print(f"AUC for fold {fold + 1}: {auc:.4f}")
    preds_clf.append(preds)
    true_label.append(y_test)
    
    preds  = clf.predict_proba(X)[:, 1]
    auc1   = roc_auc_score(y, preds)
    auc_ucsf.append(auc1)
    print(f"UCSF AUC for fold {fold + 1}: {auc1:.4f}")

    if auc > best_auc:
        best_auc = auc
        best_model_clf = clf
        best_fold      = fold + 1

# After loop
mean_auc = np.mean(auc_clf)
std_auc  = np.std(auc_clf, ddof=1)
ci95     = stats.t.interval(0.95, len(auc_clf) - 1, loc=mean_auc, scale=std_auc / np.sqrt(len(auc_clf)))
print(f"\nBest AUC: {best_auc:.4f} from Fold {best_fold}")
print(f"Average AUC: {mean_auc:.4f}")
print(f"95% CI: ({ci95[0]:.4f}, {ci95[1]:.4f})")

mean_auc = np.mean(auc_ucsf)
std_auc  = np.std(auc_ucsf, ddof=1)
ci95     = stats.t.interval(0.95, len(auc_ucsf) - 1, loc=mean_auc, scale=std_auc / np.sqrt(len(auc_ucsf)))
print(f"\nAverage AUC: {mean_auc:.4f}")
print(f"95% CI: ({ci95[0]:.4f}, {ci95[1]:.4f})")
 

#%%
def set_seed(seed=42):
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)  # for multi-GPU
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False
    os.environ['PYTHONHASHSEED'] = str(seed)

set_seed(42)

# Neural network model for binary classification
class BinaryRadiomicsNet(nn.Module):
    def __init__(self, input_dim):
        super(BinaryRadiomicsNet, self).__init__()
        self.model = nn.Sequential(
            nn.Linear(input_dim, 64),
            nn.ReLU(),
            # nn.Dropout(0.1),
            nn.Linear(64, 32),
            nn.ReLU(),
            # nn.Dropout(0.1),
            nn.Linear(32, 1),  # 1 output for binary classification
            nn.Sigmoid()       # Final sigmoid activation
        )

    def forward(self, x):
        return self.model(x)  # logits


device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
kfold  = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

auc_nn = []
auc_ucsf   = []
best_auc = -1
best_model = None
best_fold = -1

preds_nn = []

for fold, (train_idx, test_idx) in enumerate(kfold.split(X, y)):
    print(f"\nFold {fold + 1}")

    X_train, X_test = torch.tensor(X[train_idx],dtype=torch.float32), torch.tensor(X[test_idx],dtype=torch.float32)
    y_train, y_test = torch.tensor(y[train_idx],dtype=torch.float32), torch.tensor(y[test_idx],dtype=torch.float32)

    model = BinaryRadiomicsNet(input_dim=X.shape[1]).to(device)
    criterion = nn.BCELoss()
    # criterion = nn.BCEWithLogitsLoss()
    optimizer = torch.optim.SGD(model.parameters(), lr=0.1, weight_decay=1e-3)
    # optimizer = torch.optim.Adam(model.parameters(), lr=1e-3, weight_decay=1e-3)

    # Training
    model.train()
    for epoch in range(7000):  # Adjust epochs
        optimizer.zero_grad()
        outputs = model(X_train.to(device)).squeeze()
        loss = criterion(outputs, y_train.to(device))
        loss.backward()
        optimizer.step()
        if epoch % 100 == 0 or epoch == 1:
            print(f"Epoch {epoch:3d}, Loss: {loss.item():.4f}")

    # Evaluation
    model.eval()
    with torch.no_grad():
        preds = model(X_test.to(device)).squeeze().cpu().numpy()
        auc = roc_auc_score(y[test_idx], preds)
        auc_nn.append(auc)
        print(f"AUC for fold {fold + 1}: {auc:.4f}")
        
        preds_nn.append(preds)
        
    # Evaluation
    model.eval()
    with torch.no_grad():
        preds = model(torch.tensor(X_upenn, dtype=torch.float32).to(device)).squeeze().cpu().numpy()
        auc = roc_auc_score(y_upenn, preds)
        auc_ucsf.append(auc)
        print(f"UCSF AUC for fold {fold + 1}: {auc:.4f}")

# Compute mean and 95% CI
mean_auc = np.mean(auc_nn)
std_auc  = np.std(auc_nn, ddof=1)
ci95     = stats.t.interval(0.95, len(auc_nn) - 1, loc=mean_auc, scale=std_auc / np.sqrt(len(auc_nn)))
print(f"\nAverage AUC: {mean_auc:.4f}")
print(f"95% CI: ({ci95[0]:.4f}, {ci95[1]:.4f})")

mean_auc = np.mean(auc_ucsf)
std_auc  = np.std(auc_ucsf, ddof=1)
ci95     = stats.t.interval(0.95, len(auc_ucsf) - 1, loc=mean_auc, scale=std_auc / np.sqrt(len(auc_ucsf)))
print(f"\nAverage AUC: {mean_auc:.4f}")
print(f"95% CI: ({ci95[0]:.4f}, {ci95[1]:.4f})")

#%% Random Forest

from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score
from scipy import stats

auc_rf = []
auc_ucsf   = []
best_auc = -1
best_model_clf = None
best_fold = -1
preds_rf  = []
for fold, (train_idx, test_idx) in enumerate(kfold.split(X, y)):
    print(f"\nFold {fold + 1}")

    X_train, X_test = X[train_idx], X[test_idx]
    y_train, y_test = y[train_idx], y[test_idx]

    # Train Random Forest
    rf = RandomForestClassifier(n_estimators=100, max_depth=None, random_state=42)
    rf.fit(X_train, y_train)

    preds = rf.predict_proba(X_test)[:, 1]
    auc   = roc_auc_score(y_test, preds)
    auc_rf.append(auc)
    print(f"AUC for fold {fold + 1}: {auc:.4f}")
    preds_rf.append(preds)
    
    preds = rf.predict_proba(X_upenn)[:, 1]
    auc1   = roc_auc_score(y_upenn, preds)
    auc_ucsf.append(auc1)
    print(f"AUC for fold {fold + 1}: {auc1:.4f}")

    if auc > best_auc:
        best_auc = auc
        best_model_rf = rf
        best_fold     = fold + 1

mean_auc = np.mean(auc_rf)
std_auc  = np.std(auc_rf, ddof=1)
ci95     = stats.t.interval(0.95, len(auc_rf) - 1, loc=mean_auc, scale=std_auc / np.sqrt(len(auc_rf)))
print(f"\nBest AUC: {best_auc:.4f} from Fold {best_fold}")
print(f"Average AUC: {mean_auc:.4f}")
print(f"95% CI: ({ci95[0]:.4f}, {ci95[1]:.4f})")

mean_auc = np.mean(auc_ucsf)
std_auc  = np.std(auc_ucsf, ddof=1)
ci95     = stats.t.interval(0.95, len(auc_ucsf) - 1, loc=mean_auc, scale=std_auc / np.sqrt(len(auc_ucsf)))
print(f"\nAverage AUC: {mean_auc:.4f}")
print(f"95% CI: ({ci95[0]:.4f}, {ci95[1]:.4f})")

#%% TPOT

from tpot import TPOTClassifier
from sklearn.metrics import roc_auc_score
from scipy import stats

auc_scores = []
auc_ucsf   = []
best_auc = -1
best_model_tpot = None
best_fold = -1
preds_tpot = []

for fold, (train_idx, test_idx) in enumerate(kfold.split(X, y)):
    print(f"\nFold {fold + 1}")

    X_train, X_test = X[train_idx], X[test_idx]
    y_train, y_test = y[train_idx], y[test_idx]

    # TPOT classifier (generations and population size can be tuned)
    tpot = TPOTClassifier(generations=5, population_size=20, verbosity=2, scoring='roc_auc', random_state=42, n_jobs=-1)
    tpot.fit(X_train, y_train)

    preds = tpot.predict(X_test)
    auc = roc_auc_score(y_test, preds)
    auc_scores.append(auc)
    print(f"AUC for fold {fold + 1}: {auc:.4f}")
    preds_tpot.append(preds)
    
    preds = tpot.predict(X_upenn)
    auc = roc_auc_score(y_upenn, preds)
    auc_ucsf.append(auc)
    print(f"AUC for fold {fold + 1}: {auc:.4f}")

    if auc > best_auc:
        best_auc = auc
        best_model_tpot = tpot
        best_fold = fold + 1

mean_auc = np.mean(auc_scores)
std_auc = np.std(auc_scores, ddof=1)
ci95 = stats.t.interval(0.95, len(auc_scores) - 1, loc=mean_auc, scale=std_auc / np.sqrt(len(auc_scores)))

print(f"\nBest AUC: {best_auc:.4f} from Fold {best_fold}")
print(f"Average AUC: {mean_auc:.4f}")
print(f"95% CI: ({ci95[0]:.4f}, {ci95[1]:.4f})")

mean_auc = np.mean(auc_ucsf)
std_auc  = np.std(auc_ucsf, ddof=1)
ci95     = stats.t.interval(0.95, len(auc_ucsf) - 1, loc=mean_auc, scale=std_auc / np.sqrt(len(auc_ucsf)))
print(f"\nAverage AUC: {mean_auc:.4f}")
print(f"95% CI: ({ci95[0]:.4f}, {ci95[1]:.4f})")

#%%



#%%


#%%


#%% Perform SHAP Analysis
import shap
selected_feature_names = X_data.columns[selector.get_support()]

explainer = shap.DeepExplainer(model, torch.tensor(X_train))  # If model is Keras
shap_values = explainer(torch.tensor(X_test))
shap.summary_plot(shap_values, X_test, feature_names=selected_feature_names)

# selected_feature_names = X_data.columns[selector.get_support()]

explainer = shap.Explainer(clf, X_train)  # If model is Keras
shap_values = explainer(X_test)
# shap.summary_plot(shap_values, X_test, feature_names=selected_feature_names)

# Set figure size before plotting
plt.figure(figsize=(20, 10))  # Larger figure
shap.summary_plot(shap_values, X_test, feature_names=selected_feature_names)
plt.tight_layout(pad=2.0)  # Automatically adjust layout
plt.show()

# Extract SHAP values
sv = shap_values.values  # shape: (n_samples, n_features)

# Compute mean absolute SHAP per feature
shap_mean_abs = np.abs(sv).mean(axis=0)  # shape: (n_features,)

# Get sorted indices (largest importance first)
sorted_indices = np.argsort(shap_mean_abs)[::-1]

# Top N features
top_n = 20
top_indices = sorted_indices[:top_n]

print("Top feature indices:", top_indices)
print("Top feature names:", [selected_feature_names[i] for i in top_indices])

#%%

import shap
import torch
import numpy as np

# Define prediction function
def predict_fn(X_numpy):
    model.eval()
    X_tensor = torch.tensor(X_numpy, dtype=torch.float32)
    with torch.no_grad():
        logits = model(X_tensor)
        probs = torch.softmax(logits, dim=1).numpy()
    return probs

# Use a sample as the background distribution
X_background = shap.sample(X_train, 100)  # smaller sample for speed

# Create SHAP explainer
explainer = shap.DeepExplainer(model, torch.tensor(X_background))
shap_values = explainer.shap_values(torch.tensor(X_test))

# Plot
shap.summary_plot(np.squeeze(shap_values), X_test, feature_names=selected_feature_names, max_display=20)

#%% Count name numbers

count_shape = sum(1 for name in selected_feature_names if "shape" in name.lower())
count_first = sum(1 for name in selected_feature_names if "firstorder" in name.lower())
count_glcm  = sum(1 for name in selected_feature_names if "glcm" in name.lower())
count_gldm  = sum(1 for name in selected_feature_names if "gldm" in name.lower())
count_glrlm = sum(1 for name in selected_feature_names if "glrlm" in name.lower())
count_glszm = sum(1 for name in selected_feature_names if "glszm" in name.lower())
count_ngtdm = sum(1 for name in selected_feature_names if "ngtdm" in name.lower())

print(count_shape / len(selected_feature_names))
print(count_first / len(selected_feature_names))
print(count_glcm  / len(selected_feature_names))
print(count_gldm  / len(selected_feature_names))
print(count_glrlm / len(selected_feature_names))
print(count_glszm / len(selected_feature_names))
print(count_ngtdm / len(selected_feature_names))

#%% Generate Box Plot

origin_percent = [0.70,19.2,25.6,15.0,17.1,17.1,5.3]
mgmt_percent   = [9, 11.4, 33.4, 10.2, 16.8, 16.8, 2.4]
egfr_percent   = [0, 21.2, 37.2, 12.0, 9.0, 9.0, 11.6]

# X locations
labels = ['Shape','First Order','GLCM', 'GLDM', 'GLRLM', 'GLSZM', 'NGTDM']  # you can change group names
x = np.arange(len(labels))  # the label locations
width = 0.25  # width of each bar

# Create plot
fig, ax = plt.subplots(figsize=(10, 6))
rects1 = ax.bar(x - width, origin_percent, width, label='All Features')
rects2 = ax.bar(x, mgmt_percent, width, label='Top 500 features from MGMT')
rects3 = ax.bar(x + width, egfr_percent, width, label='Top 500 features from EGFR')

# Add labels, title, and legend
ax.set_ylabel('Percent (%)', fontsize=16)
# ax.set_xlabel('Groups', fontsize=16)
# ax.set_title('Comparison of Percentages')
ax.set_xticks(x)
ax.set_xticklabels(labels, fontsize=16)
ax.tick_params(axis='y', labelsize=16)
ax.legend(fontsize=16)

# Show plot
plt.tight_layout()
plt.show()

#%% Perform a DeLong Test

true = true_label
delong_pval = []

prob_A = preds_nn
prob_B = preds_clf

for i in range(len(true)):
    # Perform DeLong's test
    z_score, p_value = Delong_test(np.array(true[i]), np.array(prob_A[i]), np.array(prob_B[i]))
    delong_pval.append(p_value)

mean_pval = np.mean(delong_pval)
std_pval  = np.std(delong_pval, ddof=1)
ci95      = stats.t.interval(0.95, len(delong_pval) - 1, loc=mean_pval, scale=std_pval / np.sqrt(len(delong_pval)))
print(f"\nAverage P-Value: {mean_pval:.4f}")
print(f"95% CI: ({ci95[0]:.4f}, {ci95[1]:.4f})")


#%% Generate the ROC Curve
from sklearn.metrics import roc_curve, auc
import matplotlib.pyplot as plt

for i in range(len(true_label)):
    
    fpr_nn, tpr_nn, _ = roc_curve(true_label[i], preds_nn[i])
    roc_auc_nn = auc(fpr_nn, tpr_nn)
    
    fpr_clf, tpr_clf, _ = roc_curve(true_label[i], preds_clf[i])
    roc_auc_clf = auc(fpr_clf, tpr_clf)
    
    fpr_rf, tpr_rf, _ = roc_curve(true_label[i], preds_rf[i])
    roc_auc_rf = auc(fpr_rf, tpr_rf)
    
    fpr_tpot, tpr_tpot, _ = roc_curve(true_label[i], preds_tpot[i])
    roc_auc_tpot = auc(fpr_tpot, tpr_tpot)
    
    
    plt.figure()
    plt.plot(fpr_nn, tpr_nn, color='red', lw=1.5, 
             label=f'NN (AUC = {roc_auc_nn:.2f})')
    plt.plot(fpr_clf, tpr_clf, color='darkorange', lw=1.5, 
             label=f'LR (AUC = {roc_auc_clf:.2f})')
    plt.plot(fpr_rf, tpr_rf, color='darkgreen', lw=1.5, 
             label=f'RF (AUC = {roc_auc_rf:.2f})')
    plt.plot(fpr_tpot, tpr_tpot, color='blue', lw=1.5, 
             label=f'TPOT (AUC = {roc_auc_tpot:.2f})')
    
    plt.plot([0, 1], [0, 1], color='navy', lw=1.5, linestyle='--')  # diagonal line
    
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title(f'Receiver Operating Characteristic (ROC) on folds {i+1}')
    plt.legend(loc='lower right')
    plt.show()

#%% Plot the Heatmap for the radiomics

os_name       = gene_filter_1.columns[29]
status_name   = gene_filter_1.columns[28]

# Single
temp_df1 = ce1_radiomic_clean.add_suffix('_ce0')
temp_df1 = pd.merge(temp_df1, ce4_radiomic_clean.add_suffix('_ce1'), left_index=True, right_index=True)
temp_df1 = pd.merge(temp_df1, ce2_radiomic_clean.add_suffix('_ce2'), left_index=True, right_index=True)
# temp_df = pd.merge(temp_df, ce0_radiomic_clean, left_index=True, right_index=True, suffixes=('','_ce0'))
temp_df1 = pd.merge(temp_df1, ce_shape_clean.add_suffix('_ce'), left_index=True, right_index=True)
# temp_df = fl1_radiomic_clean_add
temp_df1 = pd.merge(temp_df1, fl1_radiomic_clean.add_suffix('_fl0'), left_index=True, right_index=True)
temp_df1 = pd.merge(temp_df1, fl4_radiomic_clean.add_suffix('_fl1'), left_index=True, right_index=True)
temp_df1 = pd.merge(temp_df1, fl2_radiomic_clean.add_suffix('_fl2'), left_index=True, right_index=True)
# temp_df = pd.merge(temp_df, fl0_radiomic_clean, left_index=True, right_index=True, suffixes=('', '_fl0'))
temp_df1 = pd.merge(temp_df1, fl_shape_clean.add_suffix('_fl'), left_index=True, right_index=True)
# temp_df = adc1_radiomic_clean_add
temp_df1 = pd.merge(temp_df1, adc1_radiomic_clean.add_suffix('_adc0'), left_index=True, right_index=True)
temp_df1 = pd.merge(temp_df1, adc4_radiomic_clean.add_suffix('_adc1'), left_index=True, right_index=True)
temp_df1 = pd.merge(temp_df1, adc2_radiomic_clean.add_suffix('_adc2'), left_index=True, right_index=True)
# temp_df = pd.merge(temp_df, adc0_radiomic_clean, left_index=True, right_index=True, suffixes=('', '_adc0'))
temp_df1 = pd.merge(temp_df1, adc_shape_clean.add_suffix('_adc'), left_index=True, right_index=True)

# Summary
final_df1  = temp_df1
merged_df1 = pd.merge(final_df1.reset_index(drop=True), gene_filter_1[['MGMT']].reset_index(drop=True), left_index=True, right_index=True)
merged_df1 = pd.merge(merged_df1.reset_index(drop=True), gene_filter_1[['EGFR']].reset_index(drop=True), left_index=True, right_index=True)
merged_df1 = pd.merge(merged_df1.reset_index(drop=True), gene_filter_1[['PTEN']].reset_index(drop=True), left_index=True, right_index=True)
merged_df1 = pd.merge(merged_df1.reset_index(drop=True), gene_filter_1[['TP53']].reset_index(drop=True), left_index=True, right_index=True)

merged_df1 = pd.merge(merged_df1, gene_filter_1['IDH'].reset_index(drop=True), left_index=True, right_index=True)
merged_df1 = pd.merge(merged_df1, gene_filter_1[os_name].reset_index(drop=True), left_index=True, right_index=True)
merged_df1 = pd.merge(merged_df1, gene_filter_1[status_name].reset_index(drop=True), left_index=True, right_index=True)
merged_df1 = merged_df1[merged_df1['IDH'] == 'negative']

merged_df1 = merged_df1.dropna(subset=['MGMT'])
merged_df1 = merged_df1.dropna(subset=['EGFR'])
merged_df1 = merged_df1.dropna(subset=['TP53'])
merged_df1 = merged_df1.dropna(subset=['PTEN'])

shap_names = [selected_feature_names[i] for i in top_indices]

# shap_names = [c for c in shap_names if c in merged_df1.columns]
feature_hm = merged_df1[shap_names]
print(feature_hm.shape)

feature_hm = merged_df1[shap_names]
mgmt_hm    = merged_df1['MGMT']
egfr_hm    = merged_df1['EGFR']
tp53_hm    = merged_df1['TP53']
pten_hm    = merged_df1['PTEN']

print("merged_df1.shape:", merged_df1.shape)
print("Number of SHAP features found:", len(shap_names))
print("feature_hm.shape:", feature_hm.shape)

# df_hm = feature_hm
# df_hm = pd.merge(df_hm, mgmt_hm, left_index=True, right_index=True) 
# df_hm = pd.merge(df_hm, egfr_hm, left_index=True, right_index=True) 
# df_hm = pd.merge(df_hm, tp53_hm, left_index=True, right_index=True) 
# df_hm = pd.merge(df_hm, pten_hm, left_index=True, right_index=True) 

#%%

gene_cols = ['MGMT','EGFR','TP53','PTEN']
feature_cols = shap_names

# Compute correlations
corr_matrix = pd.DataFrame(index=feature_cols, columns=gene_cols)

for g in gene_cols:
    for f in feature_cols:
        corr_matrix.loc[f, g] = merged_df1[f].corr(merged_df1[g])

# Convert to numeric (corr returns float or NaN)
corr_matrix = corr_matrix.astype(float)

# Plot heatmap
plt.figure(figsize=(10, 8))
sns.heatmap(corr_matrix, annot=False, cmap='coolwarm', center=0)
plt.title('Feature–Gene Correlation Heatmap')
plt.tight_layout()
plt.show()


# Hierarchical clustering heatmap
sns.clustermap(
    corr_matrix,
    cmap='coolwarm',
    center=0,
    method='average',   # linkage method: 'average', 'ward', 'single', etc.
    metric='correlation', # distance metric
    figsize=(10, 12)
)

plt.show()


#%% Generate the heatmap for AUC

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Example: AUC + 95% CI
data_auc = {
    'MGMT': [0.51, 0.74, 0.77, 0.77, 0.65, 0.67, 0.69],
    'EGFR': [0.52, 0.67, 0.68, 0.60, 0.65, 0.65, 0.58]
}
data_ci_lower = {
    'MGMT': [0.45, 0.69, 0.71, 0.71, 0.59, 0.62, 0.65],
    'EGFR': [0.42, 0.58, 0.58, 0.49, 0.55, 0.54, 0.52]
}
data_ci_upper = {
    'MGMT': [0.57, 0.8, 0.83, 0.84, 0.71, 0.73, 0.74],
    'EGFR': [0.63, 0.76, 0.79, 0.71, 0.75, 0.76, 0.64]
}

features = ['Shape', 'First Order', 'GLCM', 'GLDM', 'GLRLM', 'GLSZM', 'NGTDM']

auc_df = pd.DataFrame(data_auc, index=features)
ci_lower_df = pd.DataFrame(data_ci_lower, index=features)
ci_upper_df = pd.DataFrame(data_ci_upper, index=features)

# Create annotation text with CI
annot_df = auc_df.copy()
for gene in auc_df.columns:
    annot_df[gene] = auc_df[gene].round(2).astype(str) + "\n(" + \
                     ci_lower_df[gene].round(2).astype(str) + "-" + \
                     ci_upper_df[gene].round(2).astype(str) + ")"

# Plot heatmap
plt.figure(figsize=(8, 6))
sns.heatmap(
    auc_df,
    annot=annot_df,
    fmt="",
    cmap='coolwarm',
    vmin=0.5, vmax=1.0
)
plt.title('AUC Heatmap with 95% CI')
plt.ylabel('Radiomic Feature')
plt.xlabel('Gene')
plt.show()

#%%

#%% Radiomic Regional Analysis

# Single
temp_df = ce1_radiomic_clean_final.add_suffix('_ce0')
temp_df = pd.merge(temp_df, ce4_radiomic_clean_final.add_suffix('_ce1'), left_index=True, right_index=True)
temp_df = pd.merge(temp_df, ce2_radiomic_clean_final.add_suffix('_ce2'), left_index=True, right_index=True)
temp_df = pd.merge(temp_df, ce0_radiomic_clean, left_index=True, right_index=True, suffixes=('','_ce4'))
# temp_df = pd.merge(temp_df, ce_shape_clean.add_suffix('_ce'), left_index=True, right_index=True)
# temp_df = fl1_radiomic_clean.add_suffix('_fl1')
temp_df = pd.merge(temp_df, fl1_radiomic_clean_final.add_suffix('_fl0'), left_index=True, right_index=True)
temp_df = pd.merge(temp_df, fl4_radiomic_clean_final.add_suffix('_fl1'), left_index=True, right_index=True)
temp_df = pd.merge(temp_df, fl2_radiomic_clean_final.add_suffix('_fl2'), left_index=True, right_index=True)
temp_df = pd.merge(temp_df, fl0_radiomic_clean, left_index=True, right_index=True, suffixes=('', '_fl4'))
# temp_df = pd.merge(temp_df, fl_shape_clean.add_suffix('_fl'), left_index=True, right_index=True)
# temp_df = adc1_radiomic_clean_add
temp_df = pd.merge(temp_df, adc1_radiomic_clean_final.add_suffix('_adc0'), left_index=True, right_index=True)
temp_df = pd.merge(temp_df, adc4_radiomic_clean_final.add_suffix('_adc1'), left_index=True, right_index=True)
temp_df = pd.merge(temp_df, adc2_radiomic_clean_final.add_suffix('_adc2'), left_index=True, right_index=True)
temp_df = pd.merge(temp_df, adc0_radiomic_clean, left_index=True, right_index=True, suffixes=('', '_adc4'))
# temp_df = pd.merge(temp_df, adc_shape_clean.add_suffix('_adc'), left_index=True, right_index=True)

# temp_df = ce_shape_clean.add_suffix('_ce')
# temp_df = pd.merge(temp_df, fl_shape_clean.add_suffix('_ce'), left_index=True, right_index=True)
# temp_df = pd.merge(temp_df, adc_shape_clean.add_suffix('_adc'), left_index=True, right_index=True)

# Summary
final_df  = temp_df
merged_df = pd.merge(final_df.reset_index(drop=True), gene_filter_1[['MGMT']].reset_index(drop=True), left_index=True, right_index=True)
merged_df = pd.merge(merged_df.reset_index(drop=True), gene_filter_1[['EGFR']].reset_index(drop=True), left_index=True, right_index=True)
merged_df = pd.merge(merged_df, gene_filter_1['IDH'].reset_index(drop=True), left_index=True, right_index=True)
merged_df = pd.merge(merged_df, gene_filter_1[os_name].reset_index(drop=True), left_index=True, right_index=True)
merged_df = pd.merge(merged_df, gene_filter_1[status_name].reset_index(drop=True), left_index=True, right_index=True)
merged_df = merged_df[merged_df['IDH'] == 'negative']

merged_df = pd.merge(merged_df, gene_filter_1['Sex'].reset_index(drop=True), left_index=True, right_index=True)
merged_df = pd.merge(merged_df, gene_filter_1['Age at MRI'].reset_index(drop=True), left_index=True, right_index=True)
# merged_df = pd.merge(merged_df, gene_filter_1['Survival Status (1-dead,0-alive)'].reset_index(drop=True), left_index=True, right_index=True)
merged_df = pd.merge(merged_df, gene_filter_1['PTEN'].reset_index(drop=True), left_index=True, right_index=True)

merged_df = merged_df.dropna(subset=['MGMT'])
merged_df = merged_df.dropna(subset=['EGFR'])
file_numbers_idh = merged_df.index.tolist()

X_data = merged_df.drop(['MGMT','EGFR','IDH','Sex','Age at MRI','PTEN',os_name,status_name], axis=1)
# y_data    = merged_df[status_name]
# y_data      = merged_df['MGMT']
# y_data      = merged_df['EGFR']
# y_data      = merged_df['PTEN']
# y_data      = merged_df['Survival Status (1-dead,0-alive)']

from sklearn.preprocessing import StandardScaler, MinMaxScaler
X_data = MinMaxScaler().fit_transform(X_data)

ce0_feature = X_data[:,:1860*1]
ce1_feature = X_data[:,1860*1:1860*2]
ce2_feature = X_data[:,1860*2:1860*3]
ce4_feature = X_data[:,1860*3:1860*4]
fl0_feature = X_data[:,1860*4:1860*5]
fl1_feature = X_data[:,1860*5:1860*6]
fl2_feature = X_data[:,1860*6:1860*7]
fl4_feature = X_data[:,1860*7:1860*8]
adc0_feature = X_data[:,1860*8:1860*9]
adc1_feature = X_data[:,1860*9:1860*10]
adc2_feature = X_data[:,1860*10:1860*11]
adc4_feature = X_data[:,1860*11:]

#%% Generate the analysis plot

# First Order Feature: 0 - 17
# GLCM: 18 - 41
# GLDM: 42 - 55
# GLRLM: 56 - 71
# GLSZM: 72 - 87
# NGTDM: 88 - 92

from cliffs_delta import cliffs_delta  # pip install cliffs-delta

p_thres  = 0.05
p_thres0 = 0.05

sig01_list = []
sig12_list = []
sig24_list = []

effect_size_01 = []
effect_size_12 = []
effect_size_24 = []

delta_size_01 = []
delta_size_12 = []
delta_size_24 = []

def cohens_dz(x, y):
    diff = np.array(x) - np.array(y)
    return np.mean(diff) / np.std(diff, ddof=1)

for j in range(1):
    radiomic_index = j+8
    print(f'Working on Radiomic {j}')
    
    list_index = []
    
    for i in range(20):
        list_index.append(radiomic_index + i*93)
    
    ce0_filter = fl0_feature[:,list_index]
    ce1_filter = fl1_feature[:,list_index]
    ce2_filter = fl2_feature[:,list_index]
    ce4_filter = fl4_feature[:,list_index]
    
    sig01 = 0
    sig12 = 0
    sig24 = 0
    
    cohen_list01 = []
    cohen_list12 = []
    cohen_list24 = []
    
    delta_list01 = []
    delta_list12 = []
    delta_list24 = []
    
    for patient_id in range(ce0_filter.shape[0]):
        # patient_id = 0
        data = [ce0_filter[patient_id,:],ce1_filter[patient_id,:],ce2_filter[patient_id,:],ce4_filter[patient_id,:]]
        
        from scipy.stats import shapiro
        from scipy.stats import mannwhitneyu
        from scipy.stats import ttest_ind
        
        stat, p1 = shapiro(ce0_filter[patient_id,:])
        # print("Necrotic p =", p1)  # p > 0.05 means data looks normal
        stat, p2 = shapiro(ce1_filter[patient_id,:])
        # print("T1 p =", p2)  # p > 0.05 means data looks normal
        stat, p3 = shapiro(ce2_filter[patient_id,:])
        # print("T2 p =", p3)  # p > 0.05 means data looks normal
        stat, p4 = shapiro(ce4_filter[patient_id,:])
        # print("2cm p =", p4)  # p > 0.05 means data looks normal
        
        cohen_list01.append(cohens_dz(ce0_filter[patient_id,:], ce1_filter[patient_id,:]))
        cohen_list12.append(cohens_dz(ce1_filter[patient_id,:], ce2_filter[patient_id,:]))
        cohen_list24.append(cohens_dz(ce2_filter[patient_id,:], ce4_filter[patient_id,:]))
        
        delta_list01.append(np.abs(cliffs_delta(ce0_filter[patient_id,:], ce1_filter[patient_id,:])[0]))
        delta_list12.append(np.abs(cliffs_delta(ce1_filter[patient_id,:], ce2_filter[patient_id,:])[0]))
        delta_list24.append(np.abs(cliffs_delta(ce2_filter[patient_id,:], ce4_filter[patient_id,:])[0]))
        
        if p1>p_thres0 and p2>p_thres0 and p3>p_thres0 and p4>p_thres0:
            stat01, p01 = ttest_ind(ce0_filter[patient_id,:], ce1_filter[patient_id,:])
            stat12, p12 = ttest_ind(ce1_filter[patient_id,:], ce2_filter[patient_id,:])
            stat24, p24 = ttest_ind(ce2_filter[patient_id,:], ce4_filter[patient_id,:])
            # print(f"Patient {patient_id} has T-test p01={p01}, p12={p12}, p24={p24}")
            
            if p01 < p_thres:
                sig01 += 1
            if p12 < p_thres:
                sig12 += 1
            if p24 < p_thres:
                sig24 += 1
             
        else:
            stat01, p01 = mannwhitneyu(ce0_filter[patient_id,:], ce1_filter[patient_id,:], alternative='two-sided')
            stat12, p12 = mannwhitneyu(ce1_filter[patient_id,:], ce2_filter[patient_id,:], alternative='two-sided')
            stat24, p24 = mannwhitneyu(ce2_filter[patient_id,:], ce4_filter[patient_id,:], alternative='two-sided')
            # print(f"Patient {patient_id} has Mann-Whitney p01={p01}, p12={p12}, p24={p24}")
            
            if p01 < p_thres:
                sig01 += 1
            if p12 < p_thres:
                sig12 += 1
            if p24 < p_thres:
                sig24 += 1
                
    sig01_list.append(sig01)
    sig12_list.append(sig12)
    sig24_list.append(sig24)
    
    effect_size_01.append(np.mean(cohen_list01))
    effect_size_12.append(np.mean(cohen_list12))
    effect_size_24.append(np.mean(cohen_list24))
    
    delta_size_01.append(np.mean(delta_list01))
    delta_size_12.append(np.mean(delta_list12))
    delta_size_24.append(np.mean(delta_list24))

#%%
plt.figure()
plt.boxplot(data, labels=['Necrotic', 'T1', 'T2', '2cm'])

# plt.title("Comparison of Methods", fontsize=12)
plt.ylabel("Values")
plt.show()

#%%

indices01 = [i for i, val in enumerate(sig01_list) if val > 0]
print(indices01)
indices12 = [i for i, val in enumerate(sig12_list) if val > 0]
print(indices12)
indices24 = [i for i, val in enumerate(sig24_list) if val > 0]
print(indices24)

#%% t-SNE for individual feature

# First Order Feature: 0 - 17
# GLCM: 18 - 41
# GLDM: 42 - 55
# GLRLM: 56 - 71
# GLSZM: 72 - 87
# NGTDM: 88 - 92

sne_score = []

for j in range(93):
    radiomic_index = j+0
    print(f'Working on Radiomic {ce_list[radiomic_index]}')
    
    list_index = []
    
    for i in range(20):
        list_index.append(radiomic_index + i*93)
    
    ce0_filter = ce0_feature[:,list_index]
    ce1_filter = ce1_feature[:,list_index]
    ce2_filter = ce2_feature[:,list_index]
    ce4_filter = ce4_feature[:,list_index]
    
    X = np.array(ce0_filter, dtype=np.float64)
    y = np.array(merged_df['MGMT'], dtype=np.float64)  # binary labels

    # Run t-SNE
    tsne = TSNE(n_components=2, random_state=42, perplexity=30)
    embedding = tsne.fit_transform(X)
    
    from sklearn.metrics import silhouette_score

    score = silhouette_score(embedding, y)
    sne_score.append(score)
    
    # print("Silhouette score:", score)
    if score > 0.25:
        print(f'There is strong separation at {ce_list[radiomic_index]} with value {score}')

        # Plot
    
        plt.figure(figsize=(8,6))
        for lab, color in zip([0,1], ["steelblue", "darkorange"]):
            idx = y == lab
            plt.scatter(embedding[idx,0], embedding[idx,1], c=color, label=f"Class {lab}", alpha=0.7, edgecolors="k", s=50)
    
        plt.xlabel("t-SNE Dimension 1")
        plt.ylabel("t-SNE Dimension 2")
        plt.title("t-SNE projection of Radiomic Features")
        plt.legend()
        plt.tight_layout()
        plt.show()

#%% 2D t-SNE for multiple features
# radiomic_index = 8
radiomic_index = np.linspace(0,17,18)
# radiomic_index = np.linspace(18,41,24)
# radiomic_index = np.linspace(42,55,14)
# radiomic_index = np.linspace(56,71,16)
# radiomic_index = np.linspace(72,87,16)
# radiomic_index = np.linspace(88,92,5)

# print(f'Working on Radiomic {ce_list[radiomic_index]}')

list_index = []

for i in range(20):
    list_index.append(radiomic_index + i*93)
    
list_index = np.concatenate(list_index).astype(int).tolist()

ce0_filter = ce0_feature[:,list_index]
ce1_filter = ce1_feature[:,list_index]
ce2_filter = ce2_feature[:,list_index]
ce4_filter = ce4_feature[:,list_index]

# X = np.array(ce4_filter, dtype=np.float64)
X = np.concatenate((ce0_filter, ce1_filter, ce2_filter, ce4_filter),axis=1)
y = np.array(merged_df['MGMT'], dtype=np.float64)  # binary labels

# Run t-SNE
tsne = TSNE(n_components=2, random_state=42, perplexity=15)
embedding = tsne.fit_transform(X)

from sklearn.metrics import silhouette_score

score = silhouette_score(embedding, y)
print("Silhouette score:", score)

# Plot
plt.figure(figsize=(8,6))
for lab, color in zip([0,1], ["steelblue", "darkorange"]):
    idx = y == lab
    plt.scatter(embedding[idx,0], embedding[idx,1], c=color, label=f"Class {lab}", alpha=0.7, edgecolors="k", s=50)

plt.xlabel("t-SNE Dimension 1")
plt.ylabel("t-SNE Dimension 2")
plt.title("t-SNE projection of Radiomic Features")
plt.legend()
plt.tight_layout()
plt.show()

from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.metrics import silhouette_score
import matplotlib.pyplot as plt

lda = LinearDiscriminantAnalysis(n_components=1)   # 1D projection for binary labels
proj = lda.fit_transform(X, y)                     # shape (n_samples, 1)

# optional 2D by combining LDA + PCA
from sklearn.decomposition import PCA
pca = PCA(n_components=1).fit_transform(X)
emb2 = np.hstack([proj, pca])   # 2D: [LDA, PCA]

print("LDA silhouette:", silhouette_score(emb2, y))
plt.figure()
plt.scatter(emb2[:,0], emb2[:,1], c=y, cmap="coolwarm", alpha=0.7)
# plt.title("LDA+PCA projection")
plt.xlabel('LDA axis', fontsize=12)
plt.ylabel('PCA axis', fontsize=12)
plt.show()

import umap
reducer = umap.UMAP(n_components=2, random_state=42, target_metric='categorical')
emb_super = reducer.fit_transform(X, y)   # supply labels
print("supervised UMAP silhouette:", silhouette_score(emb_super, y))
emb_unsup = reducer.fit_transform(X)
print("unsupervised UMAP silhouette:", silhouette_score(emb_unsup, y))

#%% 3D t-SNE for multiple features
# radiomic_index = 8
# radiomic_index = np.linspace(0,17,18)
# radiomic_index = np.linspace(18,41,24)
radiomic_index = np.linspace(42,55,14)
# radiomic_index = np.linspace(56,71,16)
# radiomic_index = np.linspace(72,87,16)
# radiomic_index = np.linspace(88,92,5)
# print(f'Working on Radiomic {ce_list[radiomic_index]}')

list_index = []

for i in range(20):
    list_index.append(radiomic_index + i*93)
    
list_index = np.concatenate(list_index).astype(int).tolist()

ce0_filter = ce0_feature[:,list_index]
ce1_filter = ce1_feature[:,list_index]
ce2_filter = ce2_feature[:,list_index]
ce4_filter = ce4_feature[:,list_index]

# X = np.array(ce4_filter, dtype=np.float64)
# X = np.concatenate((ce0_filter, ce1_filter, ce2_filter, ce4_filter),axis=1)
X = np.concatenate((ce1_filter, ce2_filter, ce4_filter),axis=1)
y = np.array(merged_df['MGMT'], dtype=np.float64)  # binary labels

# Run t-SNEß
tsne = TSNE(n_components=3, random_state=42, perplexity=15)
embedding = tsne.fit_transform(X)

from sklearn.metrics import silhouette_score

score = silhouette_score(embedding, y)
print("Silhouette score:", score)

# Plot

fig = plt.figure(figsize=(8,6))
ax  = fig.add_subplot(projection='3d')
for lab, color in zip([0,1], ["steelblue", "darkorange"]):
    idx = y == lab
    ax.scatter(embedding[idx,0], embedding[idx,1], embedding[idx,2], c=color, label=f"Class {lab}", alpha=0.7, edgecolors="k", s=50)

plt.xlabel("t-SNE Dimension 1")
plt.ylabel("t-SNE Dimension 2")
plt.title("t-SNE projection of Radiomic Features")
plt.legend()
plt.tight_layout()
plt.show()

from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.metrics import silhouette_score
import matplotlib.pyplot as plt

lda = LinearDiscriminantAnalysis(n_components=1)   # 1D projection for binary labels
proj = lda.fit_transform(X, y)                     # shape (n_samples, 1)

# optional 2D by combining LDA + PCA
from sklearn.decomposition import PCA
pca  = PCA(n_components=1).fit_transform(X)
emb2 = np.hstack([proj, pca])   # 2D: [LDA, PCA]

print("LDA silhouette:", silhouette_score(emb2, y))
plt.figure()
plt.scatter(emb2[:,0], emb2[:,1], c=y, cmap="coolwarm", alpha=0.7)
# plt.title("LDA+PCA projection")
plt.show()

#%%

import numpy as np
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.decomposition import PCA
from sklearn.metrics import silhouette_score
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # needed for 3D plotting

# --- LDA (1D, since you said binary labels) ---
lda = LinearDiscriminantAnalysis(n_components=1)
proj = lda.fit_transform(X, y)   # shape (n_samples, 1)

# --- PCA (2D) ---
pca  = PCA(n_components=2).fit_transform(X)   # shape (n_samples, 2)

# --- Combine into 3D ---
emb3 = np.hstack([proj, pca])    # shape (n_samples, 3)

print("LDA+PCA silhouette:", silhouette_score(emb3, y))

# --- Plot 3D ---
fig = plt.figure(figsize=(6,6))
ax = fig.add_subplot(111, projection="3d")
sc = ax.scatter(emb3[:,0], emb3[:,1], emb3[:,2], c=y, cmap="coolwarm", alpha=0.7)
ax.set_xlabel("LDA axis", fontsize=12, labelpad=5)
ax.set_ylabel("PCA1 axis", fontsize=12, labelpad=5)
# ax.set_zlabel("PCA2 axis", labelpad=5)
ax.set_zlabel("")                  # hide default zlabel

# put a vertical label on the right side of the figure
fig.text(0.96, 0.5, "PCA2 axis", fontsize=12, rotation=270,
         va='center', ha='center')
# plt.title("LDA+PCA 3D Projection", fontsize=12)
# Adjust layout to prevent clipping
plt.tight_layout()
# Add more space on the right
plt.subplots_adjust(left=0.05, right=0.92, top=0.95, bottom=0.05)

# plt.title("LDA + PCA (3D Projection)")
plt.show()

#%% UMAP 3D
import umap
reducer = umap.UMAP(n_components=3, random_state=42, target_metric='categorical')
emb_super = reducer.fit_transform(X, y)   # supply labels
print("supervised UMAP silhouette:", silhouette_score(emb_super, y))
emb_unsup = reducer.fit_transform(X)
print("unsupervised UMAP silhouette:", silhouette_score(emb_unsup, y))

from mpl_toolkits.mplot3d import Axes3D

fig = plt.figure(figsize=(12,5))

# Supervised
ax = fig.add_subplot(121, projection='3d')
ax.scatter(emb_super[:,0], emb_super[:,1], emb_super[:,2], c=y, cmap='coolwarm', alpha=0.7)
ax.set_title("Supervised UMAP")

# Unsupervised
ax = fig.add_subplot(122, projection='3d')
ax.scatter(emb_unsup[:,0], emb_unsup[:,1], emb_unsup[:,2], c=y, cmap='coolwarm', alpha=0.7)
ax.set_title("Unsupervised UMAP")

plt.show()

#%% UMAP 2D

import umap
reducer = umap.UMAP(n_components=2, random_state=42, target_metric='categorical')
emb_super = reducer.fit_transform(X, y)   # supply labels
print("supervised UMAP silhouette:", silhouette_score(emb_super, y))
emb_unsup = reducer.fit_transform(X)
print("unsupervised UMAP silhouette:", silhouette_score(emb_unsup, y))

from mpl_toolkits.mplot3d import Axes3D

fig = plt.figure(figsize=(12,5))

# Supervised
ax = fig.add_subplot(121)
ax.scatter(emb_super[:,0], emb_super[:,1], c=y, cmap='coolwarm', alpha=0.7)
ax.set_title("Supervised UMAP")

# Unsupervised
ax = fig.add_subplot(122)
ax.scatter(emb_unsup[:,0], emb_unsup[:,1], c=y, cmap='coolwarm', alpha=0.7)
ax.set_title("Unsupervised UMAP")

plt.show()

#%% PCA
from sklearn.decomposition import PCA

# Fit PCA on your data
pca = PCA(n_components=3)
emb_pca = pca.fit_transform(X)   # shape (n_samples, 3)

print("Explained variance ratio:", pca.explained_variance_ratio_)

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

fig = plt.figure(figsize=(8,6))
ax = fig.add_subplot(111, projection='3d')

ax.scatter(
    emb_pca[:,0], emb_pca[:,1], emb_pca[:,2],
    c=y,                # color by class labels
    cmap='coolwarm',
    alpha=0.7
)

ax.set_xlabel('PC1')
ax.set_ylabel('PC2')
ax.set_zlabel('PC3')


ax.set_title('3D PCA Visualization')
plt.show()

score = silhouette_score(emb_pca, y)
print("PCA silhouette score:", score)

#%%

from sklearn.manifold import TSNE

# Convert X_data to numeric numpy array
X = np.array(X_data, dtype=np.float64)
y = np.array(merged_df['MGMT'], dtype=np.float64)  # binary labels

# Run t-SNE
tsne = TSNE(n_components=2, random_state=42, perplexity=30)
embedding = tsne.fit_transform(X)

# Plot
import matplotlib.pyplot as plt

plt.figure(figsize=(8,6))
for lab, color in zip([0,1], ["steelblue", "darkorange"]):
    idx = y == lab
    plt.scatter(embedding[idx,0], embedding[idx,1], c=color, label=f"Class {lab}", alpha=0.7, edgecolors="k", s=50)

plt.xlabel("t-SNE Dimension 1")
plt.ylabel("t-SNE Dimension 2")
plt.title("t-SNE projection of Radiomic Features")
plt.legend()
plt.tight_layout()
plt.show()

#%%

#%%

#%% Save the trained models

# best_model_tpot.export("best_tpot_pipeline.py")

# from best_tpot_pipeline import exported_pipeline
# Predict on new data
# preds = exported_pipeline.predict(X_new)
import joblib

joblib.dump(best_model_clf, "best_logistic_model.pkl")
joblib.dump(best_model_rf, "best_forest_model.pkl")

torch.save(model, "full_nn_model.pt")
# model = torch.load("full_nn_model.pt")
# model.eval()

torch.save(model.state_dict(), "best_nn_model.pt")
# model = BinaryRadiomicsNet(input_dim=X.shape[1])
# model.load_state_dict(torch.load("nn_state_dict.pt"))
# model.eval()

