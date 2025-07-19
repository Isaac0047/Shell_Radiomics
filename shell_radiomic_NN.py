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
    
    files_ce1 = [f for f in os.listdir(data_dir) if f.startswith(pattern1) and f.endswith('_CE_features_multi8_norm_update_l1.xlsx')]
    files_ce2 = [f for f in os.listdir(data_dir) if f.startswith(pattern2) and f.endswith('_CE_features_multi8_norm_update_l1.xlsx')]
    files_ce4 = [f for f in os.listdir(data_dir) if f.startswith(pattern4) and f.endswith('_CE_features_multi8_norm_update_l1.xlsx')]
    files_ce0 = [f for f in os.listdir(data_dir) if f.startswith(pattern0) and f.endswith('_CE_features_multi8_norm_update_l1.xlsx')]
    files_fl1 = [f for f in os.listdir(data_dir) if f.startswith(pattern1) and f.endswith('_FL_features_multi8_norm_update_l1.xlsx')]
    files_fl2 = [f for f in os.listdir(data_dir) if f.startswith(pattern2) and f.endswith('_FL_features_multi8_norm_update_l1.xlsx')]
    files_fl4 = [f for f in os.listdir(data_dir) if f.startswith(pattern4) and f.endswith('_FL_features_multi8_norm_update_l1.xlsx')]
    files_fl0 = [f for f in os.listdir(data_dir) if f.startswith(pattern0) and f.endswith('_FL_features_multi8_norm_update_l1.xlsx')]
    files_adc1 = [f for f in os.listdir(data_dir) if f.startswith(pattern1) and f.endswith('_ADC_features_multi8_norm_update_l1.xlsx')]
    files_adc2 = [f for f in os.listdir(data_dir) if f.startswith(pattern2) and f.endswith('_ADC_features_multi8_norm_update_l1.xlsx')]
    files_adc4 = [f for f in os.listdir(data_dir) if f.startswith(pattern4) and f.endswith('_ADC_features_multi8_norm_update_l1.xlsx')]
    files_adc0 = [f for f in os.listdir(data_dir) if f.startswith(pattern0) and f.endswith('_ADC_features_multi8_norm_update_l1.xlsx')]

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

# len_ce1t = []
# len_ce2t = []
# len_ce4t = []

# for i in range(len(stacked_data_ce1t)):
#     len_ce1t.append(stacked_data_ce1t[i].shape[0])
# for i in range(len(stacked_data_ce2t)):
#     len_ce2t.append(stacked_data_ce2t[i].shape[0])
# for i in range(len(stacked_data_ce4t)):
#     len_ce4t.append(stacked_data_ce4t[i].shape[0])

# min_len_ce1t = np.min(len_ce1t)
# min_len_ce2t = np.min(len_ce2t)
# min_len_ce4t = np.min(len_ce4t)

# print(min_len_ce1t)
# print(min_len_ce2t)
# print(min_len_ce4t)

#%%
import numpy as np

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

#%%
from collections import Counter

stacked_data_ce00  = [stacked_data_ce0[i]  for i in valid_indices]
stacked_data_ce11  = [stacked_data_ce1[i]  for i in valid_indices]
stacked_data_ce22  = [stacked_data_ce2[i]  for i in valid_indices]
stacked_data_ce44  = [stacked_data_ce4[i]  for i in valid_indices]
stacked_data_fl00  = [stacked_data_fl0[i]  for i in valid_indices]
stacked_data_fl11  = [stacked_data_fl1[i]  for i in valid_indices]
stacked_data_fl22  = [stacked_data_fl2[i]  for i in valid_indices]
stacked_data_fl44  = [stacked_data_fl4[i]  for i in valid_indices]
stacked_data_adc00 = [stacked_data_adc0[i] for i in valid_indices]
stacked_data_adc11 = [stacked_data_adc1[i] for i in valid_indices]
stacked_data_adc22 = [stacked_data_adc2[i] for i in valid_indices]
stacked_data_adc44 = [stacked_data_adc4[i] for i in valid_indices]

length = []
data_test = stacked_data_ce4

for i in range(len(data_test)):
    length.append(data_test[i].shape[1])

length_new = np.sort(length)
count_dict = Counter(length_new)
print(count_dict)

stacked_data_ce0_new  = np.array([matrix.flatten() for matrix in stacked_data_ce00])
stacked_data_fl0_new  = np.array([matrix.flatten() for matrix in stacked_data_fl00])
stacked_data_ce1_new  = np.array([matrix.flatten() for matrix in stacked_data_ce11])
stacked_data_fl1_new  = np.array([matrix.flatten() for matrix in stacked_data_fl11])
stacked_data_ce2_new  = np.array([matrix.flatten() for matrix in stacked_data_ce22])
stacked_data_fl2_new  = np.array([matrix.flatten() for matrix in stacked_data_fl22])
stacked_data_ce4_new  = np.array([matrix.flatten() for matrix in stacked_data_ce44])
stacked_data_fl4_new  = np.array([matrix.flatten() for matrix in stacked_data_fl44])
stacked_data_adc0_new = np.array([matrix.flatten() for matrix in stacked_data_adc00])
stacked_data_adc1_new = np.array([matrix.flatten() for matrix in stacked_data_adc11])
stacked_data_adc2_new = np.array([matrix.flatten() for matrix in stacked_data_adc22])
stacked_data_adc4_new = np.array([matrix.flatten() for matrix in stacked_data_adc44])

#%% Reform the matrix ########################################################################
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

min_max_scaler = MinMaxScaler()

ce0_mat_norm  = min_max_scaler.fit_transform(stacked_data_ce0_new)
fl0_mat_norm  = min_max_scaler.fit_transform(stacked_data_fl0_new)
ce1_mat_norm  = min_max_scaler.fit_transform(stacked_data_ce1_new)
fl1_mat_norm  = min_max_scaler.fit_transform(stacked_data_fl1_new)
ce2_mat_norm  = min_max_scaler.fit_transform(stacked_data_ce2_new)
fl2_mat_norm  = min_max_scaler.fit_transform(stacked_data_fl2_new)
ce4_mat_norm  = min_max_scaler.fit_transform(stacked_data_ce4_new)
fl4_mat_norm  = min_max_scaler.fit_transform(stacked_data_fl4_new)
adc0_mat_norm = min_max_scaler.fit_transform(stacked_data_adc0_new)
adc1_mat_norm = min_max_scaler.fit_transform(stacked_data_adc1_new)
adc2_mat_norm = min_max_scaler.fit_transform(stacked_data_adc2_new)
adc4_mat_norm = min_max_scaler.fit_transform(stacked_data_adc4_new)

#%%
radiomic_feature_online  = ce_list * 20
radiomic_feature_online1 = ce_list * 20
ce0_radiomic  = pd.DataFrame(ce0_mat_norm,  columns=radiomic_feature_online1)
fl0_radiomic  = pd.DataFrame(fl0_mat_norm,  columns=radiomic_feature_online1)
ce1_radiomic  = pd.DataFrame(ce1_mat_norm,  columns=radiomic_feature_online)
fl1_radiomic  = pd.DataFrame(fl1_mat_norm,  columns=radiomic_feature_online)
ce2_radiomic  = pd.DataFrame(ce2_mat_norm,  columns=radiomic_feature_online)
fl2_radiomic  = pd.DataFrame(fl2_mat_norm,  columns=radiomic_feature_online)
ce4_radiomic  = pd.DataFrame(ce4_mat_norm,  columns=radiomic_feature_online)
fl4_radiomic  = pd.DataFrame(fl4_mat_norm,  columns=radiomic_feature_online)
adc0_radiomic = pd.DataFrame(adc0_mat_norm,  columns=radiomic_feature_online1)
adc1_radiomic = pd.DataFrame(adc1_mat_norm,  columns=radiomic_feature_online1)
adc2_radiomic = pd.DataFrame(adc2_mat_norm,  columns=radiomic_feature_online1)
adc4_radiomic = pd.DataFrame(adc4_mat_norm,  columns=radiomic_feature_online1)

#%% Define the function to remove constant in the dataframe
from sklearn.feature_selection import VarianceThreshold
from statsmodels.stats.outliers_influence import variance_inflation_factor

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

#%% MGMT status check
gene_filter = gene[gene['Assigned ID'].isin(file_numbers_unique)]
gene_matrix = gene_filter.to_numpy()

# gene_matrix = gene_filter.to_numpy()
MGMT = gene_filter['MGMT']
gene_filter_1 = gene_filter
gene_filter_1['MGMT'] = gene_filter_1['MGMT'].map({'positive': 1, 'negative': 0})
gene_filter_1['EGFR'] = gene_filter_1['EGFR'].map(lambda x: 0 if x == 'negative' else 1)
#gene_filter_1['IDH']

# Apply the Min-Max normalization
def min_max_final(matrix):
    # row_min = np.min(matrix, axis=1)  # Minimum of each row
    # row_max = np.max(matrix, axis=1)  # Maximum of each row
    
    row_min = np.min(matrix, axis=0, keepdims=True)  # Minimum of each row
    row_max = np.max(matrix, axis=0, keepdims=True)  # Maximum of each row
    normalized_matrix = (matrix - row_min) / (row_max - row_min)
    
    return normalized_matrix

#%% Conduct CDKN2A/B Test
from scipy.stats import pointbiserialr
import seaborn as sns
# Predefine the dataframe needed
# OS (Days from diagnosis surgery)

gene_filter_1 = gene_filter_1[gene_filter_1['Assigned ID'].isin(common_indices)]
os_name       = gene_filter_1.columns[29]
status_name   = gene_filter_1.columns[28]

# Single
temp_df = ce1_radiomic_clean
temp_df = pd.merge(temp_df, ce2_radiomic_clean, left_index=True, right_index=True)
temp_df = pd.merge(temp_df, ce4_radiomic_clean, left_index=True, right_index=True)
temp_df = pd.merge(temp_df, ce0_radiomic_clean, left_index=True, right_index=True, suffixes=('','_ce0'))
# temp_df = pd.merge(temp_df, fl1_radiomic_clean, left_index=True, right_index=True, suffixes=('', '_fl1'))
# temp_df = pd.merge(temp_df, fl2_radiomic_clean, left_index=True, right_index=True, suffixes=('', '_fl2'))
# temp_df = pd.merge(temp_df, fl4_radiomic_clean, left_index=True, right_index=True, suffixes=('', '_fl4'))
# temp_df = pd.merge(temp_df, fl0_radiomic_clean, left_index=True, right_index=True, suffixes=('', '_fl0'))
temp_df = pd.merge(temp_df, adc1_radiomic_clean, left_index=True, right_index=True, suffixes=('', '_adc1'))
temp_df = pd.merge(temp_df, adc2_radiomic_clean, left_index=True, right_index=True, suffixes=('', '_adc2'))
temp_df = pd.merge(temp_df, adc4_radiomic_clean, left_index=True, right_index=True, suffixes=('', '_adc4'))
temp_df = pd.merge(temp_df, adc0_radiomic_clean, left_index=True, right_index=True, suffixes=('', '_adc0'))

# Summary
final_df  = temp_df
merged_df = pd.merge(final_df.reset_index(drop=True), gene_filter_1[['MGMT']].reset_index(drop=True), left_index=True, right_index=True)
merged_df = pd.merge(merged_df.reset_index(drop=True), gene_filter_1[['EGFR']].reset_index(drop=True), left_index=True, right_index=True)
merged_df = pd.merge(merged_df, gene_filter_1['IDH'].reset_index(drop=True), left_index=True, right_index=True)
merged_df = pd.merge(merged_df, gene_filter_1[os_name].reset_index(drop=True), left_index=True, right_index=True)
merged_df = pd.merge(merged_df, gene_filter_1[status_name].reset_index(drop=True), left_index=True, right_index=True)
merged_df = merged_df[merged_df['IDH'] == 'negative']

merged_df = merged_df.dropna(subset=['MGMT'])
merged_df = merged_df.dropna(subset=['EGFR'])
file_numbers_idh = merged_df.index.tolist()

file_numbers_common = [common_indices[i] for i in file_numbers_idh]

# Assuming you have 'df_radiomic' (radiomic features) and 'df_gene' (gene mutation status)
# Merge DataFrames on a common identifier, like 'patient_id'
# merged_df = pd.merge(final_df, gene_filter[['CDKN2A/B.1']], left_index=True, right_index=True)

# Extract radiomic features and CDKN2A/B loss column
# radiomic_features = merged_df.drop(['MGMT'], axis=1)
# cdkn2ab_status    = merged_df['MGMT']  # Binary status (0 or 1)

# Separate features and target variable
# X = merged_df.drop(['MGMT'], axis=1)  # Radiomic features
# y = merged_df['MGMT']  #== CDKN2A/B mutation status

#%% Change the gene name here
X    = merged_df.drop(['MGMT','EGFR','IDH',os_name,status_name], axis=1)
# y    = merged_df[status_name]
y    = merged_df['EGFR']

Labels = y.to_numpy()

# Try for predictions

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

X_numeric = X.to_numpy()
y_numeric = y.to_numpy()

# Convert your object array to float
X = np.array(X_numeric, dtype=np.float32)
y = np.array(y_numeric).astype(np.float32)

# Normalize features
scaler = StandardScaler()
X = scaler.fit_transform(X)

# selector = SelectKBest(score_func=f_classif, k=350)  # e.g., top 100
# X        = selector.fit_transform(X, y)

from sklearn.decomposition import PCA
# pca = PCA(n_components=60)  # or 20–50 depending on your variance explained
# X_scaled = pca.fit_transform(X_scaled)  # use scaled features
from sklearn.feature_selection import SelectKBest, f_classif

selector = SelectKBest(score_func=f_classif, k=350)  # e.g., top 100
X        = selector.fit_transform(X, y)

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
            nn.Linear(32, 1)  # 1 output for binary classification
            # nn.Sigmoid()       # Final sigmoid activation
        )

    def forward(self, x):
        return self.model(x)  # logits


device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
kfold  = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

auc_scores = []

for fold, (train_idx, test_idx) in enumerate(kfold.split(X, y)):
    print(f"\nFold {fold + 1}")

    X_train, X_test = torch.tensor(X[train_idx]), torch.tensor(X[test_idx])
    y_train, y_test = torch.tensor(y[train_idx]), torch.tensor(y[test_idx])

    model = BinaryRadiomicsNet(input_dim=X.shape[1]).to(device)
    # criterion = nn.BCELoss()
    criterion = nn.BCEWithLogitsLoss()
    optimizer = torch.optim.SGD(model.parameters(), lr=0.1, weight_decay=1e-3)
    # optimizer = torch.optim.Adam(model.parameters(), lr=1e-3, weight_decay=1e-3)

    # Training
    model.train()
    for epoch in range(8000):  # Adjust epochs
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
        auc_scores.append(auc)
        print(f"AUC for fold {fold + 1}: {auc:.4f}")

# Compute mean and 95% CI
mean_auc = np.mean(auc_scores)
std_auc = np.std(auc_scores, ddof=1)
ci95 = stats.t.interval(0.95, len(auc_scores) - 1, loc=mean_auc, scale=std_auc / np.sqrt(len(auc_scores)))

print(f"\nAverage AUC: {mean_auc:.4f}")
print(f"95% CI: ({ci95[0]:.4f}, {ci95[1]:.4f})")

#%% Logistic Regression

from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score

auc_scores = []
best_auc = -1
best_model = None
best_fold = -1

for fold, (train_idx, test_idx) in enumerate(kfold.split(X, y)):
    print(f"\nFold {fold + 1}")

    X_train, X_test = X[train_idx], X[test_idx]
    y_train, y_test = y[train_idx], y[test_idx]

    # Train logistic regression
    clf = LogisticRegression(penalty='l2', solver='liblinear', max_iter=1000)
    clf.fit(X_train, y_train)

    # Predict probabilities for class 1
    preds = clf.predict_proba(X_test)[:, 1]

    auc = roc_auc_score(y_test, preds)
    auc_scores.append(auc)
    print(f"AUC for fold {fold + 1}: {auc:.4f}")

    if auc > best_auc:
        best_auc = auc
        best_model_clf = clf
        best_fold = fold + 1

# After loop
mean_auc = np.mean(auc_scores)
std_auc = np.std(auc_scores, ddof=1)
ci95 = stats.t.interval(0.95, len(auc_scores) - 1, loc=mean_auc, scale=std_auc / np.sqrt(len(auc_scores)))

print(f"\nBest AUC: {best_auc:.4f} from Fold {best_fold}")
print(f"Average AUC: {mean_auc:.4f}")
print(f"95% CI: ({ci95[0]:.4f}, {ci95[1]:.4f})")

#%% Random Forest

from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score
from scipy import stats

auc_scores = []
best_auc = -1
best_model_clf = None
best_fold = -1

for fold, (train_idx, test_idx) in enumerate(kfold.split(X, y)):
    print(f"\nFold {fold + 1}")

    X_train, X_test = X[train_idx], X[test_idx]
    y_train, y_test = y[train_idx], y[test_idx]

    # Train Random Forest
    clf = RandomForestClassifier(n_estimators=100, max_depth=None, random_state=42)
    clf.fit(X_train, y_train)

    preds = clf.predict_proba(X_test)[:, 1]
    auc = roc_auc_score(y_test, preds)
    auc_scores.append(auc)
    print(f"AUC for fold {fold + 1}: {auc:.4f}")

    if auc > best_auc:
        best_auc = auc
        best_model_clf = clf
        best_fold = fold + 1

mean_auc = np.mean(auc_scores)
std_auc = np.std(auc_scores, ddof=1)
ci95 = stats.t.interval(0.95, len(auc_scores) - 1, loc=mean_auc, scale=std_auc / np.sqrt(len(auc_scores)))

print(f"\nBest AUC: {best_auc:.4f} from Fold {best_fold}")
print(f"Average AUC: {mean_auc:.4f}")
print(f"95% CI: ({ci95[0]:.4f}, {ci95[1]:.4f})")

#%% TPOT

from tpot import TPOTClassifier
from sklearn.metrics import roc_auc_score
from scipy import stats

auc_scores = []
best_auc = -1
best_model_tpot = None
best_fold = -1

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
