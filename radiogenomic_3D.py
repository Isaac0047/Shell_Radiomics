# This code performs the radiogenomic analysis

# The radiomic feature position
# Shape Feature: 0 - 13
# First Order Feature: 14 - 31
# GLCM: 32 - 55
# GLDM: 56 - 69
# GLRLM: 70 - 85
# GLSZM: 86 - 101
# NGTDM: 102 - 106

import numpy as np
import numpy as np
import matplotlib.pyplot as plt

import os
import pandas as pd

#%% Load existing files

# Load the clinical dataset
file_path = '/Users/hafeng/Documents/Research_Data/UCSF_Dataset/UCSF_PDGM/'
sheet_name = 'preop'

filename = 'preop_glioma_9-20-15--4-1-21-update-surv.xlsx'

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

#%% Check for the value of the existing files

# Specify the folder where your files are located
radiomic_path = '/Users/hafeng/Documents/Postdoc_Research_Meetings/TCGA/UCSF_PDGM/radiogenomic_analysis/'

# Initialize an empty list to store the integers
file_numbers = []

# Loop through all the files in the folder
for filename in os.listdir(radiomic_path):
    # Check if the filename is at least 4 characters long
    # print(filename)
    
    # if filename.endswith('_adc_combined_iso_256.xlsx') and filename[:4].isdigit():
    if filename.endswith('_ce4_filter_256.xlsx') and filename[:4].isdigit():
        # Extract the first 4 characters and convert to an integer
        file_number = int(filename[:4])
        # Append the integer to the list
        file_numbers.append(file_number)

# Display the resulting list of integers
# file_numbers_unique = list(set(file_numbers_adc))
file_numbers.sort()

print(file_numbers)
print(len(file_numbers))

#%% Extract the genomic matrix with the given unique values

gene_filter = gene[gene['Assigned ID'].isin(file_numbers)]
gene_matrix = gene_filter.to_numpy()

# os_time   = gene_fitler['']
# os_status = pd.read_excel(file_path+filename, sheet_name=sheet_name, usecols='AC')

#%% Load the radiomic features and stack into submatrices

# Initialize the data storage matrix
# ce1_matrix = np.zeros((145-38, gene_matrix.shape[0]))
# ce2_matrix = np.zeros((145-38, gene_matrix.shape[0]))
# ce3_matrix = np.zeros((145-38, gene_matrix.shape[0]))
# ce4_matrix = np.zeros((145-38, gene_matrix.shape[0]))

# fl1_matrix = np.zeros((145-38, gene_matrix.shape[0]))
# fl2_matrix = np.zeros((145-38, gene_matrix.shape[0]))
# fl3_matrix = np.zeros((145-38, gene_matrix.shape[0]))
# fl4_matrix = np.zeros((145-38, gene_matrix.shape[0]))

# adc1_matrix = np.zeros((145-38, gene_matrix.shape[0]))
# adc2_matrix = np.zeros((145-38, gene_matrix.shape[0]))
# adc3_matrix = np.zeros((145-38, gene_matrix.shape[0]))
# adc4_matrix = np.zeros((145-38, gene_matrix.shape[0]))

ce1_matrix = np.zeros((130-24+1, gene_matrix.shape[0]))
ce2_matrix = np.zeros((130-24+1, gene_matrix.shape[0]))
ce3_matrix = np.zeros((130-24+1, gene_matrix.shape[0]))
ce4_matrix = np.zeros((130-24+1, gene_matrix.shape[0]))

fl1_matrix = np.zeros((130-24+1, gene_matrix.shape[0]))
fl2_matrix = np.zeros((130-24+1, gene_matrix.shape[0]))
fl3_matrix = np.zeros((130-24+1, gene_matrix.shape[0]))
fl4_matrix = np.zeros((130-24+1, gene_matrix.shape[0]))

adc1_matrix = np.zeros((130-24+1, gene_matrix.shape[0]))
adc2_matrix = np.zeros((130-24+1, gene_matrix.shape[0]))
adc3_matrix = np.zeros((130-24+1, gene_matrix.shape[0]))
adc4_matrix = np.zeros((130-24+1, gene_matrix.shape[0]))

#%% Conduct data extraction
# subname_ce1  = '_ce1_filter_256.xlsx'
# subname_ce2  = '_ce2_filter_256.xlsx'
# subname_ce4  = '_ce4_filter_256.xlsx'
# subname_fl1  = '_fl1_filter_256.xlsx'
# subname_fl2  = '_fl2_filter_256.xlsx'
# subname_fl4  = '_fl4_filter_256.xlsx'
subname_ce1  = '_ce1.xlsx'
subname_ce2  = '_ce2.xlsx'
subname_ce4  = '_ce4.xlsx'
subname_fl1  = '_fl1.xlsx'
subname_fl2  = '_fl2.xlsx'
subname_fl4  = '_fl4.xlsx'
# subname_ce1  = '_ce1_filter.xlsx'
# subname_ce2  = '_ce2_filter.xlsx'
# subname_ce4  = '_ce4_filter.xlsx'
# subname_fl1  = '_fl1_filter.xlsx'
# subname_fl2  = '_fl2_filter.xlsx'
# subname_fl4  = '_fl4_filter.xlsx'
subname_adc1 = '_adc1.xlsx'
subname_adc2 = '_adc2.xlsx'
subname_adc4 = '_adc4.xlsx'

skip_row = 23

for i in range(len(file_numbers)):
    
    idx = file_numbers[i]
    folder_name = f'{idx:04d}'
    
    ce1_matrix[:,i]  = np.squeeze(pd.read_excel(radiomic_path+folder_name+subname_ce1, sheet_name='Sheet1', skiprows=skip_row, usecols='B', header=None).to_numpy())
    ce2_matrix[:,i]  = np.squeeze(pd.read_excel(radiomic_path+folder_name+subname_ce2, sheet_name='Sheet1', skiprows=skip_row, usecols='B', header=None).to_numpy())
    ce4_matrix[:,i]  = np.squeeze(pd.read_excel(radiomic_path+folder_name+subname_ce4, sheet_name='Sheet1', skiprows=skip_row, usecols='B', header=None).to_numpy())
    
    fl1_matrix[:,i]  = np.squeeze(pd.read_excel(radiomic_path+folder_name+subname_fl1, sheet_name='Sheet1', skiprows=skip_row, usecols='B', header=None).to_numpy())
    fl2_matrix[:,i]  = np.squeeze(pd.read_excel(radiomic_path+folder_name+subname_fl2, sheet_name='Sheet1', skiprows=skip_row, usecols='B', header=None).to_numpy())
    fl4_matrix[:,i]  = np.squeeze(pd.read_excel(radiomic_path+folder_name+subname_fl4, sheet_name='Sheet1', skiprows=skip_row, usecols='B', header=None).to_numpy())
    
    adc1_matrix[:,i]  = np.squeeze(pd.read_excel(radiomic_path+folder_name+subname_adc1, sheet_name='Sheet1', skiprows=skip_row, usecols='B', header=None).to_numpy())
    adc2_matrix[:,i]  = np.squeeze(pd.read_excel(radiomic_path+folder_name+subname_adc2, sheet_name='Sheet1', skiprows=skip_row, usecols='B', header=None).to_numpy())
    adc4_matrix[:,i]  = np.squeeze(pd.read_excel(radiomic_path+folder_name+subname_adc4, sheet_name='Sheet1', skiprows=skip_row, usecols='B', header=None).to_numpy())
#%% Load the radiomic features and stack into submatrices

# Initialize the data storage matrix
# ce_matrix   = np.zeros((130-24+1, gene_matrix.shape[0]))
# fl_matrix   = np.zeros((130-24+1, gene_matrix.shape[0]))
# adc_matrix  = np.zeros((130-24+1, gene_matrix.shape[0]))
# pre_matrix  = np.zeros((130-24+1, gene_matrix.shape[0]))
# subname_ce  = '_ce_combined.xlsx'
# subname_fl  = '_fl_combined.xlsx'
# subname_adc = '_adc_combined.xlsx'
# skip_row    = 23

# ce_matrix  = np.zeros((145-38, gene_matrix.shape[0]))
# fl_matrix  = np.zeros((145-38, gene_matrix.shape[0]))
# adc_matrix = np.zeros((145-38, gene_matrix.shape[0]))
# subname_ce  = '_ce_combined_preCrop.xlsx'
# subname_fl  = '_fl_combined_preCrop.xlsx'
# subname_adc = '_adc_combined_preCrop.xlsx'
# # subname_ce  = '_ce_combined_iso.xlsx'
# # subname_fl  = '_fl_combined_iso.xlsx'
# # subname_adc = '_adc_combined_iso.xlsx'
# skip_row    = 23

# for i in range(len(file_numbers)):
    
#     idx = file_numbers[i]
#     folder_name = f'{idx:04d}'
    
#     ce_matrix[:,i]  = np.squeeze(pd.read_excel(radiomic_path+folder_name+subname_ce,  sheet_name='Sheet1', skiprows=skip_row, usecols='B', header=None).to_numpy())
#     fl_matrix[:,i]  = np.squeeze(pd.read_excel(radiomic_path+folder_name+subname_fl,  sheet_name='Sheet1', skiprows=skip_row, usecols='B', header=None).to_numpy())
#     adc_matrix[:,i] = np.squeeze(pd.read_excel(radiomic_path+folder_name+subname_adc, sheet_name='Sheet1', skiprows=skip_row, usecols='B', header=None).to_numpy())
  
#%% Load Radiomic titles    
ce1_df = pd.read_excel(radiomic_path+folder_name+subname_ce1, sheet_name='Sheet1', skiprows=23, usecols='A', header=None)
ce1_list = ce1_df.iloc[:,0].tolist()
ce2_df = pd.read_excel(radiomic_path+folder_name+subname_ce2, sheet_name='Sheet1', skiprows=23, usecols='A', header=None)
ce2_list = ce2_df.iloc[:,0].tolist()
ce4_df = pd.read_excel(radiomic_path+folder_name+subname_ce4, sheet_name='Sheet1', skiprows=23, usecols='A', header=None)
ce4_list = ce4_df.iloc[:,0].tolist()

fl1_df = pd.read_excel(radiomic_path+folder_name+subname_fl1, sheet_name='Sheet1', skiprows=23, usecols='A', header=None)
fl1_list = fl1_df.iloc[:,0].tolist()
fl2_df = pd.read_excel(radiomic_path+folder_name+subname_fl2, sheet_name='Sheet1', skiprows=23, usecols='A', header=None)
fl2_list = fl2_df.iloc[:,0].tolist()
fl4_df = pd.read_excel(radiomic_path+folder_name+subname_fl4, sheet_name='Sheet1', skiprows=23, usecols='A', header=None)
fl4_list = fl4_df.iloc[:,0].tolist()

adc1_df = pd.read_excel(radiomic_path+folder_name+subname_adc1, sheet_name='Sheet1', skiprows=23, usecols='A', header=None)
adc1_list = adc1_df.iloc[:,0].tolist()
adc2_df = pd.read_excel(radiomic_path+folder_name+subname_adc2, sheet_name='Sheet1', skiprows=23, usecols='A', header=None)
adc2_list = adc2_df.iloc[:,0].tolist()
adc4_df = pd.read_excel(radiomic_path+folder_name+subname_adc4, sheet_name='Sheet1', skiprows=23, usecols='A', header=None)
adc4_list = adc4_df.iloc[:,0].tolist()

ce_df    = pd.read_excel(radiomic_path+folder_name+subname_ce, sheet_name='Sheet1', skiprows=skip_row, usecols='A', header=None)
ce_list  = ce_df.iloc[:,0].tolist()
fl_df    = pd.read_excel(radiomic_path+folder_name+subname_fl, sheet_name='Sheet1', skiprows=skip_row, usecols='A', header=None)
fl_list  = fl_df.iloc[:,0].tolist()
adc_df   = pd.read_excel(radiomic_path+folder_name+subname_adc, sheet_name='Sheet1', skiprows=skip_row, usecols='A', header=None)
adc_list = adc_df.iloc[:,0].tolist()

#%% Normalize the radiomic data

# Apply the Min-Max normalization
def min_max_norm(matrix):
    # row_min = np.min(matrix, axis=1)  # Minimum of each row
    # row_max = np.max(matrix, axis=1)  # Maximum of each row
    
    row_min = np.min(matrix, axis=1, keepdims=True)  # Minimum of each row
    row_max = np.max(matrix, axis=1, keepdims=True)  # Maximum of each row
    
    normalized_matrix = (matrix - row_min) / (row_max - row_min)
    
    return normalized_matrix

def min_max_norm_df(df):
    
    # Function for min-max normalization
    row_min = df.min(axis=1)  # Minimum of each row
    row_max = df.max(axis=1)  # Maximum of each row
    
    # Perform min-max normalization for each row
    normalized_df = (df.subtract(row_min, axis=0)).div(row_max - row_min, axis=0)
    
    return normalized_df

# df_ce1_mat_norm = min_max_norm_df(ce1_df)
# df_ce2_mat_norm = min_max_norm_df(ce2_df)
# df_ce3_mat_norm = min_max_norm_df(ce3_df)

# df_fl1_mat_norm = min_max_norm_df(fl1_df)
# df_fl2_mat_norm = min_max_norm_df(fl2_df)
# df_fl3_mat_norm = min_max_norm_df(fl3_df)

ce1_mat_norm = min_max_norm(ce1_matrix)
ce2_mat_norm = min_max_norm(ce2_matrix)
ce4_mat_norm = min_max_norm(ce4_matrix)

fl1_mat_norm = min_max_norm(fl1_matrix)
fl2_mat_norm = min_max_norm(fl2_matrix)
fl4_mat_norm = min_max_norm(fl4_matrix)

adc1_mat_norm = min_max_norm(adc1_matrix)
adc2_mat_norm = min_max_norm(adc2_matrix)
adc4_mat_norm = min_max_norm(adc4_matrix)

# pre1_mat_norm = min_max_norm(pre1_matrix)
# pre2_mat_norm = min_max_norm(pre2_matrix)
# pre4_mat_norm = min_max_norm(pre4_matrix)

# ce_mat_norm  = min_max_norm(ce_matrix)
# fl_mat_norm  = min_max_norm(fl_matrix)
# adc_mat_norm = min_max_norm(adc_matrix)

#%% Plot the box plot across different dataset

# The radiomic feature position
# Shape Feature: 0 - 13
# First Order Feature: 14 - 31
# GLCM: 32 - 55
# GLDM: 56 - 69
# GLRLM: 70 - 85
# GLSZM: 86 - 101
# NGTDM: 102 - 106

def matrix_indice(mat, row_s, row_e):
    return mat[row_s:row_e,:]

row_s = 32
row_e = 55

# Original matrices (replace these with your actual matrices)
matrix1 = ce1_mat_norm[row_s:row_e,:]
matrix2 = ce2_mat_norm[row_s:row_e,:]
matrix3 = ce4_mat_norm[row_s:row_e,:]
matrix4 = fl1_mat_norm[row_s:row_e,:]
matrix5 = fl2_mat_norm[row_s:row_e,:]
matrix6 = fl4_mat_norm[row_s:row_e,:]
matrix7 = adc1_mat_norm[row_s:row_e,:]
matrix8 = adc2_mat_norm[row_s:row_e,:]
matrix9 = adc4_mat_norm[row_s:row_e,:]

# Normalized matrices
matrix1 = ce1_matrix[row_s:row_e,:]
matrix2 = ce2_matrix[row_s:row_e,:]
matrix3 = ce4_matrix[row_s:row_e,:]
matrix4 = fl1_matrix[row_s:row_e,:]
matrix5 = fl2_matrix[row_s:row_e,:]
matrix6 = fl4_matrix[row_s:row_e,:]
matrix7 = adc1_matrix[row_s:row_e,:]
matrix8 = adc2_matrix[row_s:row_e,:]
matrix9 = adc4_matrix[row_s:row_e,:]

# Combine the matrices in a list

modality = 2

if modality == 1:
    matrices   = [matrix1, matrix2, matrix3]
    title_name = 'Box Plot for Each Row Across Matrices for T1ce'   
elif modality == 2:
    matrices = [matrix4, matrix5, matrix6]
    title_name = 'Box Plot for Each Row Across Matrices for T2flair'
elif modality == 3:
    matrices = [matrix7, matrix8, matrix9]
    title_name = 'Box Plot for Each Row Across Matrices for ADC'

# Get the number of rows
num_rows = matrix1.shape[0]

# Prepare the data for plotting (we'll flatten everything in a single axis)
box_plot_data = []
labels = []

# Iterate through each row and collect data for box plots
for i in range(num_rows):
    # Extract the row from each matrix and stack them into a list for the boxplot
    row_data = [matrix[i, :] for matrix in matrices]
    box_plot_data.extend(row_data)  # Combine all rows from matrices into one list

    # Create labels for each box
    labels.extend([f'Row {i+1} - Matrix 1', f'Row {i+1} - Matrix 2', f'Row {i+1} - Matrix 3'])

# Create the boxplot in a single figure
plt.figure(figsize=(12, 6))

plt.boxplot(box_plot_data, labels=labels, patch_artist=True)
plt.xticks(rotation=90)  # Rotate x labels for better readability
plt.title(title_name)
plt.tight_layout()

# Show the plot
plt.show()

#%% Conduct survival analysis for specific features

from lifelines import CoxPHFitter
from lifelines.statistics import logrank_test
from lifelines import KaplanMeierFitter

radiomic_feature_df          = pd.read_excel(radiomic_path+folder_name+subname_ce1, sheet_name='Sheet1', skiprows=23, usecols='A', header=None)
radiomic_feature_df_online   = pd.read_excel(radiomic_path+folder_name+subname_ce,  sheet_name='Sheet1', skiprows=skip_row, usecols='A', header=None)
radiomic_feature_list        = radiomic_feature_df.iloc[:,0].tolist()
radiomic_feature_online      = radiomic_feature_df_online.iloc[:,0].tolist()

#%% define all the dataframes
# ce1_radiomic = pd.DataFrame(ce1_matrix.T, columns=radiomic_feature_list)
# ce2_radiomic = pd.DataFrame(ce2_matrix.T, columns=radiomic_feature_list)
# ce4_radiomic = pd.DataFrame(ce4_matrix.T, columns=radiomic_feature_list)
# fl1_radiomic = pd.DataFrame(fl1_matrix.T, columns=radiomic_feature_list)
# fl2_radiomic = pd.DataFrame(fl2_matrix.T, columns=radiomic_feature_list)
# fl4_radiomic = pd.DataFrame(fl4_matrix.T, columns=radiomic_feature_list)
# adc1_radiomic = pd.DataFrame(adc1_matrix.T, columns=radiomic_feature_list)
# adc2_radiomic = pd.DataFrame(adc2_matrix.T, columns=radiomic_feature_list)
# adc4_radiomic = pd.DataFrame(adc4_matrix.T, columns=radiomic_feature_list)
# pre1_radiomic = pd.DataFrame(pre1_matrix.T, columns=radiomic_feature_list)
# pre2_radiomic = pd.DataFrame(pre2_matrix.T, columns=radiomic_feature_list)
# pre4_radiomic = pd.DataFrame(pre4_matrix.T, columns=radiomic_feature_list)

# ce_radiomic  = pd.DataFrame(ce_matrix.T,  columns=radiomic_feature_list)
# fl_radiomic  = pd.DataFrame(fl_matrix.T,  columns=radiomic_feature_list)
# adc_radiomic = pd.DataFrame(adc_matrix.T, columns=radiomic_feature_list)

ce1_radiomic = pd.DataFrame(ce1_mat_norm.T, columns=radiomic_feature_list)
ce2_radiomic = pd.DataFrame(ce2_mat_norm.T, columns=radiomic_feature_list)
ce4_radiomic = pd.DataFrame(ce4_mat_norm.T, columns=radiomic_feature_list)
fl1_radiomic = pd.DataFrame(fl1_mat_norm.T, columns=radiomic_feature_list)
fl2_radiomic = pd.DataFrame(fl2_mat_norm.T, columns=radiomic_feature_list)
fl4_radiomic = pd.DataFrame(fl4_mat_norm.T, columns=radiomic_feature_list)
adc1_radiomic = pd.DataFrame(adc1_mat_norm.T, columns=radiomic_feature_list)
adc2_radiomic = pd.DataFrame(adc2_mat_norm.T, columns=radiomic_feature_list)
adc4_radiomic = pd.DataFrame(adc4_mat_norm.T, columns=radiomic_feature_list)
pre1_radiomic = pd.DataFrame(pre1_mat_norm.T, columns=radiomic_feature_list)
pre2_radiomic = pd.DataFrame(pre2_mat_norm.T, columns=radiomic_feature_list)
pre4_radiomic = pd.DataFrame(pre4_mat_norm.T, columns=radiomic_feature_list)

# ce_radiomic  = pd.DataFrame(ce_mat_norm.T,  columns=radiomic_feature_list)
# fl_radiomic  = pd.DataFrame(fl_mat_norm.T,  columns=radiomic_feature_list)
# adc_radiomic = pd.DataFrame(adc_mat_norm.T, columns=radiomic_feature_list)

ce_radiomic  = pd.DataFrame(ce_mat_norm.T,  columns=radiomic_feature_online)
fl_radiomic  = pd.DataFrame(fl_mat_norm.T,  columns=radiomic_feature_online)
adc_radiomic = pd.DataFrame(adc_mat_norm.T, columns=radiomic_feature_online)

#%% Add survival and duration to the dataframe
# ce1_radiomic['survival_time'] = np.squeeze(os_time_matrix)
# ce1_radiomic['event']         = np.squeeze(os_status_matrix)
# ce2_radiomic['survival_time'] = np.squeeze(os_time_matrix)
# ce2_radiomic['event']         = np.squeeze(os_status_matrix)
# # ce3_radiomic['survival_time'] = np.squeeze(os_time_matrix)
# # ce3_radiomic['event']         = np.squeeze(os_status_matrix)
# ce4_radiomic['survival_time'] = np.squeeze(os_time_matrix)
# ce4_radiomic['event']         = np.squeeze(os_status_matrix)

# fl1_radiomic['survival_time'] = np.squeeze(os_time_matrix)
# fl1_radiomic['event']         = np.squeeze(os_status_matrix)
# fl2_radiomic['survival_time'] = np.squeeze(os_time_matrix)
# fl2_radiomic['event']         = np.squeeze(os_status_matrix)
# # fl3_radiomic['survival_time'] = np.squeeze(os_time_matrix)
# # fl3_radiomic['event']         = np.squeeze(os_status_matrix)
# fl4_radiomic['survival_time'] = np.squeeze(os_time_matrix)
# fl4_radiomic['event']         = np.squeeze(os_status_matrix)

# adc1_radiomic['survival_time'] = np.squeeze(os_time_matrix)
# adc1_radiomic['event']         = np.squeeze(os_status_matrix)
# adc2_radiomic['survival_time'] = np.squeeze(os_time_matrix)
# adc2_radiomic['event']         = np.squeeze(os_status_matrix)
# # adc3_radiomic['survival_time'] = np.squeeze(os_time_matrix)
# # adc3_radiomic['event']         = np.squeeze(os_status_matrix)
# adc4_radiomic['survival_time'] = np.squeeze(os_time_matrix)
# adc4_radiomic['event']         = np.squeeze(os_status_matrix)

# pre1_radiomic['survival_time'] = np.squeeze(os_time_matrix)
# pre1_radiomic['event']         = np.squeeze(os_status_matrix)
# pre2_radiomic['survival_time'] = np.squeeze(os_time_matrix)
# pre2_radiomic['event']         = np.squeeze(os_status_matrix)
# # pre3_radiomic['survival_time'] = np.squeeze(os_time_matrix)
# # pre3_radiomic['event']         = np.squeeze(os_status_matrix)
# pre4_radiomic['survival_time'] = np.squeeze(os_time_matrix)
# pre4_radiomic['event']         = np.squeeze(os_status_matrix)

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

def VT_remove(df, threshold=0.01):
    
    # features_df = df.drop(columns=['survival_time', 'event'])
    selector = VarianceThreshold(threshold=threshold)
    selector.fit(features_df)
    # print(len(selector.get_support()))
    
    columns_to_keep = features_df.columns[selector.get_support()]
    reduced_features_df = pd.DataFrame(selector.transform(features_df), columns=columns_to_keep)
    print('Reduced Feature Shape:',  reduced_features_df)
    
    final_df = pd.concat([reduced_features_df, df[['survival_time', 'event']]], axis=1)
    
    return final_df

def VIF_remove(df, threshold=10):

    # Step 1: Create a copy of the dataframe to avoid modifying the original one
    df_copy = df
    
    # Step 2: Iterate until all VIF values are below the threshold
    while True:
        # Calculate VIF for each feature
        vif = pd.DataFrame()
        vif["features"] = df_copy.columns
        vif["VIF"] = [variance_inflation_factor(df_copy.values, i) for i in range(df_copy.shape[1])]

        # Check the max VIF value
        max_vif = vif["VIF"].max()

        if max_vif > threshold:
            # Find the feature with the highest VIF and drop it
            feature_to_drop = vif.loc[vif["VIF"].idxmax(), "features"]
            print(f"Dropping feature '{feature_to_drop}' with VIF: {max_vif}")
            df_copy = df_copy.drop(columns=[feature_to_drop])
        else:
            # All VIF values are below the threshold
            break

    df_copy['survival_time'] = df['survival_time']
    df_copy['event']         = df['event']
    
    return df_copy

#%% Clean the data storage

#
# df_radiomic_clean_var = var_remove(df_radiomic_clean, 0.01)
# df_gene_clean_var     = var_remove(df_gene_clean, 0.01) 

#
# df_radiomic_clean_vt  = VT_remove(df_radiomic_clean, 0.03)
# df_gene_clean_vt      = VT_remove(df_gene_clean, 0.03)


#%% Remove constant and collineary parameters

# gene_clean         = const_remove(df_gene) 
# gene_clean_vif     = VIF_remove(gene_clean, 10)

#%%
ce1_radiomic_clean     = const_remove(ce1_radiomic)
# ce1_radiomic_clean_vif = VIF_remove(ce1_radiomic_clean, 10)
ce2_radiomic_clean     = const_remove(ce2_radiomic)
# ce2_radiomic_clean_vif = VIF_remove(ce2_radiomic_clean, 10)
# ce3_radiomic_clean     = const_remove(ce3_radiomic)
# ce3_radiomic_clean_vif = VIF_remove(ce3_radiomic_clean, 10)
ce4_radiomic_clean     = const_remove(ce4_radiomic)
# ce4_radiomic_clean_vif = VIF_remove(ce4_radiomic_clean, 10)

fl1_radiomic_clean     = const_remove(fl1_radiomic)
# fl1_radiomic_clean_vif = VIF_remove(fl1_radiomic_clean, 10)
fl2_radiomic_clean     = const_remove(fl2_radiomic)
# fl2_radiomic_clean_vif = VIF_remove(fl2_radiomic_clean, 10)
# fl3_radiomic_clean     = const_remove(fl3_radiomic)
# fl3_radiomic_clean_vif = VIF_remove(fl3_radiomic_clean, 10)
fl4_radiomic_clean     = const_remove(fl4_radiomic)
# fl4_radiomic_clean_vif = VIF_remove(fl4_radiomic_clean, 10)

adc1_radiomic_clean     = const_remove(adc1_radiomic)
# adc1_radiomic_clean_vif = VIF_remove(adc1_radiomic_clean, 10)
adc2_radiomic_clean     = const_remove(adc2_radiomic)
# adc2_radiomic_clean_vif = VIF_remove(adc2_radiomic_clean, 10)
# adc3_radiomic_clean     = const_remove(adc3_radiomic)
# adc3_radiomic_clean_vif = VIF_remove(adc3_radiomic_clean, 10)
adc4_radiomic_clean     = const_remove(adc4_radiomic)
# adc4_radiomic_clean_vif = VIF_remove(adc4_radiomic_clean, 10)

# pre1_radiomic_clean     = const_remove(pre1_radiomic)
# pre1_radiomic_clean_vif = VIF_remove(pre1_radiomic_clean, 10)
# pre2_radiomic_clean     = const_remove(pre2_radiomic)
# pre2_radiomic_clean_vif = VIF_remove(pre2_radiomic_clean, 10)
# pre3_radiomic_clean     = const_remove(pre3_radiomic)
# pre3_radiomic_clean_vif = VIF_remove(pre3_radiomic_clean, 10)
# pre4_radiomic_clean     = const_remove(pre4_radiomic)
# pre4_radiomic_clean_vif = VIF_remove(pre4_radiomic_clean, 10)

ce_radiomic_clean   = const_remove(ce_radiomic)
fl_radiomic_clean   = const_remove(fl_radiomic)
adc_radiomic_clean  = const_remove(adc_radiomic)

#%% MGMT status check

# def MGMT_process(df):
gene_filter = gene[gene['Assigned ID'].isin(file_numbers)]
# gene_matrix = gene_filter.to_numpy()
MGMT = gene_filter['MGMT']
gene_filter_1 = gene_filter
gene_filter_1['MGMT'] = gene_filter['MGMT'].map({'positive': 1, 'negative': 0})
gene_filter_1['EGFR'] = gene_filter['EGFR'].apply(lambda x: 0 if x == 'negative' else 1)

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

target_gene = 'EGFR'
# Tria

ce_radiomic = pd.merge(ce1_radiomic_clean, ce2_radiomic_clean, left_index=True, right_index=True)
ce_radiomic = pd.merge(ce_radiomic, ce4_radiomic_clean, left_index=True, right_index=True)

fl_radiomic = pd.merge(fl1_radiomic_clean, fl2_radiomic_clean, left_index=True, right_index=True)
fl_radiomic = pd.merge(fl_radiomic_clean,  fl4_radiomic_clean, left_index=True, right_index=True)

adc_radiomic = pd.merge(adc1_radiomic_clean, adc2_radiomic_clean, left_index=True, right_index=True)
adc_radiomic = pd.merge(adc_radiomic, adc4_radiomic_clean, left_index=True, right_index=True)

# Single
temp_df = ce_radiomic
# temp_df = ce4_radiomic_clean

# Double
# temp_df = pd.merge(ce1_radiomic_clean, ce2_radiomic_clean, left_index=True, right_index=True)

# Tria
# temp_df = pd.merge(fl1_radiomic_clean, fl2_radiomic_clean, left_index=True, right_index=True)
# temp_df = pd.merge(temp_df, fl4_radiomic_clean, left_index=True, right_index=True)

# All
# temp_df = pd.merge(ce1_radiomic_clean, ce2_radiomic_clean, left_index=True, right_index=True)
# temp_df = pd.merge(temp_df, ce4_radiomic_clean.add_suffix('_4'), left_index=True, right_index=True)
# temp_df = pd.merge(temp_df, fl1_radiomic_clean.add_suffix('_fl'), left_index=True, right_index=True)
# temp_df = pd.merge(temp_df, fl2_radiomic_clean.add_suffix('_fl'), left_index=True, right_index=True)
# temp_df = pd.merge(temp_df, fl4_radiomic_clean.add_suffix('_fl4'), left_index=True, right_index=True)

# Summary
final_df  = temp_df
merged_df = pd.merge(final_df, gene_filter_1[[target_gene]], left_index=True, right_index=True)
merged_df = merged_df.dropna(subset=[target_gene])

# Assuming you have 'df_radiomic' (radiomic features) and 'df_gene' (gene mutation status)
# Merge DataFrames on a common identifier, like 'patient_id'
# merged_df = pd.merge(final_df, gene_filter[['CDKN2A/B.1']], left_index=True, right_index=True)

# # Extract radiomic features and CDKN2A/B loss column
# radiomic_features = merged_df.drop(['MGMT'], axis=1)
# cdkn2ab_status    = merged_df['MGMT']  # Binary status (0 or 1)

# # Calculate point-biserial correlation between each radiomic feature and CDKN2A/B loss
# correlations = {}
# for feature in radiomic_features.columns:
#     corr, p_value = pointbiserialr(cdkn2ab_status, merged_df[feature])
#     correlations[feature] = corr

# # Convert the dictionary to a DataFrame for easier manipulation
# correlation_df = pd.DataFrame(list(correlations.items()), columns=['Feature', 'Correlation'])

# # Visualize the correlation using a bar plot
# plt.figure(figsize=(10, 8))
# sns.barplot(x='Correlation', y='Feature', data=correlation_df.sort_values(by='Correlation', ascending=False))
# plt.title('Correlation Between Radiomic Features and CDKN2A/B Loss')
# plt.show()

#%% Create New DataFrame

# Separate features and target variable
X = merged_df.drop([target_gene], axis=1)  # Radiomic features
y = merged_df[target_gene]  #== CDKN2A/B mutation status

# Copy the original DataFrame
new_df = merged_df.drop([target_gene], axis=1).copy()
# Add 0.15 to all features if the last column (Label) is equal to 1
new_df.loc[y == 1] += 0

#%% Select balanced labels
Labels = y.to_numpy()

# Step 1: Separate Positive and Negative Samples
positive_indices = np.where(Labels == 1)[0]  # Get indices of positive samples
negative_indices = np.where(Labels == 0)[0]  # Get indices of negative samples

# Step 2: Sample from Positive Indices
num_negative = len(negative_indices)
random_positive_indices = np.random.choice(positive_indices, size=num_negative, replace=False)

# Step 3: Combine the Balanced Dataset
balanced_indices = np.concatenate([negative_indices, random_positive_indices])
np.random.shuffle(balanced_indices)  # Shuffle the indices

# Step 4: Create Balanced Dataset
X_bal = X.iloc[balanced_indices]  # Select rows by balanced indices
y_bal = y.iloc[balanced_indices]  # Select corresponding labels

#%% Regression

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
#
# X = new_df
# X = X.iloc[:,14:31]

X_numeric = X.to_numpy()  # or X.values
y_numeric = y.to_numpy()   # or Y.values

# X_numeric = X_numeric[:,32:35]

# X_numeric = min_max_final(X_numeric)

# Define the number of folds for k-fold cross-validation
k = 5
# kf = KFold(n_splits=k, shuffle=True, random_state=42)
kf = StratifiedKFold(n_splits=k, shuffle=True, random_state=42)

smote = SMOTE(random_state=42)
# X_resampled, y_resampled = smote.fit_resample(X_numeric, y_numeric)

# Split the data into training and testing sets
# X_train, X_test, y_train, y_test = train_test_split(X_numeric, y_numeric, test_size=0.2, random_state=42)

########### Test models
# Store accuracy scores for each fold
accuracy_scores = []
ruc_scores      = []
accuracy_train  = []
ruc_train       = []

# Perform k-fold cross-validation
for train_index, test_index in kf.split(X_numeric, y_numeric):
    # X_train, X_test = X.iloc[train_index], X.iloc[test_index]
    # X_train, X_test = X_numeric[train_index], X_numeric[test_index]
    X_train, X_test = X_numeric[train_index], X_numeric[test_index]
    y_train, y_test = y_numeric[train_index], y_numeric[test_index]

    # Train the model using the cleaned training set
    # model = LogisticRegression(max_iter=1000, solver='saga', penalty='elasticnet', l1_ratio=0.8, C=1000)  # Use any model you prefer
    model = LogisticRegression(max_iter=1000, solver='lbfgs', penalty='l2')
    # model = RandomForestClassifier()
    # model = GradientBoostingClassifier(n_estimators=100, learning_rate=0.1, max_depth=3)
    X_resampled, y_resampled = smote.fit_resample(X_train, y_train)
    # model.fit(X_resampled, y_resampled)
    model.fit(X_train, y_train) 

    # Make predictions on the cleaned test set
    y_pred = model.predict(X_test)
    y_pred_train = model.predict(X_train)

    # Calculate accuracy
    accuracy = accuracy_score(y_test, y_pred)
    ruc_acc  = roc_auc_score(y_test, y_pred)
    accuracy_scores.append(accuracy)
    ruc_scores.append(ruc_acc)
    
    accuracy_t = accuracy_score(y_train, y_pred_train)
    ruc_acc_t  = roc_auc_score(y_train, y_pred_train)
    accuracy_train.append(accuracy_t)
    ruc_train.append(ruc_acc_t)

# Calculate the average accuracy over all folds
average_accuracy = np.mean(accuracy_scores)
average_accuracy_train = np.mean(accuracy_train)
print(f"Average Test accuracy over {k} folds: {average_accuracy:.4f}")
print(f"Average Train accuracy over {k} folds: {average_accuracy_train:.4f}")

average_ruc = np.mean(ruc_scores)
average_ruc_train = np.mean(ruc_train)
print(f"Average Test RUC over {k} folds: {average_ruc:.4f}")
print(f"Average Train RUC over {k} folds: {average_ruc_train:.4f}")

#%% Plot the box plot across different dataset

# The radiomic feature position
# Shape Feature: 0 - 13
# First Order Feature: 14 - 31
# GLCM: 32 - 55
# GLDM: 56 - 69
# GLRLM: 70 - 85
# GLSZM: 86 - 101
# NGTDM: 102 - 106

def matrix_indice(mat, row_s, row_e):
    return mat[row_s:row_e,:]

row_s = 32
row_e = 55

# Original matrices (replace these with your actual matrices)
# matrix1 = ce1_mat_norm[row_s:row_e,:]
# matrix2 = ce2_mat_norm[row_s:row_e,:]
# matrix3 = ce4_mat_norm[row_s:row_e,:]
# matrix4 = fl1_mat_norm[row_s:row_e,:]
# matrix5 = fl2_mat_norm[row_s:row_e,:]
# matrix6 = fl4_mat_norm[row_s:row_e,:]
# matrix7 = adc1_mat_norm[row_s:row_e,:]
# matrix8 = adc2_mat_norm[row_s:row_e,:]
# matrix9 = adc4_mat_norm[row_s:row_e,:]

# Normalized matrices
matrix1 = ce1_matrix[row_s:row_e,:]
matrix2 = ce2_matrix[row_s:row_e,:]
matrix3 = ce4_matrix[row_s:row_e,:]
matrix4 = fl1_matrix[row_s:row_e,:]
matrix5 = fl2_matrix[row_s:row_e,:]
matrix6 = fl4_matrix[row_s:row_e,:]
matrix7 = adc1_matrix[row_s:row_e,:]
matrix8 = adc2_matrix[row_s:row_e,:]
matrix9 = adc4_matrix[row_s:row_e,:]

# Combine the matrices in a list

modality = 2

if modality == 1:
    matrices   = [matrix1, matrix2, matrix3]
    title_name = 'Box Plot for Each Row Across Matrices for T1ce'   
elif modality == 2:
    matrices = [matrix4, matrix5, matrix6]
    title_name = 'Box Plot for Each Row Across Matrices for T2flair'
elif modality == 3:
    matrices = [matrix7, matrix8, matrix9]
    title_name = 'Box Plot for Each Row Across Matrices for ADC'

# Shape 0-13
# First 14-31
# GLCM 32-55
# GLDM 56-69
# GLRLM 70-85
# GLSZM 86-101
# NGTDM 102-106

# Assume X is your feature matrix and Y is your label vector
# Selecting the first 10 features
X = X
Y = y

features_to_plot = X.columns[:14]

# Merge X and Y into a single DataFrame for easier plotting
data = pd.concat([X[features_to_plot], Y], axis=1)
data.columns = list(features_to_plot) + ['MGMT']

# Set up subplots
plt.figure(figsize=(15, 10))
for i, feature in enumerate(features_to_plot, 1):
    plt.subplot(3, 5, i)  # Create a 2x5 grid of subplots for 10 features
    data.boxplot(column=feature, by='MGMT', ax=plt.gca())
    plt.title(feature)
    plt.suptitle('')  # Suppress the default title
    plt.xlabel('Label')
    plt.ylabel(feature)

plt.tight_layout()
plt.show()

#%% Use TPOT API
from tpot import TPOTClassifier

tpot_accuracy_scores = []
tpot_ruc_scores      = []
tpot_accuracy_train  = []
tpot_ruc_train       = []

for train_index, test_index in kf.split(X_numeric, y_numeric):
    
    X_train, X_test = X_numeric[train_index], X_numeric[test_index]
    y_train, y_test = y_numeric[train_index], y_numeric[test_index]
    
    # X_train, X_test = X.iloc[train_index], X.iloc[test_index]
    # y_train, y_test = y.iloc[train_index], y.iloc[test_index]
    X_resampled, y_resampled = smote.fit_resample(X_train, y_train)

    # Initialize TPOTClassifier
    tpot = TPOTClassifier(generations=5, population_size=20, verbosity=2, scoring='roc_auc')
    # tpot.fit(X_resampled, y_resampled)
    tpot.fit(X_train, y_train)
    print("Test Accuracy:", tpot.score(X_test, y_test))
    
    # Make predictions on the cleaned test set
    y_pred       = tpot.predict(X_test)
    y_pred_train = tpot.predict(X_train)
    
    # Calculate accuracy
    tpot_accuracy = accuracy_score(y_test, y_pred)
    tpot_ruc_acc  = roc_auc_score(y_test, y_pred)
    tpot_accuracy_scores.append(tpot_accuracy)
    tpot_ruc_scores.append(tpot_ruc_acc)
    
    tpot_accuracy_t = accuracy_score(y_train, y_pred_train)
    tpot_ruc_acc_t  = roc_auc_score(y_train, y_pred_train)
    tpot_accuracy_train.append(tpot_accuracy_t)
    tpot_ruc_train.append(tpot_ruc_acc_t)
  
# Calculate the average accuracy over all folds
tpot_average_accuracy = np.mean(tpot_accuracy_scores)
tpot_average_accuracy_train = np.mean(tpot_accuracy_train)
print(f"Average Test accuracy over {k} folds: {tpot_average_accuracy:.4f}")
print(f"Average Train accuracy over {k} folds: {tpot_average_accuracy_train:.4f}")

tpot_average_ruc = np.mean(tpot_ruc_scores)
tpot_average_ruc_train = np.mean(tpot_ruc_train)
print(f"Average Test RUC over {k} folds: {tpot_average_ruc:.4f}")
print(f"Average Train RUC over {k} folds: {tpot_average_ruc_train:.4f}")

#%% ANOVA F-test
from sklearn.feature_selection import f_classif
import numpy as np

# Assuming X_train is a feature matrix and Y_train is the target (class labels)
f_values, p_values = f_classif(X_train, y_train)

# Print the F-statistics and corresponding p-values for each feature
for i, (f_val, p_val) in enumerate(zip(f_values, p_values)):
    print(f"Feature {i+1}: F-statistic = {f_val}, P-value = {p_val}")

#%% Mann-Whitney Test

from scipy.stats import mannwhitneyu
import numpy as np

# Example setup: `y_pred` is your model's predictions, and `y_test` is the true labels.
y_pred = model.predict(X_test)  # Your model's predictions on test data
y_test = np.array(y_test)  # Ensure `y_test` is an array for easier handling

# Separate predictions into two groups based on true labels
group1 = y_pred[y_test == 0]  # Predictions where true label is 0
group2 = y_pred[y_test == 1]  # Predictions where true label is 1

# group1 = y_test[y_test == 0]  # Predictions where true label is 0
# group2 = y_test[y_test == 1]  # Predictions where true label is 1

# Perform the Mann-Whitney U test
stat, p_value = mannwhitneyu(group1, group2, alternative='two-sided')

print(f"Mann-Whitney U Test Statistic: {stat}")
print(f"P-value: {p_value}")

# Interpret the results
if p_value < 0.05:
    print("Significant difference between groups, suggesting predictive power.")
else:
    print("No significant difference between groups, indicating limited predictive power.")

#%% BH selection

import numpy as np
from scipy.stats import ttest_ind

X_train = X_train

# Assuming X_train is your feature matrix and Y_train is your target variable
p_values = []
for i in range(X_train.shape[1]):
    feature = X_train[:, i]
    group1 = feature[y_train == 0]
    group2 = feature[y_train == 1]
    _, p_value = ttest_ind(group1, group2, equal_var=False)  # p-value for the feature
    p_values.append(p_value)

print("P-values for each feature:", p_values)

from statsmodels.stats.multitest import multipletests
# Apply the Benjamini-Hochberg correction
alpha = 0.05  # Desired false discovery rate
_, corrected_p_values, _, _ = multipletests(p_values, alpha=alpha, method='fdr_bh')

print("Corrected p-values:", corrected_p_values)

#%% Plot the scatters

# Replace 'feature_name' with the name of the feature you want to plot

# Create a scatter plot with different colors for each class in y
plt.figure(figsize=(8, 6))
col_id = 3
# gene_filter_1['MGMT']
# col_name = 'original_shape_MeshVolume'
col_name = 'original_shape_SurfaceArea'
# sns.scatterplot(x=X.iloc[:,col_id], y=gene_filter_1['MGMT'], hue=y, palette={0: "blue", 1: "red"}, s=50)
plt.scatter(X[col_name], y=y)

# Add labels and title
plt.xlabel(col_name)
plt.ylabel("Survival_Label")
plt.title(f'Scatter Plot of {X.columns[col_id]} by Class')
# plt.legend(title="Class", labels=["Class 0", "Class 1"])

plt.show()

#%% Check outlier test

from scipy.stats import zscore
import numpy as np

# Assuming `X` is your feature matrix (numpy array or pandas DataFrame)
z_scores = np.abs(zscore(X))

# Set a threshold, e.g., Z-score of 3
threshold = 6

# Find the indices where the Z-score is greater than the threshold
outlier_indices = np.where(z_scores > threshold)

# Retrieve row and column names for outliers
outlier_rows = X.index[outlier_indices[0]]
outlier_columns = X.columns[outlier_indices[1]]

# Print the outlier names
# for row, col in zip(outlier_rows, outlier_columns):
    # print(f"Outlier in '{col}' for '{row}' with Z-score: {z_scores[outlier_indices[0][0], outlier_indices[1][0]]:.2f}")

import pandas as pd
import numpy as np
from scipy.stats import zscore

# Create a list to store columns to remove
columns_to_remove = []

# Identify columns with outliers
for col in X.columns:
    # Check if there are any outliers in the column
    if (np.abs(z_scores[col]) > threshold).any():
        columns_to_remove.append(col)

# Remove the identified columns from the original DataFrame
X_cleaned = X.drop(columns=columns_to_remove)

# Print the names of the columns removed
print(f"Removed columns: {columns_to_remove}")

# Display the remaining feature matrix
print("Remaining feature matrix:")
print(X_cleaned)

X_numeric = X_cleaned.to_numpy()  # or X.values
y_numeric = y.to_numpy()   # or Y.values

# Define the number of folds for k-fold cross-validation
k = 5

# Initialize KFold
kf = KFold(n_splits=k, shuffle=True, random_state=42)

# Split the data into training and testing sets
# X_train, X_test, y_train, y_test = train_test_split(X_numeric, y_numeric, test_size=0.2, random_state=42)

########### Test models
# Store accuracy scores for each fold
accuracy_scores = []
accuracy_train  = []

# Perform k-fold cross-validation
for train_index, test_index in kf.split(X):
    X_train, X_test = X.iloc[train_index], X.iloc[test_index]
    y_train, y_test = y.iloc[train_index], y.iloc[test_index]

    # Calculate Z-scores for each column in the training set
    z_scores = X_train.apply(zscore)

    # Define a threshold for identifying outliers (e.g., Z-score > 3)
    threshold = 100

    # Identify columns to remove
    columns_to_remove = []
    for col in X_train.columns:
        if (np.abs(z_scores[col]) > threshold).any():
            columns_to_remove.append(col)

    # Remove identified columns from training and test sets
    X_train_cleaned = X_train.drop(columns=columns_to_remove)
    X_test_cleaned  = X_test.drop(columns=columns_to_remove)

    # Train the model using the cleaned training set
    model = LogisticRegression(max_iter=1000,class_weight='balanced')  # Use any model you prefer
    model.fit(X_train_cleaned, y_train)

    # Make predictions on the cleaned test set
    y_pred = model.predict(X_test_cleaned)
    y_pred_train = model.predict(X_train_cleaned)

    # Calculate accuracy
    accuracy = accuracy_score(y_test, y_pred)
    accuracy_scores.append(accuracy)
    
    accuracy_t = accuracy_score(y_train, y_pred_train)
    accuracy_train.append(accuracy_t)

# Calculate the average accuracy over all folds
average_accuracy = np.mean(accuracy_scores)
average_accuracy_train = np.mean(accuracy_train)
print(f"Average Test accuracy over {k} folds: {average_accuracy:.4f}")
print(f"Average Train accuracy over {k} folds: {average_accuracy_train:.4f}")
# Define the model
# rf = RandomForestClassifier()

# # Define the hyperparameters to tune
# param_grid = {
#     'n_estimators': [50, 100, 200],
#     'max_depth': [5, 10, 15],
#     'min_samples_split': [2, 5, 10],
#     'min_samples_leaf': [1, 2, 4],
#     'max_features': ['sqrt', 'log2', None]
# }

# Grid search
# grid_search = GridSearchCV(estimator=rf, param_grid=param_grid, cv=5, scoring='roc_auc')
# grid_search.fit(X_train, y_train)

# Best parameters
# print("Best parameters:", grid_search.best_params_)

########## Evaluate model

# Choose and train the model
# best_params = grid_search.best_params_
# model       = RandomForestClassifier(**best_params)
# model = LogisticRegression(max_iter=1000, class_weight='balanced')  # Increase max_iter if convergence issues occur
# # model = GradientBoostingClassifier(n_estimators=100, learning_rate=0.1, max_depth=3)
# # model.fit(X_train, y_train)

# # Make predictions
# y_pred = model.predict(X_test)
# y_pred_proba  = model.predict_proba(X_test)[:, 1]  # Get probabilities for ROC AUC
# y_train_pred  = model.predict(X_train)
# y_train_proba = model.predict_proba(X_train)[:,1]  # 

# y_train_round = np.round(y_train_proba)
# y_pred_round = np.round(y_pred_proba)

# # Evaluate the model
# print(classification_report(y_test, y_pred))
# print("Test ROC AUC Score:",  roc_auc_score(y_test, y_pred_round))
# print("Train ROC AUC Score:", roc_auc_score(y_train, y_train_round))

# # Confusion matrix
# conf_matrix = confusion_matrix(y_test, y_pred)
# print("Confusion Matrix:\n", conf_matrix)




#%% Check the resampled results

# from sklearn.model_selection import train_test_split, GridSearchCV, KFold, cross_val_score
# from sklearn.ensemble import RandomForestClassifier
# from sklearn.linear_model import LogisticRegression
# from sklearn.metrics import classification_report, roc_auc_score, confusion_matrix
# import pandas as pd

# # Assuming final_df and df_gene are already loaded

# # Merge dataframes and separate features/target
# merged_df = pd.merge(final_df, df_gene[['CDKN2A/B.1']], left_index=True, right_index=True)

# X = merged_df.drop(['CDKN2A/B.1'], axis=1)  # Radiomic features
# y = merged_df['CDKN2A/B.1']  # CDKN2A/B mutation status

# # Define the model (Random Forest for grid search)
# rf = RandomForestClassifier()

# # Define the hyperparameters to tune
# param_grid = {
#     'n_estimators': [50, 100, 200],
#     'max_depth': [5, 10, 15],
#     'min_samples_split': [2, 5, 10],
#     'min_samples_leaf': [1, 2, 4],
#     'max_features': ['sqrt', 'log2', None]
# }

# # K-fold cross-validation (5-fold in this case)
# kf = KFold(n_splits=5, shuffle=True, random_state=42)

# # Grid search with cross-validation
# grid_search = GridSearchCV(estimator=rf, param_grid=param_grid, cv=kf, scoring='roc_auc')
# grid_search.fit(X, y)  # Use the full dataset with cross-validation

# # Best parameters
# print("Best parameters:", grid_search.best_params_)

# ########## Evaluate model with K-Fold Cross-Validation ##########
# # Use Logistic Regression as the chosen model for final evaluation
# model = LogisticRegression(max_iter=1000, class_weight='balanced')

# # Perform cross-validation and get the average score across all folds
# cv_scores = cross_val_score(model, X, y, cv=kf, scoring='roc_auc')

# # Print the cross-validation scores and the mean ROC AUC score
# print("Cross-validation ROC AUC scores:", cv_scores)
# print("Mean ROC AUC score:", cv_scores.mean())

# ########## Final Evaluation on Training/Test Split (Optional) ##########

# # Split the data into training and testing sets for a final evaluation
# X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

# # Train and evaluate the model on the split data
# model.fit(X_train, y_train)

# y_pred = model.predict(X_test)
# y_pred_proba = model.predict_proba(X_test)[:, 1]  # Get probabilities for ROC AUC
# y_train_proba = model.predict_proba(X_train)[:,1]  # 

# # Evaluate the model
# print(classification_report(y_test, y_pred))
# print("Test ROC AUC Score:", roc_auc_score(y_test, y_pred_proba))
# print("Train ROC AUC Score:", roc_auc_score(y_train, y_train_proba))

# # Confusion matrix
# conf_matrix = confusion_matrix(y_test, y_pred)
# print("Confusion Matrix:\n", conf_matrix)

#%% Test for K-Fold

import pandas as pd
import numpy as np
from scipy.stats import zscore
from sklearn.model_selection import KFold
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score

# Assuming X is your DataFrame with patients as rows and features as columns
# and y is your target variable

# Define the number of folds for k-fold cross-validation
k = 5

# Initialize KFold
kf = KFold(n_splits=k, shuffle=True, random_state=42)

# Store accuracy scores for each fold
accuracy_scores = []
ruc_scores      = []
accuracy_train  = []
ruc_train       = []

# Perform k-fold cross-validation
for train_index, test_index in kf.split(X):
    X_train, X_test = X.iloc[train_index], X.iloc[test_index]
    y_train, y_test = y.iloc[train_index], y.iloc[test_index]

    # Calculate Z-scores for each column in the training set
    z_scores = X_train.apply(zscore)

    # Define a threshold for identifying outliers (e.g., Z-score > 3)
    threshold = 3

    # Identify columns to remove
    columns_to_remove = []
    for col in X_train.columns:
        if (np.abs(z_scores[col]) > threshold).any():
            columns_to_remove.append(col)

    # Remove identified columns from training and test sets
    X_train_cleaned = X_train.drop(columns=columns_to_remove)
    X_test_cleaned = X_test.drop(columns=columns_to_remove)

    # Train the model using the cleaned training set
    model = LogisticRegression(max_iter=1000)  # Use any model you prefer
    model.fit(X_train_cleaned, y_train)

    # Make predictions on the cleaned test set
    y_pred = model.predict(X_test_cleaned)
    y_pred_train = model.predict(X_train_cleaned)

    # Calculate accuracy
    accuracy = accuracy_score(y_test, y_pred)
    ruc_acc  = roc_auc_score(y_test, y_pred)
    accuracy_scores.append(accuracy)
    ruc_scores.append(ruc_acc)
    
    accuracy_t = accuracy_score(y_train, y_pred_train)
    ruc_acc_t  = roc_auc_score(y_train, y_pred_train)
    accuracy_train.append(accuracy_t)
    ruc_train.append(ruc_acc_t)

# Calculate the average accuracy over all folds
average_accuracy = np.mean(accuracy_scores)
average_accuracy_train = np.mean(accuracy_train)
print(f"Average Test accuracy over {k} folds: {average_accuracy:.4f}")
print(f"Average Train accuracy over {k} folds: {average_accuracy_train:.4f}")

average_ruc = np.mean(ruc_scores)
average_ruc_train = np.mean(ruc_train)
print(f"Average Test RUC over {k} folds: {average_ruc:.4f}")
print(f"Average Train RUC over {k} folds: {average_ruc_train:.4f}")
    
#%% Test performance for SVM

from sklearn.model_selection import KFold
from sklearn.svm import SVC
from sklearn.metrics import accuracy_score, roc_auc_score
import numpy as np
import pandas as pd

# Merge and prepare data
# merged_df = pd.merge(final_df, df_gene[['CDKN2A/B.1']], left_index=True, right_index=True)

# Separate features and target variable
X = merged_df.drop(['MGMT'], axis=1)  # Radiomic features
y = merged_df['MGMT']  #== CDKN2A/B mutation status

X_numeric = X.to_numpy()  # or X.values
y_numeric = y.to_numpy()   # or Y.values

# Define the number of folds for k-fold cross-validation
k = 5
kf = KFold(n_splits=k, shuffle=True, random_state=42)

# Store accuracy and ROC-AUC scores for each fold
accuracy_scores = []
roc_auc_scores = []
accuracy_train = []
roc_auc_train = []

# Perform k-fold cross-validation with SVM
for train_index, test_index in kf.split(X_numeric):
    X_train, X_test = X_numeric[train_index], X.iloc[test_index]
    y_train, y_test = y.iloc[train_index], y.iloc[test_index]
    
    # Define and train the SVM model
    model = SVC(kernel='rbf', probability=True)  # Use kernel='rbf' or 'linear' for nonlinear
    model.fit(X_train, y_train)
    
    # Make predictions
    y_pred = model.predict(X_test)
    y_pred_train = model.predict(X_train)
    
    # Get probabilities for ROC-AUC calculation
    y_pred_proba = model.predict_proba(X_test)[:, 1]
    y_pred_proba_train = model.predict_proba(X_train)[:, 1]
    
    # Calculate and store accuracy and ROC-AUC scores
    accuracy_scores.append(accuracy_score(y_test, y_pred))
    roc_auc_scores.append(roc_auc_score(y_test, y_pred_proba))
    accuracy_train.append(accuracy_score(y_train, y_pred_train))
    roc_auc_train.append(roc_auc_score(y_train, y_pred_proba_train))

# Calculate the average scores over all folds
average_accuracy = np.mean(accuracy_scores)
average_accuracy_train = np.mean(accuracy_train)
average_roc_auc = np.mean(roc_auc_scores)
average_roc_auc_train = np.mean(roc_auc_train)

print(f"Average Test Accuracy over {k} folds: {average_accuracy:.4f}")
print(f"Average Train Accuracy over {k} folds: {average_accuracy_train:.4f}")
print(f"Average Test ROC-AUC over {k} folds: {average_roc_auc:.4f}")
print(f"Average Train ROC-AUC over {k} folds: {average_roc_auc_train:.4f}")

#%% Try KNN

from sklearn.model_selection import KFold
from sklearn.neighbors import KNeighborsClassifier
from sklearn.metrics import accuracy_score, roc_auc_score
import numpy as np
import pandas as pd

# Merge and prepare data
# merged_df = pd.merge(final_df, df_gene[['CDKN2A/B.1']], left_index=True, right_index=True)

# Separate features and target variable
X = merged_df.drop(['MGMT'], axis=1)  # Radiomic features
y = merged_df['MGMT']  #== CDKN2A/B mutation status
# Define the number of folds for k-fold cross-validation
k = 5
kf = KFold(n_splits=k, shuffle=True, random_state=42)

# Store accuracy and ROC-AUC scores for each fold
accuracy_scores = []
roc_auc_scores = []
accuracy_train = []
roc_auc_train = []

# Initialize the KNN model (with example of 5 neighbors)
model = KNeighborsClassifier(n_neighbors=5)

# Perform k-fold cross-validation with KNN
for train_index, test_index in kf.split(X):
    X_train, X_test = X.iloc[train_index], X.iloc[test_index]
    y_train, y_test = y.iloc[train_index], y.iloc[test_index]
    
    # Train the KNN model
    model.fit(X_train, y_train)
    
    # Make predictions
    y_pred = model.predict(X_test)
    y_pred_train = model.predict(X_train)
    
    # Get probabilities for ROC-AUC calculation
    y_pred_proba = model.predict_proba(X_test)[:, 1]
    y_pred_proba_train = model.predict_proba(X_train)[:, 1]
    
    # Calculate and store accuracy and ROC-AUC scores
    accuracy_scores.append(accuracy_score(y_test, y_pred))
    roc_auc_scores.append(roc_auc_score(y_test, y_pred_proba))
    accuracy_train.append(accuracy_score(y_train, y_pred_train))
    roc_auc_train.append(roc_auc_score(y_train, y_pred_proba_train))

# Calculate the average scores over all folds
average_accuracy = np.mean(accuracy_scores)
average_accuracy_train = np.mean(accuracy_train)
average_roc_auc = np.mean(roc_auc_scores)
average_roc_auc_train = np.mean(roc_auc_train)

print(f"Average Test Accuracy over {k} folds: {average_accuracy:.4f}")
print(f"Average Train Accuracy over {k} folds: {average_accuracy_train:.4f}")
print(f"Average Test ROC-AUC over {k} folds: {average_roc_auc:.4f}")
print(f"Average Train ROC-AUC over {k} folds: {average_roc_auc_train:.4f}")

#%% Next Tries

import numpy as np
import pandas as pd
from pingouin import intraclass_corr  # Install with `pip install pingouin`

# Assuming `X` is your feature DataFrame with repeated measurements
def calculate_icc(df):
    """
    Calculate ICC for each feature across repeated measures.
    Assume data format is suitable for ICC calculation, or reformat if needed.
    """
    icc_values = {}
    for col in df.columns:
        data = pd.DataFrame({'rater': list(range(df.shape[1])), 'subject': list(range(df.shape[0])), 'score': df[col]})
        icc = intraclass_corr(data=data, targets='subject', raters='rater', ratings='score')
        icc_values[col] = icc['ICC'][0]  # Use ICC value directly

    # Filter features with ICC >= 0.75
    repeatable_features = [feature for feature, icc_value in icc_values.items() if icc_value >= 0.75]
    return repeatable_features

# repeatable_features = calculate_icc(X)
# X_repeatable = X[repeatable_features]  # Keep only repeatable features

X_repeatable = X  # Keep only repeatable features

# Calculate Pearson correlation matrix
corr_matrix = X_repeatable.corr().abs()

# Select upper triangle of correlation matrix
upper_triangle = corr_matrix.where(np.triu(np.ones(corr_matrix.shape), k=1).astype(bool))

# Find features with correlation greater than 0.9
high_corr_features = [column for column in upper_triangle.columns if any(upper_triangle[column] > 0.9)]

# Drop highly correlated features
X_uncorrelated = X_repeatable.drop(high_corr_features, axis=1)

import pymrmr  # Install with `pip install pymrmr`

# Convert DataFrame to a format suitable for mRMR
# Make sure the label column is included for supervised selection
y = y # Replace with actual label data
X_with_label = pd.concat([y, X_uncorrelated], axis=1)

# Apply mRMR and select top 20% features
num_features_to_select = int(0.2 * X_uncorrelated.shape[1])
selected_features = pymrmr.mRMR(X_with_label, 'MIQ', num_features_to_select)

# Filter to the selected features
X_selected = X_uncorrelated[selected_features]












