#%% This code focuses on analyzing the prediction accuracy for PDGM dataset

import numpy as np
import matplotlib.pyplot as plt

import os
import pandas as pd

import nibabel as nib
#import radiomics
from radiomics import featureextractor

#%% 

multi  = 8
scalar = 2

def generate_uniform_sphere_points(n_points=8000):
    """
    Generate uniformly distributed points on a unit sphere using the Fibonacci spiral method.
    
    Parameters:
    n_points (int): Number of points to generate
    
    Returns:
    tuple: Arrays of phi (polar) and theta (azimuthal) angles
    """
    indices = np.arange(0, n_points, dtype=float) + 0.5
    
    phi   = np.arccos(1 - 2*indices/n_points)  # polar angle
    theta = np.pi * (1 + 5**0.5) * indices   # azimuthal angle
    
    # Normalize theta to [-pi, pi]
    theta = np.mod(theta + np.pi, 2*np.pi) - np.pi
    
    return phi, theta

###############################################################################
def spherical_to_cartesian(r, phi, theta):
    """
    Convert spherical coordinates to Cartesian coordinates.
    
    Parameters:
    r (float): Radius
    phi (array): Polar angles [0, pi]
    theta (array): Azimuthal angles [-pi, pi]
    
    Returns:
    tuple: x, y, z coordinates
    """
    x = r * np.sin(phi) * np.cos(theta)
    y = r * np.sin(phi) * np.sin(theta)
    z = r * np.cos(phi)
    return x, y, z

###############################################################################
def sample_sphere_shell(image, center, radius, n_points=1000):
    """
    Sample image values on a spherical shell using uniform spherical coordinates.
    
    Parameters:
    image (np.ndarray): 3D image data
    center (tuple): (x,y,z) coordinates of sphere center
    radius (float): Radius of the shell
    n_points (int): Number of points to sample on the sphere
    
    Returns:
    tuple: phi angles, theta angles, and sampled values
    """
    # Generate uniform points on unit sphere
    phi, theta = generate_uniform_sphere_points(n_points)
    
    # Convert to Cartesian coordinates with desired radius
    x, y, z = spherical_to_cartesian(radius, phi, theta)
    
    # Shift by center coordinates
    x += center[0]
    y += center[1]
    z += center[2]
    
    # Round to nearest integer for coordinates
    x_idx = np.clip(np.round(x).astype(int), 0, image.shape[0]-1)
    y_idx = np.clip(np.round(y).astype(int), 0, image.shape[1]-1)
    z_idx = np.clip(np.round(z).astype(int), 0, image.shape[2]-1)
    
    # Sample image values
    values = image[x_idx, y_idx, z_idx]
    
    return phi, theta, values

###############################################################################
def create_shell_projection(phi, theta, values, n_phi=90*multi, n_theta=180*multi):
    """
    Create a 2D projection map from spherical samples.
    
    Parameters:
    phi (array): Polar angles
    theta (array): Azimuthal angles
    values (array): Sampled intensity values
    n_phi (int): Number of bins for polar angle
    n_theta (int): Number of bins for azimuthal angle
    
    Returns:
    tuple: 2D intensity map and weight map
    """
    # Create 2D histogram
    intensity_map = np.zeros((n_phi, n_theta))
    weight_map    = np.zeros((n_phi, n_theta))
    
    # Convert angles to bin indices
    phi_bins   = np.linspace(0, np.pi, n_phi+1)
    theta_bins = np.linspace(-np.pi, np.pi, n_theta+1)
    
    phi_idx   = np.digitize(phi, phi_bins) - 1
    theta_idx = np.digitize(theta, theta_bins) - 1
    
    # Accumulate values
    for i in range(len(values)):
        intensity_map[phi_idx[i], theta_idx[i]] += values[i]
        weight_map[phi_idx[i], theta_idx[i]] += 1
    
    # Normalize
    mask = weight_map > 0
    intensity_map[mask] /= weight_map[mask]
    
    return intensity_map, weight_map

###############################################################################
def analyze_mask_region(mask, center):
    """
    Analyze mask region to determine radial boundaries.
    
    Parameters:
    mask (np.ndarray): Binary mask of the region
    center (tuple): (x,y,z) coordinates of center point
    
    Returns:
    tuple: min_radius, max_radius, and distance map
    """
    shape = mask.shape
    x, y, z = np.meshgrid(np.arange(shape[0]), 
                         np.arange(shape[1]), 
                         np.arange(shape[2]), 
                         indexing="ij")
    
    # Calculate distance map
    r = np.sqrt((x - center[0])**2 + 
                (y - center[1])**2 + 
                (z - center[2])**2)
    
    valid_mask = mask > 0
    valid_r = r[valid_mask]
    
    min_radius = np.min(valid_r)
    max_radius = np.max(valid_r)
    
    return min_radius, max_radius, r

###############################################################################
def generate_shells_from_mask(image, mask, center, num_shells=50, points_per_shell=8000):
    """
    Generate spherical shells based on mask region boundaries.
    
    Parameters:
    image (np.ndarray): 3D image data
    mask (np.ndarray): Binary mask defining the region
    center (tuple): (x,y,z) coordinates of center point
    num_shells (int): Number of shells to generate
    points_per_shell (int): Number of points to sample per shell
    
    Returns:
    list: List of (intensity_map, weight_map) tuples for each shell
    np.ndarray: Array of shell radii
    """
    # Get mask region boundaries
    min_radius, max_radius, r = analyze_mask_region(mask, center)
    # print(f"Region boundaries - Min radius: {min_radius:.2f}, Max radius: {max_radius:.2f}")
    
    # Generate shell radii using cubic root spacing for uniform volume distribution
    #alpha = np.linspace(0, 1, num_shells) ** (1/3)
    alpha = np.linspace(0, 1, num_shells)
    shell_radii = min_radius + (max_radius - min_radius) * alpha
    
    # Generate uniform points on unit sphere
    phi, theta = generate_uniform_sphere_points(points_per_shell)
    
    shell_results = []
    for i, radius in enumerate(shell_radii):
        # print(f"Processing shell {i+1}/{num_shells} at radius {radius:.2f}")
        
        # Convert to Cartesian coordinates for this radius
        x, y, z = spherical_to_cartesian(radius, phi, theta)
        
        # Shift by center coordinates
        x += center[0]
        y += center[1]
        z += center[2]
        # print(f'The values for x is {x}')
        
        # Round to nearest integer for coordinates
        x_idx = np.clip(np.round(x).astype(int), 0, image.shape[0]-1)
        y_idx = np.clip(np.round(y).astype(int), 0, image.shape[1]-1)
        z_idx = np.clip(np.round(z).astype(int), 0, image.shape[2]-1)
        # print(f'The indices are {x_idx, y_idx, z_idx}')
        
        # Only sample points that fall within the mask
        valid_points = mask[x_idx, y_idx, z_idx] >= 0
        # print(valid_points)
        
        if np.sum(valid_points) > 0:
            # Sample image values for valid points
            values = image[x_idx[valid_points], 
                         y_idx[valid_points], 
                         z_idx[valid_points]]
            
            # Create 2D projection for this shell
            intensity_map, weight_map = create_shell_projection(
                phi[valid_points], 
                theta[valid_points], 
                values
            )
            
            shell_results.append((intensity_map, weight_map))
        else:
            print(f"Warning: No valid points found in shell {i+1}")
            shell_results.append((None, None))
    
    return shell_results, shell_radii

###############################################################################
def visualize_shell_sequence(shell_results, shell_radii, region_name):
    """
    Visualize a sequence of shell projections.
    
    Parameters:
    shell_results: List of (intensity_map, weight_map) tuples
    shell_radii: Array of shell radii
    region_name: Name of the region being analyzed
    """
    valid_shells = []
    for i, (radius, result) in enumerate(zip(shell_radii, shell_results)):
        if result is not None and result[0] is not None:
            valid_shells.append((i, radius, result[0], result[1]))
    
    if not valid_shells:
        print("No valid shells to visualize")
        return
    
    n_shells = len(valid_shells)
    n_cols   = min(6, n_shells)
    n_rows   = (n_shells + n_cols - 1) // n_cols
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(4*n_cols, 4*n_rows))
    if n_rows == 1:
        axes = axes.reshape(1, -1)
    
    for idx, (i, radius, intensity_map, weight_map) in enumerate(valid_shells):
        row = idx // n_cols
        col = idx % n_cols
        
        im = axes[row, col].imshow(intensity_map, 
                                 extent=[-np.pi, np.pi, 0, np.pi],
                                 aspect='auto',
                                 cmap='viridis')
        axes[row, col].set_title(f'Shell {i}\nRadius: {radius:.2f}')
        plt.colorbar(im, ax=axes[row, col])
    
    # Clear unused subplots
    for idx in range(len(valid_shells), n_rows * n_cols):
        row = idx // n_cols
        col = idx %  n_cols
        axes[row, col].axis('off')
    
    plt.suptitle(f'Shell Projections for {region_name}')
    plt.tight_layout()
    plt.show()

###############################################################################
def get_tumor_center(mask):
    """
    Compute the center of mass of the tumor mask.
    Returns center coordinates in the original image space.
    """
    if not np.any(mask):
        raise ValueError("Empty mask")
        
    # Get coordinates of all non-zero voxels
    coords = np.array(np.where(mask>0)).T
    # Compute center of mass
    center = np.mean(coords, axis=0)
    
    return center

###############################################################################¥¥¥¥¥¥¥¥¥¥¥¥¥¥¥¥¥¥¥¥¥¥¥¥¥¥¥¥¥
# Function to process a specific region
def process_region(image, mask, region_label, region_name, num_shell):
    """
    Process a specific region in the mask with 3D visualization.
    """
    # Create binary mask for this region
    region_mask = (mask == region_label).astype(np.uint8)
    
    if not np.any(region_mask):
        print(f"No voxels found for region {region_name}")
        return
    
    # Get region center
    tumor_center = get_tumor_center(region_mask)
    
    # First, create single shell visualization
    min_radius, max_radius, _ = analyze_mask_region(region_mask, tumor_center)
    middle_radius = (min_radius + max_radius) / 2
    
    # Visualize single shell
    # print('Working on Singles')
    # visualize_shell_3d(image, tumor_center, middle_radius, n_points=8000)
    # plt.show()
    
    # Visualize multiple shells
    # print('Working on Multiples')
    # visualize_multiple_shells_3d(image, region_mask, tumor_center, 
                               # num_shells=50, points_per_shell=10000)
    # plt.show()
    
    # Generate and visualize original 2D projections
    shell_results, shell_radii = generate_shells_from_mask(
        image, 
        region_mask, 
        tumor_center,
        num_shells=num_shell,
        points_per_shell=scalar * 10000 * multi**2
    )
    
    if shell_results is not None:
    #     # shell_results = zoom(shell_results, zoom=2, order=3)
    #     # shell_results = gaussian_filter(shell_results, sigma=2)
        
        # Visualize 2D projection results
        # visualize_shell_sequence(shell_results, shell_radii, region_name)
        
        return shell_results
    
    else:
        print('No valid mask')
        return []
    
#%%
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def visualize_shell_3d(image, center, radius, n_points=1000, cmap='viridis', 
                      alpha=0.8, point_size=50, elev=30, azim=45):
    """
    Create an interactive 3D visualization of the spherical shell with intensity values.
    
    Parameters:
    image (np.ndarray): 3D image data
    center (tuple): (x,y,z) coordinates of sphere center
    radius (float): Radius of the shell
    n_points (int): Number of points to sample on the sphere
    cmap (str): Colormap for intensity values
    alpha (float): Transparency of points
    point_size (float): Size of scatter points
    elev (float): Initial elevation viewing angle
    azim (float): Initial azimuthal viewing angle
    
    Returns:
    tuple: Figure and Axes objects for further customization
    """
    # Generate uniform points on unit sphere
    phi, theta = generate_uniform_sphere_points(n_points)
    
    # Convert to Cartesian coordinates with desired radius
    x, y, z = spherical_to_cartesian(radius, phi, theta)
    
    # Shift by center coordinates
    x += center[0]
    y += center[1]
    z += center[2]
    
    # Round to nearest integer for coordinates
    x_idx = np.clip(np.round(x).astype(int), 0, image.shape[0]-1)
    y_idx = np.clip(np.round(y).astype(int), 0, image.shape[1]-1)
    z_idx = np.clip(np.round(z).astype(int), 0, image.shape[2]-1)
    
    # Sample image values
    values = image[x_idx, y_idx, z_idx]
    print(f'Max value is: {np.max(values)}')
    
    # Create figure
    fig = plt.figure(figsize=(12, 12))
    ax = fig.add_subplot(111, projection='3d')
    
    # Create scatter plot
    scatter = ax.scatter(x, y, z, c=values, cmap=cmap, alpha=alpha, s=point_size)
    
    # Add colorbar
    plt.colorbar(scatter, ax=ax, label='Intensity')
    
    # Set initial view
    ax.view_init(elev=elev, azim=azim)
    
    # Set labels and title
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title(f'3D Shell Visualization\nRadius: {radius:.2f}')
    
    # Make the plot more cubic
    max_range = np.array([x.max()-x.min(), y.max()-y.min(), z.max()-z.min()]).max() / 2.0
    mid_x = (x.max()+x.min()) * 0.5
    mid_y = (y.max()+y.min()) * 0.5
    mid_z = (z.max()+z.min()) * 0.5
    ax.set_xlim(mid_x - max_range, mid_x + max_range)
    ax.set_ylim(mid_y - max_range, mid_y + max_range)
    ax.set_zlim(mid_z - max_range, mid_z + max_range)
    
    return fig, ax

###############################################################################
def visualize_multiple_shells_3d(image, mask, center, num_shells=50, 
                            points_per_shell=10000, min_radius=None, max_radius=None):
    """
    Visualize each shell in 3D individually with intensity values.
    
    Parameters:
    image (np.ndarray): 3D image data
    mask (np.ndarray): Binary mask defining the region
    center (tuple): (x,y,z) coordinates of center point
    num_shells (int): Number of shells to visualize
    points_per_shell (int): Number of points per shell
    min_radius (float): Optional minimum radius (computed from mask if None)
    max_radius (float): Optional maximum radius (computed from mask if None)
    """
    # Get mask region boundaries if not provided
    # region_mask = (mask == region_label).astype(np.uint8)
    
    if min_radius is None or max_radius is None:
        min_r, max_r, _ = analyze_mask_region(mask, center)
        min_radius = min_radius or min_r
        max_radius = max_radius or max_r
    
    # print(f'Current min radius is {min_radius}')
    # print(f'Current max radius is {max_radius}')
    
    # Generate shell radii using cubic root spacing
    alpha = np.linspace(0, 1, num_shells) ** (1/3)
    shell_radii = min_radius + (max_radius - min_radius) * alpha
    
    for i, radius in enumerate(shell_radii):
        # Generate points for this shell
        phi, theta = generate_uniform_sphere_points(points_per_shell)
        x, y, z    = spherical_to_cartesian(radius, phi, theta)
        
        # Shift by center coordinates
        x += center[0]
        y += center[1]
        z += center[2]
        
        # Round to nearest integer for coordinates
        x_idx = np.clip(np.round(x).astype(int), 0, image.shape[0]-1)
        y_idx = np.clip(np.round(y).astype(int), 0, image.shape[1]-1)
        z_idx = np.clip(np.round(z).astype(int), 0, image.shape[2]-1)
        
        # Sample image values
        values = image[x_idx, y_idx, z_idx]
        
        # Create a figure for this shell
        fig = plt.figure(figsize=(8, 8))
        ax = fig.add_subplot(111, projection='3d')
        
        # Create scatter plot for this shell
        scatter = ax.scatter(x, y, z, c=values, cmap='viridis', 
                              alpha=0.8, s=30)
        
        # Add colorbar
        plt.colorbar(scatter, ax=ax, label='Intensity')
        
        # Set labels and title
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        ax.set_title(f'Shell Visualization: R = {radius:.1f}')
        
        # Make the plot more cubic
        all_coords = np.vstack((x, y, z))
        max_range  = np.max(all_coords.max(axis=1) - all_coords.min(axis=1)) / 2.0
        mid_x = (x.max() + x.min()) * 0.5
        mid_y = (y.max() + y.min()) * 0.5
        mid_z = (z.max() + z.min()) * 0.5
        ax.set_xlim(mid_x - max_range, mid_x + max_range)
        ax.set_ylim(mid_y - max_range, mid_y + max_range)
        ax.set_zlim(mid_z - max_range, mid_z + max_range)
        
        # Display the figure
        # plt.tight_layout()
        plt.show()

#%% Load the file number IDs

import glob
import re

file_loc   = '/Volumes/ResearchPhysics/zz_Haotian_Feng/GBM/Pre-Op/NF1-MT/seg_result/'
output_loc = '/Users/hafeng/Documents/Postdoc_Research_Meetings/TCGA/UCSF_IMAGING/Shell_Analysis/NF1_MT/'
num_blocks = 8

number_list = []

for filename in os.listdir(file_loc):
    
    number_list.append(filename)
        
number_list.sort()

#%% Read the files
import SimpleITK as sitk
from scipy.ndimage import distance_transform_edt
output_folder = output_loc
num_shell = 20

extractor = featureextractor.RadiomicsFeatureExtractor(preCrop=True)
extractor.settings['label'] = 1  # Set the label for tumor regions
extractor.settings['resampledPixelSpacing'] = None
extractor.settings['force2D'] = True
# extractor.settings['normalize'] = True
# extractor.settings['normalizeScale'] = 100

# For demonstration, process only the first item
for idx in number_list:
    
    if int(idx) > 146:
    
        folder_name = idx
        item_path   = os.path.join(file_loc, folder_name)
        
        # (1) Get the Pred File
        pred_path = '/' + idx + '/' + idx + '_pred.nii.gz'
        nii_pred  = nib.load(file_loc + pred_path)
        # Get the image data as a numpy array
        pred_data = nii_pred.get_fdata()
        
        # Assuming `mask_data` is the 3D mask matrix with values 0, 1, 2, 3
        mask_region_1 = (pred_data == 1)  # Mask for region 1
        mask_region_2 = (pred_data == 2)  # Mask for region 2
        mask_region_3 = (pred_data == 3)  # Mask for region 3
        mask = pred_data
        
        mask_quality = False
    
        if np.count_nonzero(mask_region_1) > 10 and np.count_nonzero(mask_region_2) > 10 and np.count_nonzero(mask_region_3) > 10:
            mask_quality = True
            
        if mask_quality:
            
            ### (2) Get the Flair File
            fla_path = '/' + idx + '/' + idx + '_fla_bet.nii.gz'
            nii_fla  = nib.load(file_loc + fla_path)
            
            # Get the image data as a numpy array
            image_fl = nii_fla.get_fdata()
            
            # Print some information about the data
            print("Shape of the flair data:", image_fl.shape)
            
            ### (3) Get the CE File
            ce_path = '/' + idx + '/' + idx + '_t1ca_bet.nii.gz'
            nii_ce  = nib.load(file_loc + ce_path)
            
            # Get the image data as a numpy array
            image_ce = nii_ce.get_fdata()
            
            # Print some information about the data
            print("Shape of the T1ce data:", image_ce.shape)
            
            # Region labels in the tumor mask
            region_labels = {
                "necrotic": 1,
                "enhancing": 2,
                "lesion": 3
            }
    
            # Process each region
            for region_name, label in region_labels.items():
                print(f"\nProcessing patient {idx} On {region_name} region (label {label})...")
                print(f"\nWorking on CE Data")
                results_ce = process_region(image_ce, mask, label, region_name, num_shell)   # num_shells, intensity, weight
                print(f"\nWorking on FL Data")
                results_fl = process_region(image_fl, mask, label, region_name, num_shell)   # num_shells, intensity, weight
                
                for ii in range(num_shell):
                    try:
                        if results_ce[ii] is None or results_fl[ii] is None:
                            print(f"Skipping shell {ii} in {region_name} due to missing data.")
                            continue
        
                        ce_shell_intensity = results_ce[ii][0]
                        fl_shell_intensity = results_fl[ii][0]
                        
                        ce_mask = np.ones_like(ce_shell_intensity, dtype=np.uint8)
                        ce_mask[0, :]  = 0  # Top edge
                        ce_mask[-1, :] = 0  # Bottom edge
                        ce_mask[:, 0]  = 0  # Left edge
                        ce_mask[:, -1] = 0  # Right edge
                        fl_mask = np.ones_like(fl_shell_intensity, dtype=np.uint8)
                        fl_mask[0, :]  = 0  # Top edge
                        fl_mask[-1, :] = 0  # Bottom edge
                        fl_mask[:, 0]  = 0  # Left edge
                        fl_mask[:, -1] = 0  # Right edge
                    
                        ce_sitk   = sitk.GetImageFromArray(ce_shell_intensity)
                        fl_sitk   = sitk.GetImageFromArray(fl_shell_intensity)
                        mask_sitk = sitk.GetImageFromArray(ce_mask)
                        
                        ce_result  = extractor.execute(ce_sitk, mask_sitk)
                        fl_result  = extractor.execute(fl_sitk, mask_sitk)
                        
                        # Save results to DataFrame
                        ce_features_df = pd.DataFrame(list(ce_result.items()), columns=['Feature', 'Value'])
                        fl_features_df = pd.DataFrame(list(fl_result.items()), columns=['Feature', 'Value'])
                        
                        # Write results to Excel files
                        ce_file_name = f"{idx}_label_{label}_shell_{ii}_CE_features_multi8_norm.xlsx"
                        ce_excel_path = os.path.join(output_folder, ce_file_name)
                        ce_features_df.to_excel(ce_excel_path, index=False)
                        
                        fl_file_name = f"{idx}_label_{label}_shell_{ii}_FL_features_multi8_norm.xlsx"
                        fl_excel_path = os.path.join(output_folder, fl_file_name)
                        fl_features_df.to_excel(fl_excel_path, index=False)
                        
                        print(f"Features saved for item {idx}, label {label}, shell {ii}")
                
                    except (TypeError, ValueError) as e:
                        print(f"Skipping shell {ii} in {region_name} due to error: {e}")
                        continue
        
        # for ii in range(num_shell):
        #     ce_shell_intensity = results_ce[ii][0]
        #     fl_shell_intensity = results_fl[ii][0]
            
        #     # Create a uniform mask of ones
        #     ce_mask = np.ones_like(ce_shell_intensity, dtype=np.uint8)
        #     fl_mask = np.ones_like(fl_shell_intensity, dtype=np.uint8)
            
        #     # Define file paths
        #     ce_path   = os.path.join(output_folder, f"ce_shell_{ii}_image.nii.gz")
        #     fl_path   = os.path.join(output_folder, f"fl_shell_{ii}_image.nii.gz")
        #     mask_path_ce = os.path.join(output_folder, f"shell_{ii}_mask_ce.nii.gz")
        #     mask_path_fl = os.path.join(output_folder, f"shell_{ii}_mask_fl.nii.gz")
            
        #     # Save as .nii.gz using nibabel
        #     ce_nifti   = nib.Nifti1Image(ce_shell_intensity, affine=np.eye(4))
        #     fl_nifti   = nib.Nifti1Image(fl_shell_intensity, affine=np.eye(4))
        #     mask_ce_nifti = nib.Nifti1Image(ce_mask, affine=ce_nifti.affine, dtype=np.uint8)
        #     mask_fl_nifti = nib.Nifti1Image(fl_mask, affine=fl_nifti.affine, dtype=np.uint8)
            
        #     nib.save(ce_nifti,   ce_path)
        #     nib.save(fl_nifti,   fl_path)
        #     nib.save(mask_ce_nifti, mask_path_ce)
        #     nib.save(mask_fl_nifti, mask_path_fl)
        
        #     ce_result = extractor.execute(ce_path, mask_path_ce)
        #     fl_result = extractor.execute(fl_path, mask_path_fl)

        #     # Save results to DataFrame
        #     ce_features_df = pd.DataFrame(list(ce_result.items()), columns=['Feature', 'Value'])
        #     fl_features_df = pd.DataFrame(list(fl_result.items()), columns=['Feature', 'Value'])
            
        #     # Write results to Excel files
        #     ce_file_name = f"{item}_label_{label}_shell_{ii}_CE_features.xlsx"
        #     ce_excel_path = os.path.join(output_folder, ce_file_name)
        #     ce_features_df.to_excel(ce_excel_path, index=False)
            
        #     fl_file_name = f"{item}_label_{label}_shell_{ii}_FL_features.xlsx"
        #     fl_excel_path = os.path.join(output_folder, fl_file_name)
        #     fl_features_df.to_excel(fl_excel_path, index=False)
        
#%%
# fig = plt.figure(figsize=(8, 8))
# ax = fig.add_subplot(111, projection='3d')

# ax.voxels(image, facecolors='blue', edgecolors=None, alpha=0.1)
# ax.voxels(mask,  facecolors='red',  edgecolors=None, alpha=0.7)
# # Set axis labels
# ax.set_xlabel('X-axis')
# ax.set_ylabel('Y-axis')
# ax.set_zlabel('Z-axis')

# plt.show()

# %%
