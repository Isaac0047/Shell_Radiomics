# This code focuses on analyzing the prediction accuracy for PDGM dataset

from scipy.ndimage import distance_transform_edt
import SimpleITK as sitk
import re
import glob
from itertools import cycle
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt

import os
import pandas as pd

import nibabel as nib
import radiomics
from radiomics import featureextractor
import io
import matplotlib.cm as cm

# %%

multi = 10
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

    phi = np.arccos(1 - 2*indices/n_points)  # polar angle
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
    weight_map = np.zeros((n_phi, n_theta))

    # Convert angles to bin indices
    phi_bins = np.linspace(0, np.pi, n_phi+1)
    theta_bins = np.linspace(-np.pi, np.pi, n_theta+1)

    phi_idx = np.digitize(phi, phi_bins) - 1
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
    min_radius = max(min_radius, 2)
    # print(f"Region boundaries - Min radius: {min_radius:.2f}, Max radius: {max_radius:.2f}")

    # Generate shell radii using cubic root spacing for uniform volume distribution
    alpha = np.linspace(0, 1, num_shells) ** (1/3)
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

def visualize_shell_sequence(shell_results, shell_radii, region_name, v_min, v_max):
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
    n_cols = min(5, n_shells)
    n_rows = (n_shells + n_cols - 1) // n_cols

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(4*n_cols, 4*n_rows))
    if n_rows == 1:
        axes = axes.reshape(1, -1)

    for idx, (i, radius, intensity_map, weight_map) in enumerate(valid_shells):
        row = idx // n_cols
        col = idx % n_cols

        im = axes[row, col].imshow(intensity_map,
                                   extent=[-np.pi, np.pi, 0, np.pi],
                                   aspect='auto',
                                   cmap='jet',
                                   vmin = v_min,
                                   vmax = v_max)
        axes[row, col].set_title(f'Shell {i}\nRadius: {radius:.2f}', fontsize=16)
        # Larger tick labels
        axes[row, col].tick_params(axis='x', labelsize=14)
        axes[row, col].tick_params(axis='y', labelsize=14)
        cbar = plt.colorbar(im, ax=axes[row, col])
        cbar.ax.tick_params(labelsize=12)  # <-- control fontsize
# cbar.set_label("Intensity", fontsize=14)  # optional: add label

    # Clear unused subplots
    for idx in range(len(valid_shells), n_rows * n_cols):
        row = idx // n_cols
        col = idx % n_cols
        axes[row, col].axis('off')

    # plt.suptitle(f'Shell Projections for {region_name}')
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
    coords = np.array(np.where(mask == 1)).T
    # Compute center of mass
    center = np.mean(coords, axis=0)

    return center

# Function to process a specific region


def process_region(image, mask, region_label, region_name, num_shell, vmin, vmax):
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
    print(min_radius)
    print(max_radius)
    middle_radius = (min_radius + max_radius) / 2

    # Generate and visualize original 2D projections
    shell_results, shell_radii = generate_shells_from_mask(
        image,
        region_mask,
        tumor_center,
        num_shells=num_shell,
        points_per_shell=scalar * 10000 * multi**2
    )

    # print('Start to visualize the layered shells')
    plot_static_shells(image, mask, tumor_center, num_shells=6,
                        points_per_shell=5000, min_radius=np.max([min_radius,3]), max_radius=max_radius)
    
    # plot_surface_shells(image, mask, tumor_center, num_shells=1,
                        # points_per_shell=20000, min_radius=np.max([min_radius,3]), max_radius=max_radius)

    if shell_results is not None:
        # shell_results = zoom(shell_results, zoom=2, order=3)
        # shell_results = gaussian_filter(shell_results, sigma=2)

        # Visualize 2D projection results
        # visualize_shell_sequence(shell_results, shell_radii, region_name, vmin, vmax)
        # print(f'Start working on visualization')
        # sliding_window_glrlm(shell_results, shell_radii, window_size=15, feature='ShortRunEmphasis')
        # sliding_window_glcm(shell_results, shell_radii, window_size=15, feature='Contrast')

        return shell_results

    else:
        print('No valid mask')
        return []


# %% Generate 3D plot visualizations ###########################################
def set_axes_equal(ax):
    """
    Make axes of 3D plot have equal scale so that spheres appear as spheres.
    Expands the limits to a cube that encloses the current data limits.
    """
    x_limits = ax.get_xlim3d()
    y_limits = ax.get_ylim3d()
    z_limits = ax.get_zlim3d()

    x_range = abs(x_limits[1] - x_limits[0])
    y_range = abs(y_limits[1] - y_limits[0])
    z_range = abs(z_limits[1] - z_limits[0])

    # The plot bounding box is a cube with sides = max range
    max_range = max([x_range, y_range, z_range]) / 2.0

    x_middle = np.mean(x_limits)
    y_middle = np.mean(y_limits)
    z_middle = np.mean(z_limits)

    ax.set_xlim3d([x_middle - max_range, x_middle + max_range])
    ax.set_ylim3d([y_middle - max_range, y_middle + max_range])
    ax.set_zlim3d([z_middle - max_range, z_middle + max_range])
    
def plot_static_shells(image, mask, center, num_shells=20, points_per_shell=5000, min_radius=None, max_radius=None):
    """
    Generate a static 3D plot of shells:
    - Plot only the left half of each shell.
    - Arrange shells from large to small, left to right.
    - Use MRI intensity values for coloring.
    """
    # Get mask region boundaries if not provided
    if min_radius is None or max_radius is None:
        min_r, max_r, _ = analyze_mask_region(mask, center)
        min_radius = min_radius or min_r
        max_radius = max_radius or max_r

    # Generate shell radii using cubic root spacing
    alpha0 = np.linspace(0, 1, num_shells) ** 1
    # shell_radii = min_radius + (max_radius - min_radius) * alpha
    shell_radii = max_radius - (max_radius - min_radius) * alpha0

    # Calculate the minimum and maximum pixel values in the tumor region
    # Extract pixel values within the tumor mask
    tumor_values = image[mask > 0]
    vmin = np.min(tumor_values)  # Minimum pixel value in the tumor region
    vmax = np.max(tumor_values)  # Maximum pixel value in the tumor region

    # Create a figure and 3D axes
    fig = plt.figure(figsize=(12, 8))
    ax = fig.add_subplot(111, projection='3d')
    
    alpha = [0.2, 0.3, 0.4, 0.5, 0.6, 0.7]
    offset_dist = [0, 1, 0.9, 0.85, 0.8, 0.75]

    # Plot each shell
    for i, radius in enumerate(shell_radii):
        # Generate points for this shell
        phi, theta = generate_uniform_sphere_points(points_per_shell)
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

        # Create a mask for the left half of the shell
        # left_half_mask = (theta >= -np.pi) & (theta <= 0)   # Only plot points where theta < 180 degrees
        left_half_mask = y < center[1]  # keep only points left of center

        # Apply the mask to the coordinates and values
        x_left = x[left_half_mask]
        y_left = y[left_half_mask]
        z_left = z[left_half_mask]
        values_left = values[left_half_mask]

        # Add an offset to separate the shells along the x-axis
        x_offset = 2.7 * i * offset_dist[i] * (max_radius - min_radius) / 3  # Adjust spacing as needed
        y_left += x_offset

        print(f'Current radiu is {radius} and offset is {x_offset}')
        
        from matplotlib.colors import ListedColormap, BoundaryNorm
        cmap = ListedColormap(["#90EE90", "#FFD580", "#FFB6C1"])
        # cmap = ListedColormap(["green", "yellow", "red"])

        # Create scatter plot for the left half of the shell
        scatter = ax.scatter(x_left, y_left, z_left, c=values_left,
                                cmap=cmap, alpha=alpha[i], s=5, vmin=vmin, vmax=vmax)
        
        # scatter = ax.scatter(x_left, y_left, z_left, c=values_left,
                                # cmap='jet', alpha=0.8, s=5, vmin=vmin, vmax=vmax)

    # Set labels and title
    # ax.set_xlabel('X')
    # ax.set_ylabel('Y')
    # ax.set_zlabel('Z')
    # ax.set_title(f'Tumor Region Shells (Left Half, {num_shells} Shells)')
    
    # Calculate the total offset for the outermost shell
    total_offset = (num_shells - 1) * (max_radius - min_radius) / num_shells

    # Set consistent view limits to fully contain the largest shell
    max_range = max_radius  # Use the largest radius to set the axis limits
    # ax.set_xlim([center[0] - max_range, center[0] + max_range + (num_shells - 1) * (max_radius - min_radius) / num_shells * 1.5])
    # ax.set_xlim([center[0] - max_range, center[0] + max_range] + (num_shells - 1) * (max_radius - min_radius) / num_shells * 1)
    # ax.set_ylim([center[1] - max_range, center[1] + max_range] + (num_shells - 1) * (max_radius - min_radius) / num_shells * 1)
    # ax.set_zlim([center[2] - max_range, center[2] + max_range] + (num_shells - 1) * (max_radius - min_radius) / num_shells * 1)

    # Fix the viewing angle
    # ax.view_init(elev=30, azim=60)

    # Ensure equal scaling on all axes
    ax.set_box_aspect([1, 1, 1])
    
    # ax.set_xticks([])           # Remove tick marks
    ax.set_xticklabels([])      # Remove tick labels
    # ax.set_yticks([])           # Remove tick marks
    ax.set_yticklabels([])      # Remove tick labels
    # ax.set_zticks([])           # Remove tick marks
    ax.set_zticklabels([])      # Remove tick labels
    
    # Equal Axis
    set_axes_equal(ax)
    ax.dist = 5   # default is 10, smaller = zoom in
    ax.axis("off")
    ax.view_init(elev=30, azim=20)
    
    # Add colorbar
    cbar = plt.colorbar(scatter, ax=ax, shrink=0.8)
    cbar.set_label('MRI Intensity')

    # Adjust layout
    plt.tight_layout()

    # Show the plot
    plt.show()
    
def plot_surface_shells(image, mask, center, num_shells=1, points_per_shell=5000, min_radius=None, max_radius=None):
    """
    Generate a static 3D plot of shells:
    - Plot only the left half of each shell.
    - Arrange shells from large to small, left to right.
    - Use MRI intensity values for coloring.
    """
    # Get mask region boundaries if not provided
    if min_radius is None or max_radius is None:
        min_r, max_r, _ = analyze_mask_region(mask, center)
        min_radius = min_radius or min_r
        max_radius = max_radius or max_r

    # Generate shell radii using cubic root spacing
    alpha = np.linspace(0, 1, num_shells) ** 1
    # shell_radii = min_radius + (max_radius - min_radius) * alpha
    shell_radii = max_radius - (max_radius - min_radius) * np.array([0.25])
    # shell_radii = max_radius - (max_radius - min_radius) * alpha

    # Calculate the minimum and maximum pixel values in the tumor region
    # Extract pixel values within the tumor mask
    tumor_values = image[mask > 0]
    vmin = np.min(tumor_values)  # Minimum pixel value in the tumor region
    vmax = np.max(tumor_values)  # Maximum pixel value in the tumor region

    # Create a figure and 3D axes
    fig = plt.figure(figsize=(12, 8))
    ax = fig.add_subplot(111, projection='3d')

    # Plot each shell
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

        # Create a mask for the left half of the shell
        # left_half_mask = (theta >= -np.pi) & (theta <= 0)   # Only plot points where theta < 180 degrees

        # Apply the mask to the coordinates and values
        x_full = x
        y_full = y
        z_full = z
        values_full = values

        # Add an offset to separate the shells along the x-axis
        x_offset = 1 * i * (max_radius - min_radius) / 3  # Adjust spacing as needed
        y_full += x_offset

        print(f'Current radiu is {radius} and offset is {x_offset}')
        
        from matplotlib.colors import ListedColormap, BoundaryNorm
        # cmap = ListedColormap(["#90EE90", "#FFD580", "#FFB6C1"])
        # cmap = ListedColormap(["green", "yellow", "red"])

        # Create scatter plot for the left half of the shell
        scatter = ax.scatter(x_full, y_full, z_full, c=values_full,
                                cmap='jet', alpha=0.8, s=10, vmin=vmin, vmax=vmax)
        
        # Plot the full shell as a surface
        # surf = ax.plot_trisurf(
        #     x_full, y_full, z_full, 
        #     triangles=None,       # let matplotlib auto-triangulate
        #     cmap='jet', 
        #     linewidth=0, 
        #     antialiased=True, 
        #     alpha=0.8,
        #     shade=True,
        #     facecolors=plt.cm.jet((values_full - vmin)/(vmax - vmin))  # map values to colormap
        # )
        
        # # Normalize values for colormap
        # vmin = np.min(image)
        # vmax = np.max(image)
        # normalized = (values - vmin) / (vmax - vmin + 1e-8)
        
        # # Plot smooth surface
        # ax.plot_surface(
        #     x, y, z, 
        #     rstride=1, cstride=1, facecolors=plt.cm.jet(normalized),
        #     linewidth=0, antialiased=True, alpha=0.8
        # )


    # Set labels and title
    # ax.set_xlabel('X')
    # ax.set_ylabel('Y')
    # ax.set_zlabel('Z')
    # ax.set_title(f'Tumor Region Shells (Left Half, {num_shells} Shells)')
    
    # Calculate the total offset for the outermost shell
    total_offset = (num_shells - 1) * (max_radius - min_radius) / num_shells

    # Set consistent view limits to fully contain the largest shell
    max_range = max_radius  # Use the largest radius to set the axis limits
    ax.set_xlim([center[0] - max_range, center[0] + max_range])
    ax.set_ylim([center[1] - max_range, center[1] + max_range + (num_shells - 1) * (max_radius - min_radius) / num_shells * 1.5])
    ax.set_zlim([center[2] - max_range, center[2] + max_range])

    # Fix the viewing angle
    ax.view_init(elev=10, azim=45)

    # Ensure equal scaling on all axes
    ax.set_box_aspect([1, 1, 1])
    ax.set_xticklabels([])      # Remove tick labels
    ax.set_yticklabels([])      # Remove tick labels
    ax.set_zticklabels([])      # Remove tick labels
    
    # Add colorbar
    cbar = plt.colorbar(scatter, ax=ax, shrink=0.8)
    cbar.set_label('MRI Intensity')

    # Adjust layout
    plt.tight_layout()

    # Show the plot
    plt.show()

# Sliding-window GLRLM visualization
def sliding_window_glrlm(shell_results, shell_radii, window_size=15, feature='ShortRunEmphasis'):
    """
    Compute and visualize sliding-window GLRLM feature maps for 2D shells.
    
    Args:
        shell_results: list of tuples/lists with (intensity_map, weight_map) or similar
        shell_radii: list of radii corresponding to each shell
        window_size: int, size of sliding window
        feature: str, GLRLM feature name (e.g., 'ShortRunEmphasis', 'LongRunEmphasis')
    """
    
    valid_shells = []
    for i, (radius, result) in enumerate(zip(shell_radii, shell_results)):
        if result is not None and result[0] is not None:
            valid_shells.append((i, radius, result[0], result[1]))

    if not valid_shells:
        print("No valid shells to visualize")
        return

    n_shells = len(valid_shells)
    n_cols = min(5, n_shells)
    n_rows = (n_shells + n_cols - 1) // n_cols

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(4*n_cols, 4*n_rows))
    if n_rows == 1:
        axes = np.array([axes])  # force 2D array for consistent indexing

    for idx, (i, radius, intensity_map, weight_map) in enumerate(valid_shells):
        print(f'Working on {i}')
        row = idx // n_cols
        col = idx % n_cols
        
        pad         = window_size // 2
        padded      = np.pad(intensity_map, pad, mode='reflect')
        feature_map = np.zeros_like(intensity_map, dtype=float)
        
        # configure PyRadiomics to compute only GLRLM features
        params = {
            'binWidth': 25,
            'resampledPixelSpacing': None,
            'interpolator': 'sitkBSpline',
            'enableCExtensions': True
        }
        # extractor = featureextractor.RadiomicsFeatureExtractor(**params)
        extractor = featureextractor.RadiomicsFeatureExtractor()
        extractor.disableAllFeatures()
        extractor.enableFeatureClassByName('glrlm')
        
        for x in range(intensity_map.shape[0]):
            for y in range(intensity_map.shape[1]):
                
                print(f'Working on {i} window {x,y}')
                patch = padded[x:x+window_size, y:y+window_size]
                patch_mask = np.ones_like(patch, dtype=np.uint8)
                patch_mask[:,0] = 0
        
                img_sitk = sitk.GetImageFromArray(patch.astype(np.float32))
                mask_sitk = sitk.GetImageFromArray(patch_mask)
        
                result = extractor.execute(img_sitk, mask_sitk, label=1)
                feature_value = result.get(f'original_glrlm_{feature}', 0.0)
        
                feature_map[x, y] = feature_value

        # Plot the GLRLM feature heatmap, not the raw intensity map
        im = axes[row, col].imshow(feature_map,
                                   extent=[-np.pi, np.pi, 0, np.pi],
                                   aspect='auto',
                                   cmap='hot')
        axes[row, col].set_title(f'Shell {i}\nRadius: {radius:.2f}', fontsize=16)
        axes[row, col].tick_params(axis='x', labelsize=14)
        axes[row, col].tick_params(axis='y', labelsize=14)
        plt.colorbar(im, ax=axes[row, col])

    # Clear unused subplots
    for idx in range(len(valid_shells), n_rows * n_cols):
        row = idx // n_cols
        col = idx % n_cols
        axes[row, col].axis('off')

    plt.tight_layout()
    plt.show()
    
    print(f'Finish one shell')

# Sliding-window GLRLM visualization
def sliding_window_glcm(shell_results, shell_radii, window_size=15, feature='Contrast'):
    """
    Compute and visualize sliding-window GLRLM feature maps for 2D shells.
    
    Args:
        shell_results: list of tuples/lists with (intensity_map, weight_map) or similar
        shell_radii: list of radii corresponding to each shell
        window_size: int, size of sliding window
        feature: str, GLRLM feature name (e.g., 'ShortRunEmphasis', 'LongRunEmphasis')
    """
    
    valid_shells = []
    for i, (radius, result) in enumerate(zip(shell_radii, shell_results)):
        if result is not None and result[0] is not None:
            valid_shells.append((i, radius, result[0], result[1]))

    if not valid_shells:
        print("No valid shells to visualize")
        return

    n_shells = len(valid_shells)
    n_cols = min(5, n_shells)
    n_rows = (n_shells + n_cols - 1) // n_cols

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(4*n_cols, 4*n_rows))
    if n_rows == 1:
        axes = np.array([axes])  # force 2D array for consistent indexing

    for idx, (i, radius, intensity_map, weight_map) in enumerate(valid_shells):
        print(f'Working on {i}')
        row = idx // n_cols
        col = idx % n_cols
        
        pad         = window_size // 2
        padded      = np.pad(intensity_map, pad, mode='reflect')
        feature_map = np.zeros_like(intensity_map, dtype=float)
        
        # configure PyRadiomics to compute only GLRLM features
        params = {
            'binWidth': 25,
            'resampledPixelSpacing': None,
            'interpolator': 'sitkBSpline',
            'enableCExtensions': True
        }
        # extractor = featureextractor.RadiomicsFeatureExtractor(**params)
        extractor = featureextractor.RadiomicsFeatureExtractor()
        extractor.disableAllFeatures()
        extractor.enableFeatureClassByName('glcm')
        
        for x in range(intensity_map.shape[0]):
            for y in range(intensity_map.shape[1]):
                
                print(f'Working on {i} window {x,y}')
                patch = padded[x:x+window_size, y:y+window_size]
                patch_mask = np.ones_like(patch, dtype=np.uint8)
                patch_mask[:,0] = 0
        
                img_sitk = sitk.GetImageFromArray(patch.astype(np.float32))
                mask_sitk = sitk.GetImageFromArray(patch_mask)
        
                result = extractor.execute(img_sitk, mask_sitk, label=1)
                feature_value = result.get(f'original_glcm_{feature}', 0.0)
        
                feature_map[x, y] = feature_value

        # Plot the GLRLM feature heatmap, not the raw intensity map
        im = axes[row, col].imshow(feature_map,
                                   extent=[-np.pi, np.pi, 0, np.pi],
                                   aspect='auto',
                                   cmap='hot')
        axes[row, col].set_title(f'Shell {i}\nRadius: {radius:.2f}', fontsize=16)
        axes[row, col].tick_params(axis='x', labelsize=14)
        axes[row, col].tick_params(axis='y', labelsize=14)
        plt.colorbar(im, ax=axes[row, col])

    # Clear unused subplots
    for idx in range(len(valid_shells), n_rows * n_cols):
        row = idx // n_cols
        col = idx % n_cols
        axes[row, col].axis('off')

    plt.tight_layout()
    plt.show()
    
    print(f'Finish one shell')


# %% Load the file number IDs


file_loc = '/Users/hafeng/Documents/Research_Data/UCSF_Dataset/UCSF_PDGM/PDGM/PKG-UCSF-PDGM-v3/UCSF-PDGM-v3/'
output_loc = '/Users/hafeng/Documents/Postdoc_Research_Meetings/TCGA/UCSF_PDGM/radiogenomic_batch/'
num_blocks = 8

number_list = []

for filename in os.listdir(file_loc):

    match = re.search(r'PDGM-(\d{4})_nifti', filename)
    if match:
        number = match.group(1)  # Get the first capturing group
        number_list.append(int(number))

number_list.sort()

# %% Read the files
output_folder = "/Users/hafeng/Documents/Postdoc_Research_Meetings/TCGA/UCSF_PDGM/radiogenomic_shell_worldmap/"
num_shell = 1

extractor = featureextractor.RadiomicsFeatureExtractor(preCrop=True)
extractor.settings['label'] = 1  # Set the label for tumor regions
extractor.settings['resampledPixelSpacing'] = [1, 1]
extractor.settings['force2D'] = True
# For demonstration, process only the first item
for idx in range(1):

    item = number_list[idx]
    item = str(64).zfill(4)
    folder_name = 'UCSF-PDGM-' + item + '_nifti/'
    item_path = os.path.join(file_loc, folder_name)

    adc_path = 'UCSF-PDGM-' + item + '_ADC.nii.gz'
    ce_path = 'UCSF-PDGM-' + item + '_T1c.nii.gz'
    fl_path = 'UCSF-PDGM-' + item + '_FLAIR.nii.gz'
    brain_mask = 'UCSF-PDGM-' + item + '_brain_segmentation.nii.gz'
    tumor_mask = 'UCSF-PDGM-' + item + '_tumor_segmentation.nii.gz'

    nii_adc = nib.load(os.path.join(item_path, adc_path))
    nii_ce = nib.load(os.path.join(item_path, ce_path))
    nii_fl = nib.load(os.path.join(item_path, fl_path))
    nii_brain = nib.load(os.path.join(item_path, brain_mask))
    nii_tumor = nib.load(os.path.join(item_path, tumor_mask))

    # Get the image data as a numpy array
    adc_data = nii_adc.get_fdata()
    ce_data = nii_ce.get_fdata()
    fl_data = nii_fl.get_fdata()
    tumor_data = nii_tumor.get_fdata()

    # Load the medical image and tumor mask
    image_ce = ce_data
    image_fl = fl_data
    image_adc = adc_data
    image_path = os.path.join(item_path, ce_path)
    mask = tumor_data
    # mask       = np.ones((240,240,155))
    # Compute tumor center
    # tumor_center = get_tumor_center(mask==1)

    # Region labels in the tumor mask
    region_labels = {
        "necrotic": 1,
        "enhancing": 4,
        "lesion": 2
    }

    # Process each region
    for region_name, label in region_labels.items():
        print(f"\nProcessing {region_name} region (label {label})...")
        print(f"\nWorking on CE Data")
        # num_shells, intensity, weight
        
        mask_ce = mask
        # results_ce = process_region(image_ce, mask, label, region_name, num_shell, 0, 6000)
        results_ce = process_region(mask_ce, mask, label, region_name, num_shell, 0, 6000)
        # print(f"\nWorking on FL Data")
        # results_fl = process_region(image_fl, mask, label, region_name, num_shell, 0, 2500)   # num_shells, intensity, weight
        # results_adc = process_region(image_adc, mask, label, region_name, num_shell, 0, 0.005)
        
#%% Plot out the original model

image_norm = image_ce

fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')

ax.voxels(image_ce, facecolors=plt.cm.viridis(image_ce), edgecolor='none')

#%%

from skimage import measure
import matplotlib.pyplot as plt

# Create isosurface
threshold = np.percentile(image_ce[image_ce > 0], 0.01)  # Use 70th percentile
verts, faces, _, _ = measure.marching_cubes(image_ce, level=threshold)

fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')

# Plot the surface
ax.plot_trisurf(verts[:, 0], verts[:, 1], faces, verts[:, 2], 
                cmap='viridis', alpha=0.5, linewidth=0)

# ax.set_xticks([])           # Remove tick marks
ax.set_xticklabels([])      # Remove tick labels
# ax.set_yticks([])           # Remove tick marks
ax.set_yticklabels([])      # Remove tick labels
# ax.set_zticks([])           # Remove tick marks
ax.set_zticklabels([])      # Remove tick labels