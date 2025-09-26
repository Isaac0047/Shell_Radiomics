%% Load your NIfTI files
% You can use niftiread and niftiinfo


adc_data   = niftiread('/Users/hafeng/Documents/Research_Data/UCSF_Dataset/UCSF_PDGM/PDGM/PKG-UCSF-PDGM-v3/UCSF-PDGM-v3/UCSF-PDGM-0010_nifti/UCSF-PDGM-0010_ADC.nii.gz');
ce_data    = niftiread('/Users/hafeng/Documents/Research_Data/UCSF_Dataset/UCSF_PDGM/PDGM/PKG-UCSF-PDGM-v3/UCSF-PDGM-v3/UCSF-PDGM-0010_nifti/UCSF-PDGM-0010_T1c.nii.gz');
fl_data    = niftiread('/Users/hafeng/Documents/Research_Data/UCSF_Dataset/UCSF_PDGM/PDGM/PKG-UCSF-PDGM-v3/UCSF-PDGM-v3/UCSF-PDGM-0010_nifti/UCSF-PDGM-0010_FLAIR.nii.gz');
brain_mask = niftiread('/Users/hafeng/Documents/Research_Data/UCSF_Dataset/UCSF_PDGM/PDGM/PKG-UCSF-PDGM-v3/UCSF-PDGM-v3/UCSF-PDGM-0010_nifti/UCSF-PDGM-0010_brain_segmentation.nii.gz');
tumor_mask = niftiread('/Users/hafeng/Documents/Research_Data/UCSF_Dataset/UCSF_PDGM/PDGM/PKG-UCSF-PDGM-v3/UCSF-PDGM-v3/UCSF-PDGM-0010_nifti/UCSF-PDGM-0010_tumor_segmentation.nii.gz');

brain_mask = double(brain_mask > 0);
tumor_mask = double(tumor_mask > 0);

%% Generate spherical shells
function [shells, centroid] = generate_solid_shells(tumor_mask, radii)
    [Z, Y, X] = ndgrid(1:size(tumor_mask,1), 1:size(tumor_mask,2), 1:size(tumor_mask,3));
    coords = find(tumor_mask);
    [zc, yc, xc] = ind2sub(size(tumor_mask), coords);
    centroid = [mean(zc), mean(yc), mean(xc)];
    
    dist = sqrt( (X - centroid(3)).^2 + (Y - centroid(2)).^2 + (Z - centroid(1)).^2 );
    
    n_shells = length(radii);
    shells = cell(1, n_shells);
    for i = 1:n_shells
        if i == 1
            shell = dist <= radii(i);
        else
            shell = (dist > radii(i-1)) & (dist <= radii(i));
        end
        shells{i} = shell;
    end
end

radii = [4 8 12 16 20 24];
[shells, centroid] = generate_solid_shells(tumor_mask, radii);

%% Plot shells with distinct colors
figure('Position',[100 100 1200 800]);
colors = [1 0 0; 1 0 1; 1 1 0; 0 1 0; 1 0.75 0.8; 0 0 1; 0.5 0 0.5]; % RGB

% Plot brain outline (optional)
p = patch(isosurface(brain_mask,0.5));
isonormals(brain_mask,p);
p.FaceColor = [0.8 0.8 0.8]; % light gray
p.FaceAlpha = 0.1;
p.EdgeColor = 'none';
hold on;

% Plot tumor mask
p = patch(isosurface(tumor_mask,0.5));
isonormals(tumor_mask,p);
p.FaceColor = [0.5 0.8 1]; % light blue
p.FaceAlpha = 0.3;
p.EdgeColor = 'none';

% Plot shells
for i = 1:length(shells)
    shell = shells{i};
    p = patch(isosurface(shell,0.5));
    isonormals(shell,p);
    p.FaceColor = colors(mod(i-1,size(colors,1))+1,:);
    p.FaceAlpha = 1;  % fully opaque to prevent blending issues
    p.EdgeColor = 'k';
    p.LineWidth = 0.5;
end

axis equal;
axis off;
xlabel('X'); ylabel('Y'); zlabel('Z');
view(3);
camlight; lighting gouraud;
%title('3D Tumor Shells with Brain Outline');

%% Plot shells with distinct colors
% Define shell colors (adjust or add more)
shell_colors = [1 0 0;     % red
                1 0.5 0;   % orange
                1 1 0;     % yellow
                0 1 0;     % green
                1 0.75 0.8; % pink
                0 0 1;     % blue
                0.5 0 0.5]; % purple

% Cutout Visualization centered on tumor
figure('Position',[100 100 1800 600]);

% Cut distance from tumor centroid (half of cut length)
cut_dist = 0;  % adjust as needed

% Loop over cutout directions: front(Z), side(Y), top(X)
cut_axes = {'Z','Y','X'};
for k = 1:3
    subplot(1,3,k);
    hold on;

    % Copy masks
    brain_cut = brain_mask;
    tumor_cut = tumor_mask;

    % Apply cut relative to tumor centroid
    switch k
        case 1 % Front cut along Z
            cut_idx = round(centroid(3) - cut_dist);
            brain_cut(1:cut_idx,:,:) = 0;
            tumor_cut(1:cut_idx,:,:) = 0;
        case 2 % Side cut along Y
            cut_idx = round(centroid(2) - cut_dist);
            brain_cut(:,1:cut_idx,:) = 0;
            tumor_cut(:,1:cut_idx,:) = 0;
        case 3 % Top cut along X
            cut_idx = round(centroid(1) - cut_dist);
            brain_cut(:,:,1:cut_idx) = 0;
            tumor_cut(:,:,1:cut_idx) = 0;
    end

    % Plot brain outline
    p = patch(isosurface(brain_cut,0.5));
    isonormals(brain_cut,p);
    p.FaceColor = [0.8 0.8 0.8];
    p.FaceAlpha = 0.1;
    p.EdgeColor = 'none';

    % Plot tumor
    p = patch(isosurface(tumor_cut,0.5));
    isonormals(tumor_cut,p);
    p.FaceColor = [0.5 0.8 1];
    p.FaceAlpha = 0.3;
    p.EdgeColor = 'none';

    % Plot shells with individual colors and optional offset
    for i = 1:length(shells)
        shell_cut = shells{i};
        % Apply cut
        switch k
            case 1, shell_cut(1:cut_idx,:,:) = 0;
            case 2, shell_cut(:,1:cut_idx,:) = 0;
            case 3, shell_cut(:,:,1:cut_idx) = 0;
        end
        % Optional: slightly offset inner shells for visibility
        % shell_cut = circshift(shell_cut, [0, 0, i-1]); % uncomment if needed

        p = patch(isosurface(shell_cut,0.5));
        isonormals(shell_cut,p);
        color_idx = mod(i-1,size(shell_colors,1)) + 1;
        p.FaceColor = shell_colors(color_idx,:);
        p.FaceAlpha = 0.7;  % adjust transparency
        p.EdgeColor = 'none';
        p.LineWidth = 0.5;
    end

    axis equal; axis tight;
    xlabel('X'); ylabel('Y'); zlabel('Z');
    camlight; lighting gouraud;

    switch k
        case 1, title('Front Cutout (Z-axis, tumor center)');
        case 2, title('Side Cutout (Y-axis, tumor center)');
        case 3, title('Top Cutout (X-axis, tumor center)');
    end
end

%% 3D view of the cutout region

% Assume tumor_mask, brain_mask are loaded (3D binary masks)

% Compute tumor centroid
[idxY, idxX, idxZ] = ind2sub(size(tumor_mask), find(tumor_mask>0));
centroid = [mean(idxY), mean(idxX), mean(idxZ)];  % [Y X Z] in MATLAB

% Define shell radii
radii = [4, 8, 12, 16, 20, 24];
colors = [1 0 0; 1 0.5 0; 1 1 0; 0 1 0; 1 0 1; 0 0 1]; % RGB colors

[dimY, dimX, dimZ] = size(tumor_mask);
[YY, XX, ZZ] = ndgrid(1:dimY, 1:dimX, 1:dimZ);

% Distance from tumor centroid
dist = sqrt((XX - centroid(2)).^2 + (YY - centroid(1)).^2 + (ZZ - centroid(3)).^2);

% Generate shells
shells = false([dimY, dimX, dimZ, numel(radii)]);
for i = 1:numel(radii)
    if i == 1
        shells(:,:,:,i) = dist <= radii(i);
    else
        shells(:,:,:,i) = (dist > radii(i-1)) & (dist <= radii(i));
    end
end

% Function to plot cutout view
function plot_cutout_shells(shells, brain_mask, tumor_mask, cut_axis, centroid, colors)
    figure; hold on;
    alpha_inner = linspace(0.8, 0.3, size(shells,4)); % transparency for shells
    for i = 1:size(shells,4)
        mask = shells(:,:,:,i);

        % Apply cutout
        switch cut_axis
            case 'Z' % front cut
                mask(:,:,1:round(centroid(3))) = 0;
            case 'Y' % side cut
                mask(1:round(centroid(1)),:,:) = 0;
            case 'X' % top cut
                mask(:,1:round(centroid(2)),:) = 0;
        end

        % Smooth mask
        mask_smooth = smooth3(mask,'gaussian',3);

        % Generate isosurface
        f = patch(isosurface(mask_smooth,0.5));
        isonormals(mask_smooth,f)
        f.FaceColor = colors(i,:);
        f.EdgeColor = 'none';
        f.FaceAlpha = alpha_inner(i);
    end

    % Plot brain outline (optional)
    brain_smooth = smooth3(brain_mask,'gaussian',3);
    fBrain = patch(isosurface(brain_smooth,0.5));
    isonormals(brain_smooth,fBrain)
    fBrain.FaceColor = [0.8 0.8 0.8];
    fBrain.EdgeColor = 'none';
    fBrain.FaceAlpha = 0.1;

    % Plot brain outline (optional)
    tumor_smooth = smooth3(tumor_mask,'gaussian',3);
    fTumor = patch(isosurface(tumor_smooth,0.5));
    isonormals(tumor_smooth,fTumor)
    fTumor.FaceColor = [0 0 1];
    fTumor.EdgeColor = 'none';
    fTumor.FaceAlpha = 0.1;

    axis equal; 
    axis off;
    view(3); camlight; lighting gouraud;
    xlabel('X'); ylabel('Y'); zlabel('Z');
end

% Plot 3 cutout views
plot_cutout_shells(shells, brain_mask, tumor_mask, 'Z', centroid, colors) % front
plot_cutout_shells(shells, brain_mask, tumor_mask, 'Y', centroid, colors) % side
plot_cutout_shells(shells, brain_mask, tumor_mask, 'X', centroid, colors) % top

%% quarter_cut

function plot_shell_quarter_cut(shells, brain_mask, tumor_mask, centroid, colors)

figure('Position',[100 100 1200 800]);

% Plot brain outline (optional)
p = patch(isosurface(brain_mask,0.5));
isonormals(brain_mask,p);
p.FaceColor = [0.8 0.8 0.8]; % light gray
p.FaceAlpha = 0.1;
p.EdgeColor = 'none';
hold on;

% Plot tumor mask
p = patch(isosurface(tumor_mask,0.5));
isonormals(tumor_mask,p);
p.FaceColor = [0.5 0.8 1]; % light blue
p.FaceAlpha = 0.1;
p.EdgeColor = 'none';

% Define cut planes (through tumor centroid)
xc = centroid(1);
yc = centroid(2);

% Plot shells with quarter cut
for i = 1:length(shells)
    shell = shells{i};

    % Create cut mask (remove 1/4)
    [X,Y,Z] = ndgrid(1:size(shell,1), 1:size(shell,2), 1:size(shell,3));
    cutMask = (X >= xc | Y >= yc);  % keep 3/4, cut 1/4

    cutShell = shell .* cutMask;

    % Plot surface
    p = patch(isosurface(cutShell,0.5));
    isonormals(cutShell,p);
    p.FaceColor = colors(mod(i-1,size(colors,1))+1,:);
    p.FaceAlpha = 1.0;
    p.EdgeColor = 'none';
    p.LineWidth = 0.5;
end

axis equal off;
view(3);
camlight; lighting gouraud;
title('3D Tumor Shells with Quarter Cutout');
end

% Plot 3 cutout views
plot_shell_quarter_cut(shells, brain_mask, tumor_mask, centroid, colors) % front

%% 5/8 Cutout

function plot_shell_cut_5over8(shells, brain_mask, tumor_mask, centroid, colors)

figure('Position',[100 100 1200 800]);

% Plot brain outline (optional)
p = patch(isosurface(brain_mask,0.5));
isonormals(brain_mask,p);
p.FaceColor = [0.8 0.8 0.8]; % light gray
p.FaceAlpha = 0.1;
p.EdgeColor = 'none';
hold on;

% Plot tumor mask
p = patch(isosurface(tumor_mask,0.5));
isonormals(tumor_mask,p);
p.FaceColor = [0.5 0.8 1]; % light blue
p.FaceAlpha = 0.2;
p.EdgeColor = 'none';

% Define cut planes through tumor centroid
xc = centroid(1);
yc = centroid(2);
zc = centroid(3);

% Grid
[X,Y,Z] = ndgrid(1:size(brain_mask,1), 1:size(brain_mask,2), 1:size(brain_mask,3));

% Define full hemisphere: keep everything above centroid in Z
hemiMask = (X >= xc);

% Define extra 1/8 chunk: below centroid, and X >= xc, Y >= yc
extraEighth = (X < xc & Y >= yc & Z <= zc);

% Combine: hemisphere + extra 1/8
finalMask = hemiMask | extraEighth;

% Plot shells
for i = 1:length(shells)
    shell = shells{i};
    cutShell = shell .* finalMask;

    p = patch(isosurface(cutShell,0.5));
    isonormals(cutShell,p);
    p.FaceColor = colors(mod(i-1,size(colors,1))+1,:);
    p.FaceAlpha = 1.0;
    p.EdgeColor = 'none';
    p.LineWidth = 0.5;
end

axis equal off;
%axis equal;
%xlabel('X'); ylabel('Y'); zlabel('Z');
view(3);
camlight; lighting gouraud;
title('3D Tumor Shells with Hemisphere + Extra 1/8 Cutout');
end

% Plot 3 cutout views
plot_shell_cut_5over8(shells, brain_mask, tumor_mask, centroid, colors) % front


%% Another test

function plot_shell_cutout_view(shells, brain_mask, tumor_mask, centroid, colors, view_dir)
    figure; hold on;

    nShells = numel(shells);
    alpha_inner = linspace(0.8,0.6,nShells);

    % Normalize view direction
    view_dir = view_dir / norm(view_dir);

    % Match ndgrid indexing
    centroid_grid = [centroid(1), centroid(2), centroid(3)];

    % Precompute voxel grid
    [YY, XX, ZZ] = ndgrid(1:size(brain_mask,1), 1:size(brain_mask,2), 1:size(brain_mask,3));
    dp_grid = (YY - centroid_grid(1)) * view_dir(1) + ...
              (XX - centroid_grid(2)) * view_dir(2) + ...
              (ZZ - centroid_grid(3)) * view_dir(3);

    % --- Plot tumor mask (cutout) ---
    tumor_cut = tumor_mask;
    tumor_cut(dp_grid > 0) = 0;   % remove back half
    tumor_smooth = smooth3(tumor_cut,'gaussian',3);
    fTumor = patch(isosurface(tumor_smooth,0.5));
    isonormals(tumor_smooth,fTumor)
    fTumor.FaceColor = [1 0 0];   % bright red
    fTumor.EdgeColor = 'none';
    fTumor.FaceAlpha = 0;

    % --- Plot shells (cutout) ---
    for i = 1:nShells
        mask = shells{i};
        mask(dp_grid < 0) = 0;   % cut by plane

        mask_smooth = smooth3(mask,'gaussian',3);
        f = patch(isosurface(mask_smooth,0.5));
        isonormals(mask_smooth,f)
        f.FaceColor = colors(i,:);
        f.EdgeColor = 'none';
        f.FaceAlpha = alpha_inner(i);
    end

    % --- Brain outline (not cut, just reference) ---
    brain_cut = brain_mask;
    brain_cut(dp_grid > 0) = 0; 
    brain_smooth = smooth3(brain_cut,'gaussian',3);
    fBrain = patch(isosurface(brain_smooth,0.5));
    isonormals(brain_smooth,fBrain)
    fBrain.FaceColor = [0.8 0.8 0.8];
    fBrain.EdgeColor = 'none';
    fBrain.FaceAlpha = 0;

    axis equal; camlight; lighting gouraud;
    xlabel('X'); ylabel('Y'); zlabel('Z');
end

% Suppose you already have: shells, brain_mask, tumor_mask, centroid
%colors = lines(numel(shells));  % distinct colors
colors = [
    1.0, 0.0, 0.0;   % red
    1.0, 0.5, 0.0;   % orange
    1.0, 1.0, 0.0;   % yellow
    0.0, 1.0, 0.0;   % green
    1.0, 0.0, 1.0;   % pink/magenta
    0.5, 0.0, 0.5    % purple
];
% Cut along X, Y, Z directions (all through tumor center)
plot_shell_cutout_view(shells, brain_mask, tumor_mask, centroid, colors, [1 0 0]); % cut along X
plot_shell_cutout_view(shells, brain_mask, tumor_mask, centroid, colors, [0 1 0]); % cut along Y
plot_shell_cutout_view(shells, brain_mask, tumor_mask, centroid, colors, [0 0 1]); % cut along Z

%%
function plot_shell_cutaway(shells, brain_mask, tumor_mask, centroid, cut_dir)
    % shells: cell array of logical 3D masks (inner→outer)
    % brain_mask: logical 3D mask
    % tumor_mask: logical 3D mask
    % centroid: [x y z]
    % cut_dir: [1 0 0] (X), [0 1 0] (Y), or [0 0 1] (Z)

    % High-contrast shell colors
    shell_colors = {
        [1 0 0],       ... % red
        [1 0.5 0],     ... % orange
        [1 1 0],       ... % yellow
        [0 1 0],       ... % green
        [1 0.4 0.7],   ... % pink
        [0.6 0 0.8]    ... % purple
    };

    figure; hold on;
    set(gcf,'Color','w'); axis equal; axis off;

    % Grid
    [X,Y,Z] = ndgrid(1:size(brain_mask,1), ...
                     1:size(brain_mask,2), ...
                     1:size(brain_mask,3));

    % Signed distance to cut plane
    dp = (X-centroid(1))*cut_dir(1) + ...
         (Y-centroid(2))*cut_dir(2) + ...
         (Z-centroid(3))*cut_dir(3);

    cut_mask = dp >= 0; % keep only one half

    % Brain
    brain_cut = brain_mask & cut_mask;
    if any(brain_cut(:))
        p = patch(isosurface(brain_cut,0.5));
        set(p,'FaceColor',[0.8 0.8 0.8],'FaceAlpha',0.1,'EdgeColor','none');
    end

    % Tumor
    tumor_cut = tumor_mask & cut_mask;
    if any(tumor_cut(:))
        p = patch(isosurface(tumor_cut,0.5));
        set(p,'FaceColor',[0 0 1],'FaceAlpha',0.5,'EdgeColor','none');
    end

    % Shells (each with distinct color, cut open)
    for i = 1:numel(shells)
        shell_cut = shells{i} & cut_mask;
        if any(shell_cut(:))
            p = patch(isosurface(shell_cut,0.5));
            set(p,'FaceColor',shell_colors{mod(i-1,numel(shell_colors))+1}, ...
                  'FaceAlpha',1.0, 'EdgeColor','none');
        end
    end

    camlight; lighting gouraud;
    view(3);
    title('Tumor-centered Cut-away Shell View', 'FontSize', 14);
end

plot_shell_cutaway(shells, brain_mask, tumor_mask, centroid, [1 0 0]) % cut along X
plot_shell_cutaway(shells, brain_mask, tumor_mask, centroid, [0 1 0]) % cut along Y
plot_shell_cutaway(shells, brain_mask, tumor_mask, centroid, [0 0 1]) % cut along Z