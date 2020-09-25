% ORIGINAL CODE FROM Linden et al. 2020 (Go with the FLOW: Visualizing spatiotemporal dynamics in1optical widefield calcium imaging)

% This script recreates all the panes in Figure 1
% 
% NOTE: The colormaps here use system colormaps and differ from those in than the manuscript. See colormaps section below for details.

clear all; close all; clc

%%%% ADD REQUIRED CODE %%%%
addpath('../flow_portraits/')
addpath('~/Dropbox/Science/DrosteEffect-BrewerMap-5b84f95/')
mkdir('./figure1/')

%%%% COLORMAPS %%%%
%   This script uses MATLAB colormaps by defualt, however these were not used in the manuscript.
cmap_svd = colormap('winter');
cmap_nnmf = colormap('cool');
cmap_dfof = colormap('gray');

% To exactly recreate the images in the manuscript please download and add ColorBrewer and Colorcet to the path:
%   https://www.mathworks.com/matlabcentral/fileexchange/45208-colorbrewer-attractive-and-distinctive-colormaps
%   https://peterkovesi.com/projects/colourmaps/
%   
%   AND uncomment below
%
% cmap_svd = flip(brewermap(64,'rdbu'),1);
% cmap_nnmf = colorcet('l12');
% cmap_dfof = colorcet('l5');
%
%

%%%% CONSTANTS %%%%
num_factors = 4;
figsize = [0 0 1000 1000]; 

%%%% LOAD and Preprocess DATA %%%

%%%% Developing Mouse Data (P7/8 mouse) %%%%
% Data from Dennis Tabuena and William J. dev (Tabuena et al. 2019) - University of Washington
load('../data/developing_wave2.mat');
dfof_dev = dat.dfof; mask_dev = dat.mask;
clear dat

dims_dev = size(dfof_dev);
dfof_dev_2d = reshape(dfof_dev, [dims_dev(1)*dims_dev(2),dims_dev(3)]);

% compute the SVD of the data
[u_dev, s_dev, v_dev] = svd(dfof_dev_2d - mean(dfof_dev_2d,2), 'econ'); % economy and mean-centered

% compute the NNMF
[w_dev, h_dev] = nnmf(dfof_dev_2d, num_factors);

% compute and save FLOW portrait
% this is the thrid column of Figure 1B for the develping data
filename = './figure1/developing_flowPortrait.png';
integrationLength = 40;
flowPortrait(dfof_dev, integrationLength, 'filename', filename,'white_background', true, 'mask', mask_dev);


%%%% Adult Mouse Data  %%%%
% Data from Zatka-Haas et al. 2020 - University College London

load('../data/fig1_adult.mat')
dfof_adult = dat.dfof; mask_adult = dat.mask;
clear dat

dims_adult = size(dfof_adult);
dfof_adult_2d = reshape(dfof_adult, [dims_adult(1)*dims_adult(2),dims_adult(3)]); 

% compute the SVD of the data
[u_adult, s_adult, v_adult] = svd(dfof_adult_2d - mean(dfof_adult_2d,2), 'econ'); % economy and mean-centered

% compute the NNMF
[w_adult, h_adult] = nnmf(dfof_adult_2d, num_factors);

% compute and save FLOW portrait
% this is the thrid column of Figure 1B for the develping data
filename = './figure1/adult_flowPortrait.png';
integrationLength = 20;
history_delay = 15;
flowPortrait(dfof_adult, integrationLength, 'filename', filename,'history_delay', history_delay, 'white_background', true, 'mask', mask_adult);

%%%% PLOT SVD SPATIAL MODES %%%%
% This is the 1st column of Figure 1B
dev_name = './figure1/svd_dev_';
adult_name = './figure1/svd_steinmetz_';

for idx = 1:4
    % dev
    dat = reshape(u_dev(:,idx), [dims_dev(1),dims_dev(2)]);
    [dat, imAlpha] = whiteBackground(dat, mask_dev);
    figure('Renderer', 'painters', 'Position', figsize,'color', 'w');
    imagesc(dat, 'AlphaData',imAlpha);
    ax = gca; ax.CLim = [-1*max(max(abs(dat))) max(max(abs(dat)))];
    ax.Position = ax.OuterPosition;
    colormap(cmap_svd); axis off; axis equal
    saveas(gcf, strcat(dev_name, num2str(idx),'.png'));
    close all

    % adult
    dat = reshape(u_adult(:,idx), [dims_adult(1),dims_adult(2)]);
    figure('Renderer', 'painters', 'Position', figsize,'color', 'w');
    [dat, imAlpha] = whiteBackground(dat, mask_adult);
    figure('Renderer', 'painters', 'Position', figsize,'color', 'w');
    imagesc(dat, 'AlphaData',imAlpha);
    ax = gca; ax.CLim = [-1*max(max(abs(dat))) max(max(abs(dat)))];
    ax.Position = ax.OuterPosition;
    colormap(cmap_svd); axis off; axis equal
    saveas(gcf, strcat(adult_name, num2str(idx),'.png'));
    close all
end

%%%% PLOT NNMF SPATIAL MODES %%%%
dev_name = './figure1/nnmf_dev_';
adult_name = './figure1/nnmf_adult_';

for idx = 1:4
    % dev
    dat = reshape(w_dev(:,idx), [dims_dev(1),dims_dev(2)]);
    [dat, imAlpha] = whiteBackground(dat, mask_dev);
    figure('Renderer', 'painters', 'Position', figsize,'color', 'w');
    clim = max(max(abs(dat)))*[-1 1];
    imagesc(dat, 'AlphaData',imAlpha);
    ax = gca; ax.CLim = clim;
    ax.Position = ax.OuterPosition;
    % cb = colorbar; cb.Ticks = [];
    colormap(cmap_svd); axis off; axis equal
    saveas(gcf, strcat(dev_name, num2str(idx),'.png'));
    close all

    % adult
    dat = reshape(w_adult(:,idx), [dims_adult(1),dims_adult(2)]);
    [dat, imAlpha] = whiteBackground(dat, mask_adult);
    figure('Renderer', 'painters', 'Position', figsize,'color', 'w');
    clim = max(max(abs(dat)))*[-1 1];
    imagesc(dat, 'AlphaData',imAlpha);
    ax = gca; ax.CLim = clim;
    ax.Position = ax.OuterPosition;
    colormap(cmap_svd); axis off; axis equal
    saveas(gcf, strcat(adult_name, num2str(idx),'.png'));
    close all
end

%%%% SAVE INDIVDUAL PANES OF DATA %%%%%
% Figure 1A
dev_name = './figure1/dev_';
adult_name = './figure1/steinmetz_';

stride_dev = 20;
stride_adult = 17;

idxs = [];
for i = 1:7
    idx = i + ((i-1)*stride_dev);
    idxs = [idx, idxs];

    % devloping mouse
    dat = dfof_dev(:,:,idx); 
    [dat, imAlpha] = whiteBackground(dat, mask_dev);
    figure('Renderer', 'painters', 'Position', figsize,'color', 'w');
    imagesc(dat,'AlphaData',imAlpha); hold on
    if i==1
        plot([1, 14.08], [dims_dev(1)-2 dims_dev(1)-2],'k','LineWidth',5); % add scale bar to first image
    end
    colormap(cmap_dfof)
    ax = gca; 
    ax.Position = ax.OuterPosition;
    axis off; axis equal
    saveas(gcf, strcat(dev_name, num2str(i),'.png'));
    close all

    % adult mosue
    idx = i + ((i-1)*stride_adult);
    dat = double(dFoF_adult(:,:,idx));
    figure('Renderer', 'painters', 'Position', figsize,'color', 'w');
    [dat, imAlpha] = whiteBackground(dat, mask_adult);
    imagesc(dat, 'AlphaData',imAlpha); hold on
    if i==1
        val = 7;
        plot([1, 18.52], [dims_adult(1)-val dims_adult(1)-val],'k','LineWidth',12); % add scale bar to first image
    end
    ax = gca;
    ax.Position = ax.OuterPosition;
    colormap(cmap_dfof); axis off; axis equal
    saveas(gcf, strcat(adult_name, num2str(i),'.png'));
    close all
end