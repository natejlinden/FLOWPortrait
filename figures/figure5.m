% ORIGINAL CODE FROM Linden et al. 2020 (Go with the FLOW: Visualizing spatiotemporal dynamics in1optical widefield calcium imaging)

% This script recreates all the panes in Figure 5
 
% NOTE: The colormaps here use system colormaps and differ from those in than the manuscript. See colormaps section below for details.

clear all; close all; clc

%%%% ADD REQUIRED CODE %%%% 
addpath('../flow_portraits/')
mkdir('./figure5/')


%%%% COLORMAPS %%%%
%   This script uses MATLAB colormaps by defualt, however these were not used in the manuscript.
cmap = colormap('gray');

% To exactly recreate the images in the manuscript please download and add ColorBrewer and Colorcet to the path:
%   https://www.mathworks.com/matlabcentral/fileexchange/45208-colorbrewer-attractive-and-distinctive-colormaps
%   https://peterkovesi.com/projects/colourmaps/
%   
%   AND uncomment below
%
% cmap = colorcet('L5');

%%%% load data %%%%
% Pup Wave
load('../data/fig5_pupPlaneWave.mat')
% Traveling wave:
load('../data/PlaneWave.mat');


%---- Save panes of data ----%
planeFilename = './figure5/planeWave_';
pupFilename = './figure5/pupWave_';

planeStride = 5; planeIdxs = 15:planeStride:40;
pupStride = 5; pupIdxs = 10:pupStride:size(pupWave,3); 
num_panes = 6;


for i = 1:num_panes

    %%%% plane wave
    figure('Renderer', 'painters', 'Position', [0 0 1000 1000],'color', 'w');
    imagesc(planeWave(:,:,planeIdxs(i))); 
    axis off; ax = gca;
    ax.Position = ax.OuterPosition;
    if ~sum(sum(planeWave(:,:,planeIdxs(i))))
        ax.CLim = [0 1];
    end
    colormap(cmap)
    saveas(gcf, strcat(planeFilename, num2str(planeIdxs(i)), '.png'));
    close all

    %%%% pup wave
    figure('Renderer', 'painters', 'Position', [0 0 1000 1000],'color', 'w');
    imagesc(pupWave(:,:,pupIdxs(i)));
    axis off; ax = gca;
    scale = [min(min(min(pupWave))), max(max(max(pupWave)))];
    ax.Position = ax.OuterPosition;
    ax.CLim = scale;
    colormap(cmap)
    saveas(gcf, strcat(pupFilename, num2str(pupIdxs(i)), '.png'));
    close all

end

%---- FLOW Portraits ----%
%%%% traveling wave
filename = './figure5/planeWave_';
int_len = 15;
flowPortraitR2020b(planeWave(:,:,10:end), int_len, 'history_delay',0, 'filename',filename,'img_format','.png', 'save_flow',false, 'thresh_quantile', 0.91);

%%%% pup wave
filename = './figure5/pupWave_';
int_len = 5;
hist_delay = 10;
flowPortraitR2020b(pupWave, int_len, 'history_delay',hist_delay,'filename',filename,'img_format','.png', 'save_flow',false, 'thresh_quantile', 0.90,...
'figsize', [0 0 500 500]);