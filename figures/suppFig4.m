% ORIGINAL CODE FROM Linden et al. 2020 (Go with the FLOW: Visualizing spatiotemporal dynamics in1optical widefield calcium imaging)

% This script recreates all the panes in Supplementary Figure 4
 
% NOTE: The colormaps here use system colormaps and differ from those in than the manuscript. See colormaps section below for details.

clear all; close all; clc

% Add required code
% Load dependendies
addpath('../flow_portraits/')
mkdir('./suppFig4/')


% load traveling Gaussian data
load('../data/TravelingGaussian.mat')

%% Panes with varrying integration length
% Compute FLOW portrait
filenameBase = './suppFig4/gaussian_flow_port_int_len';
int_lens = [5:5:35];
for i = 1:numel(int_lens)
    fname = strcat(filenameBase, '_',num2str(int_lens(i)))
    flowPortraitR2020b(travelingGaussian, int_lens(i), 'history_delay',0, 'filename',fname,'img_format','.png', 'thresh_quantile',0.85);
end

%% Panes with varrying threshold
% Compute FLOW portrait
filenameBase = './suppFig4/gaussian_flow_port_thresh';
thresholds = [0.75, 0.80, 0.85, 0.90, 0.95];
int_len = 10;
for i = 1:numel(thresholds)
    thresh = strsplit(num2str(thresholds(i)), '.');
    fname = [filenameBase, '_0_',thresh{2}];
    flowPortraitR2020b(travelingGaussian, int_len, 'history_delay',0, 'filename',fname,'img_format','.png', 'thresh_quantile',thresholds(i));
end