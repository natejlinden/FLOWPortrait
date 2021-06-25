% ORIGINAL CODE FROM Linden et al. 2020 (Go with the FLOW: Visualizing spatiotemporal dynamics in1optical widefield calcium imaging)

% This script recreates all the panes in Figure 5
 
% NOTE: The colormaps here use system colormaps and differ from those in than the manuscript. See colormaps section below for details.

clear all; close all; clc

%%%% COLORMAPS %%%%
%   This script uses MATLAB colormaps by defualt, however these were not used in the manuscript.
cmap_calc = colormap('winter');
spec = colormap('parula');
l1 = colormap('gray');

% To exactly recreate the images in the manuscript please download and add ColorBrewer and Colorcet to the path:
%   https://www.mathworks.com/matlabcentral/fileexchange/45208-colorbrewer-attractive-and-distinctive-colormaps
%   https://peterkovesi.com/projects/colourmaps/
%   
%   AND uncomment below
%
% cmap_calc = flip(brewermap(64,'rdbu'),1);
% spec = brewermap(64,'prgn');
% l1 = colorcet('L1');
spec2 = [ones(32,3); spec(33:end,:)];
forward_color = [255 165 0]/255; backward_color = [134,0,212]/255;
figsize = [0 0 348 348];

%%%% ADD REQUIRED CODE %%%% 
addpath('../flow_portraits/')
mkdir('./figure4/')

% load data
load('../data/figure7_ftle_data/masks/SamplePupTrace_14_mask.mat');
Mask = discard_zero_cols(Mask);
load('../data/figure4_ftle_data.mat')
ftle = data.ftle;
mean_dFoF = mean(data.dFof,3);

clearvars data 

% compute mean FTLE > 0
idx = ftle.f > 0;
thresh_3d = ftle.f .* idx;
mean_forward = mean(thresh_3d,3);
idx = ftle.r > 0;
thresh_3d = ftle.r .* idx;
mean_reverse = mean(thresh_3d,3);

% create FLOW portrait w/ morphological image processing
[for_FLOW, rev_FLOW, for_skel, rev_skel] = createFLOWPortrait(mean_forward,mean_reverse, 0.95);

% save mean ftle images
filename_mean_for = './figure4/A_for_mean.png';
filename_mean_rev = './figure4/A_backward_mean.png';


% saveFlowPortrait(zeros(size(mean_forward)), zeros(size(mean_forward)),ind2rgb(mean_forward, spec2), filename_mean_for, true, Mask, forward_color, backward_color, spec2, figsize)
% saveFlowPortrait(zeros(size(mean_reverse)), zeros(size(mean_reverse)),ind2rgb(mean_reverse, spec2), filename_mean_rev, true, Mask, forward_color, backward_color, spec2, figsize)

saveFlowPortrait(zeros(size(mean_forward)), zeros(size(mean_forward)),mean_forward, filename_mean_for, true, Mask, forward_color, backward_color, spec2, figsize)
saveFlowPortrait(zeros(size(mean_reverse)), zeros(size(mean_reverse)),mean_reverse, filename_mean_rev, true, Mask, forward_color, backward_color, spec2, figsize)

figure('Renderer', 'painters', 'Position', figsize,'color','w');
scale = [min(min(mean_forward)), max(max(mean_forward))];

newImg = mean_forward; % imadjust(mean_forward, [0 max(max(mean_forward))]);
[newImg, imAlpha] = whiteBackground(newImg, Mask);

a = imagesc(newImg,[-1*max(max(mean_forward)) max(max(mean_forward))]);
a.AlphaData = imAlpha;

ax = gca; 
axis off; axis image 
ax.Position = ax.OuterPosition;
colormap(spec)
saveas(gcf, filename_mean_for);
close all

figure('Renderer', 'painters', 'Position', figsize,'color','w');
scale = [min(min(mean_reverse)), max(max(mean_reverse))];

newImg = mean_reverse; % imadjust(mean_forward, [0 max(max(mean_forward))]);
[newImg, imAlpha] = whiteBackground(newImg, Mask);

a = imagesc(newImg,[-1*max(max(mean_reverse)) max(max(mean_reverse))]);
a.AlphaData = imAlpha;

ax = gca; 
axis off; axis image 
ax.Position = ax.OuterPosition;
colormap(spec)
saveas(gcf, filename_mean_rev);
close all
%%
% Standard Colors
forward_color = [255 166 0]/255;
backward_color = [134 0 212]/255;
figsize = [0 0 348 348];

filename_processed = './figure4/D_flow_port.png';
saveFlowPortrait(for_FLOW, rev_FLOW, mean_dFoF, filename_processed, true, Mask, forward_color, backward_color, l1, figsize)


filename_skel = './figure4/C_skel.png';
saveFlowPortrait(for_skel, rev_skel, mean_dFoF, filename_skel, true, Mask, forward_color, backward_color, l1, figsize)


 % perform thresholding of mean FTLE (if < thresh set to 0; else set to 1)
 thresh_for = mean_forward; thresh_rev = mean_reverse;
 % forward
 for_thresh_val = quantile(thresh_for(:),0.95);
 thresh_for(thresh_for < for_thresh_val) = 0;
 thresh_for(thresh_for >= for_thresh_val) = 1;
 % reverse
 rev_thresh_val = quantile(thresh_rev(:),0.95);
 thresh_rev(thresh_rev < rev_thresh_val) = 0;
 thresh_rev(thresh_rev >= rev_thresh_val) = 1;

filename_thresh = './figure4/B_thresholded.png';
saveFlowPortrait(thresh_for, thresh_rev, mean_dFoF, filename_thresh, true, Mask, forward_color, backward_color, l1, figsize)

%%%% FUNCTIONS %%%%
function [modified_data, zero_idxs] = discard_zero_cols(data,discard_rows)
	% convert matrix to a logical array and take the sum of a row
	% a full row of zeros will come out as 0
	if nargin < 2
		discard_rows=false;
	end

	% discard cols
	zero_cols = sum(logical(data),1); 
	always_zero_cols = sum(zero_cols,3);  % add up in time to ensure only discarding constanly zero cols
	modified_data = data(:,find(always_zero_cols),:);  % extract non-zero columns and return

	% discard rows
	if discard_rows
		zero_rows = sum(logical(data),2);
		always_zero_rows = sum(zero_rows,3);  % add up in time to ensure only discarding constanly zero cols
		modified_data = modified_data(find(always_zero_rows),:,:);  % extract non-zero columns and return
	end
end