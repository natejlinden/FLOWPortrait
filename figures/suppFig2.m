% ORIGINAL CODE FROM Linden et al. 2020 (Go with the FLOW: Visualizing spatiotemporal dynamics in1optical widefield calcium imaging)

% This script recreates all the panes in Supplementary Figure 2
 
% NOTE: The colormaps here use system colormaps and differ from those in than the manuscript. See colormaps section below for details.

clear all; close all; clc

% Add required code
% Load dependendies
addpath('../flow_portraits/')
mkdir('./suppFig2/')
mkdir('./suppFig2/panes/')

grassGreen = [63, 155, 11]/255;

% To exactly recreate the images in the manuscript please download and add ColorBrewer and Colorcet to the path:
%   https://www.mathworks.com/matlabcentral/fileexchange/45208-colorbrewer-attractive-and-distinctive-colormaps
%   https://peterkovesi.com/projects/colourmaps/
%   
%   AND uncomment below
%
% 
% cmap = colorcet('L5');
cmap = colormap('gray');

%---- Load synthetic data ----%
% See: xxx.m for code to generate data
%%%% Traveling wave:
load('../data/PlaneWave.mat')
load('../data/planeWaveSourceSink.mat');

%%% Circular wave:
load('../data/CircularWave.mat')
load('../data/circularWaveSourceSink.mat');


%%%% Traveling Gaussian (No noise)
load('../data/TravelingGaussian.mat')
load('../data/TravelingGaussianSourceSink.mat');

%---- Compute FLOW portraits for each dataset ----%
%%%% traveling wave
filename = './suppFig2/FLOW_plane_wave';
int_len = 15;
flowPortraitR2020b(planeWave(:,:,10:end), int_len, 'history_delay',0, 'filename',filename,'img_format','.png', 'save_flow',false, 'thresh_quantile', 0.91);

%%%% circular wave
filename = './suppFig2/FLOW_circular_wave';
int_len = 12;
flowPortraitR2020b(circularWave, int_len, 'history_delay',0, 'filename',filename,'img_format','.png', 'save_flow',false);

%%%% Traveling Gaussian (No noise)
filename = './suppFig2/FLOW_taveling_gaussian';
int_len = 10;
flowPortraitR2020b(travelingGaussian, int_len, 'history_delay',0, 'filename',filename,'img_format','.png', 'save_flow',false, 'thresh_quantile',0.85);

%---- Save panes of data ----%
planeFilename = './suppFig2/panes/planeWave_';
circFilename = './suppFig2/panes/circularWave_';
gaussFilename = './suppFig2/panes/travelingGaussian_';


planeStride = 5; planeIdxs = 15:planeStride:40;
circularStride = 5; circularIdxs = 8:circularStride:size(circularWave,3);
gaussStride = 10; gaussIdxs = 0:gaussStride:size(travelingGaussian,3); gaussIdxs(1) = 1;
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

    %%%% circular wave
    figure('Renderer', 'painters', 'Position', [0 0 1000 1000],'color', 'w');
    imagesc(circularWave(:,:,circularIdxs(i)));
    axis off; ax = gca;
    ax.Position = ax.OuterPosition;
    if ~sum(sum(circularWave(:,:,circularIdxs(i))))
        ax.CLim = [0 1];
    end
    colormap(cmap)
    saveas(gcf, strcat(circFilename, num2str(circularIdxs(i)), '.png'));
    close all

    %%%% traveling gaussian
    figure('Renderer', 'painters', 'Position', [0 0 1000 1000],'color', 'w');
    imagesc(travelingGaussian(:,:,gaussIdxs(i)));
    axis off; ax = gca;
    ax.Position = ax.OuterPosition;
    if ~sum(sum(travelingGaussian(:,:,gaussIdxs(i))))
        ax.CLim = [0 1];
    end
    colormap(cmap)
    saveas(gcf, strcat(gaussFilename, num2str(gaussIdxs(i)), '.png'));
    close all
end

%---- overlay OFAMM source/sink ----%
gaussFilename = './suppFig2/sourcesink_travelingGaussian';
planeFilename = './suppFig2/sourcesink_planeWave';
circleFilename = './suppFig2/sourcesink_circleWave';

% gaussian
filename = strcat(gaussFilename, '__.png');
source_sink_options.source_color = grassGreen; source_sink_options.sink_color = 'cyan'; source_sink_options.source_marker = '.'; source_sink_options.sink_marker = '.';
plotSourcSink(travelingGaussian, gaussianSourceSink.source(1:2,:),gaussianSourceSink.sink(1:2,:), false, [], filename, false, source_sink_options)

% plane wave
filename = strcat(planeFilename, '__.png');
source_sink_options.source_color = grassGreen; source_sink_options.sink_color = 'cyan'; source_sink_options.source_marker = '.'; source_sink_options.sink_marker = '.';
plotSourcSink(planeWave, planeSourceSink.source, planeSourceSink.sink, false, [], filename, false, source_sink_options)

% circle wave
filename = strcat(circleFilename, '__.png');
source_sink_options.source_color = grassGreen; source_sink_options.sink_color = 'cyan'; source_sink_options.source_marker = '.'; source_sink_options.sink_marker = '.';
plotSourcSink(circularWave, circularSourceSink.source, circularSourceSink.sink, false, [], filename, false, source_sink_options)

%---- Function to plot with overlaid source/sink
function plotSourcSink(data, source, sink, white_background, mask, filename, adjust, source_sink_options)
    cmap = colormap('gray');
    source_color = source_sink_options.source_color;
    source_marker = source_sink_options.source_marker;
    sink_color = source_sink_options.sink_color;
    sink_marker = source_sink_options.sink_marker;
    figsize = [0 0 600 600];
    mean_dFoF = mean(data,3);

    % aggreate source/sink into binary image
    sourceBinary = zeros(size(mean_dFoF));
    sinkBinary = zeros(size(mean_dFoF));

    for i = 1:size(source,2)
        sourceBinary(source(2,i), source(1,i))=1;
    end
    for i = 1:size(sink,2)
        sinkBinary(sink(2,i), sink(1,i))=1;
    end

    figure('Renderer', 'painters','position', figsize);
     % adjust range of backgorund image
     if adjust
        scale = [min(min(mean_dFoF)), max(max(mean_dFoF))];
        upper_val = 0.65; 
        if scale(2) < upper_val
            mean_dFoF = imadjust(mean_dFoF, scale, [0 upper_val]);
        else
            mean_dFoF = imadjust(mean_dFoF);
        end
    end
     
     % white_background
     if white_background
         [newImg, imAlpha] = whiteBackground(mean_dFoF, mask);
         mean_dFoF = newImg;
     else
         imAlpha = ones(size(mean_dFoF));
     end

     % plot background image
    img1 = imoverlay(mean_dFoF, sourceBinary, source_color);
    img2 = imoverlay(img1, sinkBinary, sink_color);
    imagesc(img2, 'AlphaData',imAlpha);
    hold on
    
    % sources
    plot(source(1,:), source(2,:), source_marker,'color',source_color,'MarkerSize', 60);
    % sinks
    plot(sink(1,:), sink(2,:), sink_marker,'color',sink_color,'MarkerSize', 60);

    % Format and Save Image
    ax = gca;
    if adjust
        ax.CLim = scale;
    end
    axis off; axis tight equal 
    colormap(cmap);
    ax.Position = ax.OuterPosition;
    saveas(gcf, filename);
end

