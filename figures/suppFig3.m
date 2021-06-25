% ORIGINAL CODE FROM Linden et al. 2021 (Go with the FLOW: Visualizing spatiotemporal dynamics in1optical widefield calcium imaging)
%
% This script recreates FLOW portrait and Source Sink panes in Supplemental Figure 3 (Panel B & C)
% Panel A is the same as Figure 1
% 
% NOTE: The colormaps here use system colormaps and differ from those in than the manuscript. See colormaps section below for details.

% add required code
addpath('../flow_portraits/');
mkdir('./suppFig3/')

load('../data/developing_wave2.mat');
load('../data/pupSourceSink.mat');
mean_dfof = mean(dat.dfof, 3); mask = dat.mask;
dFoF = dat.dfof;


% To exactly recreate the images in the manuscript please download and add ColorBrewer and Colorcet to the path:
%   https://www.mathworks.com/matlabcentral/fileexchange/45208-colorbrewer-attractive-and-distinctive-colormaps
%   https://peterkovesi.com/projects/colourmaps/
%   
%   AND uncomment below
%
% cmap = colorcet('l5');
cmap = colormap('gray');
grassGreen = [63, 155, 11]/255;

%---- FLOW portraits of each dataset ----%
%%%% Pup Data
pupFilename = './suppFig3/FLOW';
int_len = 20;
flowPortraitR2020b(dFoF, int_len,'filename',pupFilename,'img_format','.png', 'save_flow',false, 'white_background',true,...
    'mask', mask, 'save_flow', true, 'thresh_quantile', 0.91, 'figsize', [0 0 600 600]);

%---- Overlay source/sinks ----%
pupFilename = './suppFig3/sourcesink';

pupSource = sourceSink.source; pupSource = pupSource(1:2,pupSource(3,:)>40&pupSource(3,:)<200);
pupSink = sourceSink.sink; pupSink = pupSink(1:2,pupSink(3,:)>40&pupSink(3,:)<200);

filename = strcat(pupFilename, '__.png');
source_sink_options.source_color = grassGreen; source_sink_options.sink_color = 'cyan'; source_sink_options.source_marker = '.'; source_sink_options.sink_marker = '.';
plotSourcSink(dFoF, pupSource, pupSink, true, mask, filename, true, source_sink_options)

close all

%---- Function to plot FLOW portrait with overlaid source/sink
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
