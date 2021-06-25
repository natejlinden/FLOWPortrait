% ORIGINAL CODE FROM Linden et al. 2020 (Go with the FLOW: Visualizing spatiotemporal dynamics in1optical widefield calcium imaging)

% This script recreates all the panes in Figure 8
% The script does not perform the event extraction or plot traces of active area and movement scores as defined in the methods

% NOTE: The colormaps here use system colormaps and differ from those in than the manuscript. See colormaps section below for details.

clear all; close all; clc

%%%% ADD REQUIRED CODE %%%%
addpath('../flow_portraits/')
mkdir('./figure7/')

%%%% COLORMAPS %%%%
%   This script uses MATLAB colormaps by defualt, however these were not used in the manuscript.
l1 = colormap('gray');
l5 = colormap('gray');

% To exactly recreate the images in the manuscript please download and add ColorBrewer and Colorcet to the path:
%   https://www.mathworks.com/matlabcentral/fileexchange/45208-colorbrewer-attractive-and-distinctive-colormaps
%   https://peterkovesi.com/projects/colourmaps/
%   
%   AND uncomment below
%
% l1 = colorcet('L1');
% l5 = colorcet('L5');

%%%% CONSTANTS %%%%%
flow_size = [0 0 508 499];
pane_size = [0 0 300 300];
trace_size = [0 0 400 200]
files = {'../data/adult_event1.mat', '../data/adult_event2.mat'};
numframes=12;
pad = 140;

%%%%  LOAD AND PROCESS DATA %%%%
% loop over the 2 wave events
for idx = 2
    load(files{idx});
    dfof = dat.dfof; mask = dat.mask;
    relArea = dat.relArea; movement_score = dat.movement_score;
    bout_move = dat.bout_idxs;
    clear dat

    dims = size(dfof);

    % compute flow portrait
    % NOTE: Brain Atlas was overlaid manually using Adobe Illustrator.
    % meanWithOutlines.pdf has the aligned Atlas. Just replace mean background image with the FLOW portrait

    history_delay = 15; int_len = 20;
    filename = strcat('./figure7/adult_flow_event_', num2str(idx), '.png');
    flowPortrait(dfof(:,:,(pad+1):(end-pad)).*mask, int_len, 'filename', filename, 'history_delay', history_delay,...
                'white_background', true, 'mask', mask, 'figsize', flow_size);

    % save panes of data
    stride=floor(size(dfof,3)/numframes)
    frm=1;
    
    for j = 1:numframes
        figure('Renderer', 'painters', 'Position', pane_size);
        frame = dfof(:,:,frm).*mask;
        [newImg, imAlpha] = whiteBackground(frame, mask);
        imagesc(newImg,'AlphaData',imAlpha); 
        colormap(l5)
        axis off; axis image
        ax = gca; ax.Position = ax.OuterPosition;
        filename = strcat('./figure7/bout_',num2str(idx),'_frame_',num2str(j),'.png');
        saveas(gcf, filename);
        disp(frm)
        frm=frm+stride;
        close gcf
    end
end