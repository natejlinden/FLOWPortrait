% ORIGINAL CODE FROM Linden et al. 2020 (Go with the FLOW: Visualizing spatiotemporal dynamics in1optical widefield calcium imaging)

% This script recreates all the panes in Figure 5

% NOTE: The colormaps here use system colormaps and differ from those in than the manuscript. See colormaps section below for details.

clear all; close all; clc

%%%% ADD REQUIRED CODE %%%%
addpath('../flow_portraits/')
mkdir('./figure5/')


%%%% COLORMAPS %%%%
%   This script uses MATLAB colormaps by defualt, however these were not used in the manuscript.
cmapl1 = colormap('gray');
cmapl5 = colormap('gray');

% To exactly recreate the images in the manuscript please download and add ColorBrewer and Colorcet to the path:
%   https://www.mathworks.com/matlabcentral/fileexchange/45208-colorbrewer-attractive-and-distinctive-colormaps
%   https://peterkovesi.com/projects/colourmaps/
%   
%   AND uncomment below
%
% cmapl1 = colorcet('L1');
% cmapl5 = colorcet('L5');


%%%% Two files corresponding the two waves in figure 5 %%%%
files = {'../data/developing_wave1.mat', '../data/developing_wave2.mat'};

%%%% PROCESS AND PLOT EACH WAVE %%%%
for wave = 1:2
    % load data
    load(files{wave});
    dfof = dat.dfof; mask = dat.mask;
    clear dat;

    % Create flow portrait to overlay
    integrationLength = 40;
    filename = strcat('./figure5/wave_', num2str(wave), '_flowPortrait.png');
    flowPortrait(dfof, integrationLength, 'filename', filename, 'mask', mask, 'background_cmap', cmapl1, 'white_background', true, 'thresh_quantile', 0.92);

    %% Save panes 
    fig_size = [0 0 180 173];

    temp = dfof(:,:,70:end-39);
    num_frames = 12;
    idxs = linspace(1, size(temp,3),num_frames);

    for j = 1:num_frames
        frm = idxs(j);
        figure('Renderer', 'painters', 'Position', fig_size);
        frame = temp(:,:,round(frm));

        [newImg, imAlpha] = whiteBackground(frame, mask);
        imagesc(newImg, 'AlphaData',imAlpha); 
        
        % plot formatting
        colormap(cmapl5)
        % gc = gcf; gc.Position = fig_size;
        axis off;
        
        ax = gca; 
        ax.Position = ax.OuterPosition;
        % ax1 = hi.Parent;
        % ax1.Position = ax.OuterPosition;
        filename = strcat('./figure5/wave_',num2str(wave),'_frame_',num2str(round(frm)),'.png');
        saveas(gcf, filename);
    
        disp(frm)
        close gcf

    end

    close all

end