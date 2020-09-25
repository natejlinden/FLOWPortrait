% ORIGINAL CODE FROM Linden et al. 2020 (Go with the FLOW: Visualizing spatiotemporal dynamics in1optical widefield calcium imaging)

% This script recreates all the panes in Figure 3 and Supplementary Figure 1
 
% NOTE: The colormaps here use system colormaps and differ from those in than the manuscript. See colormaps section below for details.

% clear all; close all; clc

%%%% ADD REQUIRED CODE %%%%
addpath('../flow_portraits/')
mkdir('./figure3_sup1/')


%%%% COLORMAPS %%%%
%   This script uses MATLAB colormaps by defualt, however these were not used in the manuscript.
cmap_calc = colormap('winter')
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

%%%% CONSTANTS %%%%
num_clips = 5; step_size = 10;
times = {'0','0.5','1','1.5','2'};

%%%% LOAD and Preprocess DATA %%%%
% load('../data/fig3_developing.mat')
% dfof = dat.dfof; mask = dat.mask;
% clear dat;
% dims = size(dfof);
    
% % compute optical flow
% [x,y,u,v,subDat] = smoothScaledHornSchunck(dfof,30); 

% % compute FTLE
% integrationLength = 40;
% % ftle = FTLECompute(1, integrationLength, [dims(2), dims(1)], u, v);  

% datShow = subDat(:,:,180:220);
% [datShow, imAlpha] = whiteBackground(datShow, mask);
% uShow = u(:,:,180:220);
% vShow = v(:,:,180:220); 
% mean_dFoF = mean(datShow,3);

% %%%% FIGURE 3 %%%%

% % Save panes of raw dfof data (Figure 3A)
% base_file_name = './figure3_sup1/A_';
% idx = 1; 

% for clip = 1:num_clips
%     % create square figure with white background
%     figure('Renderer', 'painters', 'Position', [0 0 900 900],'color', 'w');
    
%     % show image
%     imagesc(datShow(:,:,idx), 'AlphaData', imAlpha); hold on
%     rectangle('Position', [40 60 60 60],'EdgeColor','yellow','LineWidth',10)
    
%     % plot formatting
%     colormap(l1); 
%     axis off; axis equal
%     ax = gca; ax.Position = [-0.04 -0.01 1.08 1.02];
%     ax.Position;
%     fig = gcf; fig.Position = fig.OuterPosition;
%     fig.Position(3) = fig.Position(3)-4;
    
%     % save
%     filename = strcat(base_file_name,times{clip},'.png');
%     saveas(gcf, filename);
%     close all
   
%     % update idx
%     idx = idx + step_size;
% end


% %% Save panes of quiver plots (Figure 3B)
% base_file_name = './figure3_sup1/B_';
% idx = 1;
% for clip = 1:num_clips
%     % create square figure with white background
%     figure('Renderer', 'painters', 'Position', [0 0 900 900],'color', 'w');
    
%     % show image
%     imagesc(datShow(60:120,40:100,idx)); hold on; 
%     quiver(x(60:6:120,40:6:100)-39, y(60:6:120,40:6:100)-59, ...
%         uShow(60:6:120,40:6:100,idx), vShow(60:6:120,40:6:100,idx),1,'LineWidth',9,'MaxHeadSize',1,'color','yellow','AlignVertexCenters','on');
    
%     % plot formatting
%     colormap(l1);
%     axis off; axis equal
%     xlim([1 59]); ylim([1 59]);
%     ax = gca; ax.Position = [-0.04 -0.01 1.08 1.02];
%     ax.Position;
%     fig = gcf; fig.Position = fig.OuterPosition;
%     fig.Position(3) = fig.Position(3)-4;
    
%     % save
%     filename = strcat(base_file_name,times{clip},'.png');
%     saveas(gcf, filename);
%     close all
%     % update frames
%     idx = idx + step_size;
% end


%% Save panes of Forward FTLE (Figure 3C)
% also supplementary 1D
base_file_name = './figure3_sup1/C_';

idx=180;
for clip = 1:num_clips
    % create square figure with white background
    figure('Renderer', 'painters', 'Position', [0 0 900 900],'color', 'w');
    
    % show image
    dat = ftle.f(:,:,idx);
    normalize = [-1*max(max(abs(dat))), max(max(abs(dat)))];
    imagesc(dat,normalize); hold on
    
    % plot formatting
    colormap(spec);
    axis off; axis equal
    ax = gca; ax.Position = [-0.04 -0.01 1.08 1.02];
    ax.Position;
    fig = gcf; fig.Position = fig.OuterPosition;
    fig.Position(3) = fig.Position(3)-4;
    
    % save
    filename = strcat(base_file_name, times{clip}, '.png');
    saveas(gcf, filename);
    close all
    
    % update idx
    idx = idx + step_size;
end

%% Save panes of Backward FTLE (Figure 3D)
% also supplemntary 1E
base_file_name = './figure3_sup1/D_';
idx=180;
for clip = 1:num_clips
    % create square figure with white background
    figure('Renderer', 'painters', 'Position', [0 0 900 900],'color', 'w');
    
    % show image
    dat = ftle.b(:,:,idx);
    normalize = [-1*max(max(abs(dat))), max(max(abs(dat)))];
    imagesc(dat,normalize); hold on
    
    % plot formatting
    colormap(spec);
    axis off; axis equal
    ax = gca; ax.Position = [-0.04 -0.01 1.08 1.02];
    ax.Position;
    fig = gcf; fig.Position = fig.OuterPosition;
    fig.Position(3) = fig.Position(3)-4;
    
    % save
    filename = strcat(base_file_name, times{clip}, '.png');
    saveas(gcf, filename);
    close all
    
    % update idx
    idx = idx + step_size;
end


%% Save FLOW Portrait (Figure 3E)

base_file_name = './figure3_sup1/E_';
i = 180;
for clip = 1:num_clips
    figure('Renderer', 'painters', 'Position', [0 0 900 900],'color', 'w');
    idx = ftle.f(:,:,i) > 0;
    thresh_3d = ftle.f(:,:,i) .* idx;
    mean_for = mean(thresh_3d,3);
    
    idx = ftle.b(:,:, i) > 0;
    thresh_3d = ftle.b(:,:,i) .* idx;
    mean_rev = mean(thresh_3d,3);

    [for_FLOW, rev_FLOW] = createFLOWPortrait(mean_for, mean_rev, 0.92);

    forward_color = [255 166 0]/255; reverse_color = [134 0 212]/255;
    fig_size = [0 0 900 900];

    filename = strcat(base_file_name,times{clip},'.png');
    saveFlowPortrait(for_FLOW, rev_FLOW, mean_dFoF, filename, true, mask, forward_color, ...
                                    reverse_color, l1, fig_size);

    % update frames
    i = i + step_size;
end


% %%%% SUPPLEMENTARY FIGURE 1 %%%%

% % compute div and curl
% div = zeros(size(uShow));
% cur = zeros(size(uShow));
% for i = 1:size(uShow,3)
%     div(:,:,i)=divergence(uShow(:,:,i), vShow(:,:,i));
%     cur(:,:,i)=curl(uShow(:,:,i), vShow(:,:,i));
% end
% [div imAlpha] = whiteBackground(div,mask);
% [cur imAlpha] = whiteBackground(cur,mask);

% % Save panes of raw dfof data (Supplementary Figure 1A)
% base_file_name = './figure3_sup1/supA_';
% idx = 1; 

% for clip = 1:num_clips
%     % create square figure with white background
%     figure('Renderer', 'painters', 'Position', [0 0 900 900],'color', 'w');
    
%     % show image
%     imagesc(datShow(:,:,idx), 'AlphaData', imAlpha); hold on
    
%     % plot formatting
%     colormap(l1); 
%     axis off; axis equal
%     ax = gca; ax.Position = [-0.04 -0.01 1.08 1.02];
%     ax.Position;
%     fig = gcf; fig.Position = fig.OuterPosition;
%     fig.Position(3) = fig.Position(3)-4;
    
%     % save
%     filename = strcat(base_file_name,times{clip},'.png');
%     saveas(gcf, filename);
%     close all
   
%     % update idx
%     idx = idx + step_size;
% end

% %% Save panes of Divergence (Supplementary Figure 1B)
% base_file_name = './figure3_sup1/supB_';
% current_frame = 1;
% idx =1;
% for clip = 1:num_clips
%     % create square figure with white background
%     figure('Renderer', 'painters', 'Position', [0 0 900 900],'color', 'w');
    
%     % show image
%     dat = div(:,:,idx);
%     normalize = [-1*max(max(abs(dat))), max(max(abs(dat)))];
%     imagesc(dat,'AlphaData',imAlpha); hold on
    
%     % plot formatting
%     colormap(cmap_calc);
%     axis off; axis equal
%     ax = gca; 
%     ax.CLim = normalize; 
%     ax.Position = [-0.04 -0.01 1.08 1.02];
%     ax.Position;
%     fig = gcf; fig.Position = fig.OuterPosition;
%     fig.Position(3) = fig.Position(3)-4;
    
%     % save
%     filename = strcat(base_file_name, times{clip}, '.png');
%     saveas(gcf, filename);
%     close all
    
%     % update idx
%     idx = idx + step_size;
% end

%% Save panes of Curl (Supplementary Figure 1C)
base_file_name = './figure3_sup1/supC_';
current_frame = 1;
idx=1;
for clip = 1:num_clips
    % create square figure with white background
    figure('Renderer', 'painters', 'Position', [0 0 900 900],'color', 'w');
    
    % show image
    dat = cur(:,:,clip);
    normalize = [-1*max(max(abs(dat))), max(max(abs(dat)))];
    imagesc(dat,'AlphaData',imAlpha); hold on
    
    % plot formatting
    colormap(cmap_calc);
    axis off; axis equal
    ax = gca;
    ax.CLim = normalize; 
     ax.Position = [-0.04 -0.01 1.08 1.02];
    ax.Position;
    fig = gcf; fig.Position = fig.OuterPosition;
    fig.Position(3) = fig.Position(3)-4;
    
    % save
    filename = strcat(base_file_name, times{clip}, '.png');
    saveas(gcf, filename);
    close all
    
    % update idx
    idx = idx + step_size;
end