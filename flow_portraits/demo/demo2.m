% Nathaniel Linden
% Bing Brunton Lab - University of Washington
% June 29th, 202

% This script provides a tutorial for choosing the correct integration length for the FTLE computation
% The data is provided by Professor William Moody (Biology, University of Washington)

% clear all; close all; clc
% addpath('../flow_portraits')

        % Preprocessing (Ignore me)
        % load data
        load('../../Data/DemoData.mat');
        dims = size(dat.dfof);
        % compute optical flow
        [x, y, u, v, ~] = smoothScaledHornSchunk(dat.dfof, 30);

%%%%% Example of choosing an FTLE integration length %%%%
% 1) To begin choose a range of integration values %
startval = 2; % the start length must be greater than 1
endval = 80;  % chose a relatively long end. We chose 80 since we know the ideal value is 40. But if we did not know thise choose a larger value

    % Parameter and variable intialization
    filename      = './example2';
    intLengths    = startval:endval;
    resolution    = [dims(2), dims(1)];
    ftle_forwards = zeros(dims(1), dims(2), numel(intLengths)); % empty matrix to store outputs
    ftle_backward = zeros(dims(1), dims(2), numel(intLengths)); % empty matrix to store outputs


% 2) Next loop over all integration lenghts and perform a single FTLE integration 
%       - perform a single int by setting the singleInt parameter to true in the FTLE_compute() func
for i = 1:numel(intLengths)
    [ftle] = FTLECompute(1, intLengths(i), resolution, u, v, true);
    ftle_forwards(:,:,i) = ftle.f(:,:,1);
    ftle_backward(:,:,i) = ftle.b(:,:,1);
    disp(num2str(intLengths(i)));
end

% 3) Visualize results in a video where the first FTLE frame is seen at increasing integration lengths
% The ideal integration length occurs when the FTLE image begins to "crispen". Eg when there are strong 
% details in the image. As the integration in increased past this the image becomes less clear.
% In the video below the titles are red when the integration length is not optimal, and becomes
% green when the integration length is ideal.

% Play a looping video of the different FTLEs
% ask user if they would like to save the video
saveVid = input('Press 1 if you would like to view and save a video (Only works on Windows/Mac). Or press 0 if you would only like to view it.');
% change img for white background
[ftle_forwards, imAlpha] = whiteBackground(ftle_forwards, dat.mask);
[ftle_backward, imAlpha] = whiteBackground(ftle_backward, dat.mask);

figure('Renderer', 'painters', 'Position', [0 0 800 400],'color', 'w');

if saveVid
    vidObj = VideoWriter('./integration_length_video', 'MPEG-4');
    vidObj.FrameRate = 15;
    open(vidObj);
end

for i = 1:size(ftle_forwards,3)
    if i > 35 && i < 50
        col = 'g';
    else
        col = 'r';
    end

    subplot(1,2,1);
    imagesc(ftle_forwards(:,:,i), 'AlphaData', imAlpha);
    colormap gray
    titl_ = strcat('Forward FTLE: Integration Length', {' '}, num2str(intLengths(i)));
    title(titl_, 'FontWeight', 'bold', 'Color', col);
    
    subplot(1,2,2);
    imagesc(ftle_backward(:,:,i), 'AlphaData', imAlpha);
    colormap gray
    titl_ = strcat('Backward FTLE: Integration Length', {' '}, num2str(intLengths(i)));
    title(titl_, 'FontWeight', 'bold', 'Color', col);

    if saveVid
        f = getframe(gcf);
        writeVideo(vidObj, f);
    end

    pause(0.05)
end

close(vidObj); close all
