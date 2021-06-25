% AUTHOR:
%   Nathaniel Linden
%
% FUNCTION:
%   flowPortraitR2020b
%%% FOR USE WITH MATLAB R2020b %%%%%%
%% DESCRIPTION: This function computes a FLOW portrait for a segment of widefield imaging data %%
%
% INPUT:
%   data (required) - (double) MxNxT input data matrix to compute FLOW portrait
%   integrationLen (required) - (double) integration length for FTLE computation, must not be longer than T
%   alpha_hs (optional) - (double) smoothness parameter for Optical Flow, default 0.75
%   max_iter (optional) - (double) maximum number of Optical Flow iterations, default 100
%   avg_win (optional) - (double) window length (# of frames) for optical flow windowed averaging, default 5
%   smooth_win (optional) - (double) kernel length (# of frames) for optical flow Gaussian smoothing, default 5
%   history_delay (optional) - (double) smoothness parameter for Gaussian Smoothing, default 1.25
%   tresh_quantile (optional) - (double) quantile for thresholding to create FLOW portrait, default 0.93
%   for_color (optional) - (double) [R G B] matrix for forward FLOW color, default [1 0.647058823529412 0] 
%   back_color (optional) - (double) [R G B] matrix for backward FLOW color, default [0.525490196078431 0 0.831372549019608] 
%   background_cmap (optional) - (double) colormap matrix for FLOW portraits, default gray
%   figsize (optional) - (double) [X Y W H] matrix for figure size, default [0 0 1000 500]
%   filename (optional) - (char) string for output file, default './flow_image.png'
%   start_frame (optional) - (double) index of starting frame for computation, default 1
%   save_optical_flow (optional) - (logical) boolean to indicate if user wishes to save optical flow output data
%   save_ftle (optional) - (logical) boolean to indicate if user wishes to save FTLE output data
%   img_format (optional) - (char) string indicating desired file extension, default '.png'
%   white_background (optional) - (logical) boolean to indicate is user would like to use a mask for making the non-cortical regions white, 
%                                  requires mask to be specified
%   mask (optional) - (logical) - MxN logical matrix to mask non-brain areas, defualt [] 
%
% OUTPUT:
%   none
%
% ASSUMPITONS and LIMITATIONS:
%   none
%
% SYNTAX:
%   flowPortrait(data, integration_len, 'Optional Name-Value Pair Arguments')

function flowPortraitR2020b(data, integration_len, options)
    arguments
        data (:,:,:) double
        integration_len (1,1) double
        % optional name value args
        options.alpha_hs (1,1) double = 0.75;
        options.max_iter (1,1) double = 100;
        options.avg_win (1,1) double = 5;
        options.smooth_win (1,1) double = 5;
        options.alpha_smooth (1,1) double = 1.25;
        options.history_delay (1,1) double = 30;
        options.thresh_quantile (1,1) double = 0.93;
        options.for_color (1,3) double {validColor(options.for_color)} = [1 0.647058823529412 0];
        options.back_color (1,3) double {validColor(options.back_color)} = [0.525490196078431 0 0.831372549019608];
        options.background_cmap (:,:) double = colormap('gray');
        options.figsize (1,4) double = [0 0 1000 500];
        options.filename char = './output';
        options.save_optical_flow (1,1) logical = false;
        options.save_ftle (1,1) logical = false;
        options.save_flow (1,1) logical = true;
        options.start_frame (1,1) double = 1;
        options.img_format char = '.png';
        options.white_background (1,1) logical = false;
        options.mask (:,:) = [];
    end

    % preprocessing
    dims = size(data);
    mean_data = mean(data,3);
    disp('Computing Optical Flow...');

    % compute Optical Flow using the Horn-Schunck (hs) method
    [x, y, u, v, ~] = smoothScaledHornSchunck(data, options.history_delay, options.alpha_hs, ...
        options.max_iter, options.avg_win, options.smooth_win, options.alpha_smooth);

    if options.save_optical_flow  % save optical flow output if the user specified
        optical_flow.x =     x; 
        optical_flow.y =     y; 
        optical_flow.u =     u; 
        optical_flow.v =     v;
        optical_filename =   strcat(options.filename,'_optical_flow.mat');
        save(optical_filename, 'optical_flow', '-v7.3');
    end
    clear data; % clear raw data to save space in RAM
    
    disp('Optical Flow Complete. Computing FTLE...');

    % compute Finite Time Lyapunov Exponent from the vector field 
    resolution = [dims(2), dims(1)];
    [ftle] = FTLECompute(options.start_frame, integration_len, resolution, u, v);

    if options.save_ftle  % save optical flow output if the user specified
        ftle_filename = strcat(options.filename,'_ftle.mat');
        save(ftle_filename, 'ftle', '-v7.3');
    end

    disp('FTLE Complete. Computing and Saving FLOW portrait...');

    % remove negative FTLE values
    idx = ftle.f > 0;
    ftle.ff = ftle.f .* idx;
    idx = ftle.b > 0;
    ftle.r = ftle.b .* idx;

    % compute FLOW portrait from the ftle
    [forward_flow, backward_flow] = createFLOWPortrait(mean(ftle.f,3), ...
        mean(ftle.b,3), options.thresh_quantile);

    % save raw flow data
    if options.save_flow
        flOW_filename = strcat(options.filename, '_raw_FLOW.mat')
        flOW.forward = forward_flow; flOW.backward = backward_flow;
        save(flOW_filename,'flOW', '-v7.3');
        clear flOW
    end

    % save output image
    flow_filename = strcat(options.filename, options.img_format);
    saveFlowPortrait(forward_flow, backward_flow, mean_data, flow_filename, options.white_background, ...
        options.mask, options.for_color, options.back_color, options.background_cmap, options.figsize);

end

% additional functions to validate inputs
function validColor(a)
    % Test for color propertiesw
    if sum(a)>3
        eid = 'Color:notValid';
        msg = 'Color values must be between 0 and 1.';
        throwAsCaller(MException(eid,msg))
    end
end
