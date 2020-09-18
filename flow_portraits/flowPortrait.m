% AUTHOR:
%   Nathaniel Linden
%
% FUNCTION:
%   flowPortrait
%
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
%   flowPortrait(data, integration_length, 'Optional Name-Value Pair Arguments')

function flowPortrait(data, integrationLen, varargin)
  
    % default values
    default_alpha_hs=       0.75;
    default_max_iter=       100;
    default_avg_win=        5;
    default_smooth_win=     5;
    default_alpha_smooth=   1.25;
    defaultHistoryDelay=    30;
    defaultTreshQuantile=   0.93;
    defaultForColor=        [1 0.647058823529412 0];
    defaultBackColor=       [0.525490196078431 0 0.831372549019608];
    defaultBackgroundCmap=  colormap('gray');
    defaultFigSize=         [0 0 1000 500];
    defaultFilename=        './flow_image.png';
    defaultMask=            [];
    
    % input parsing
    validColor = @(x) numel(x)==3 && sum(x)<3 && ndims(x)==2;
    validCMAP = @(x) isstring(x) || isnumeric(x);
    validFigSize = @(x) isnumeric(x) && ndims(x)==4;
    validMask = @(x) (size(x,1)==size(data,1) && size(x,2)==size(data,2)) || x==[];

    p = inputParser;  % input parser for optional parameters
    addOptional(p, 'alpha_hs',          default_alpha_hs,      @isnumeric);
    addOptional(p, 'max_iter',          default_max_iter,      @isnumeric);
    addOptional(p, 'avg_win',           default_avg_win,       @isnumeric);
    addOptional(p, 'smooth_win',        default_smooth_win,    @isnumeric);
    addOptional(p, 'alpha_smooth',      default_alpha_smooth,  @isnumeric);
    addOptional(p, 'history_delay',     defaultHistoryDelay,   @isnumeric);
    addOptional(p, 'thresh_quantile',   defaultTreshQuantile,  @isnumeric);
    addOptional(p, 'for_color',         defaultForColor,       validColor);
    addOptional(p, 'back_color',        defaultBackColor,      validColor);
    addOptional(p, 'background_cmap',   defaultBackgroundCmap, validCMAP);
    addOptional(p, 'figsize',           defaultFigSize,        validFigSize);
    addOptional(p, 'filename',          defaultFilename,       @isstring);
    addOptional(p, 'start_frame',       1,                     @isnumeric);
    addOptional(p, 'save_optical_flow', false,                 @islogical);
    addOptional(p, 'save_ftle',         false,                 @islogical);
    addOptional(p, 'img_format',        '.png',                @isstring);
    addOptional(p, 'white_background',  false,                 @islogical); 
    addOptional(p, 'mask',              defaultMask,           validMask);
    

    parse(p);
    alpha_hs =          p.Results.alpha_hs;
    max_iter =          p.Results.max_iter;
    avg_win =           p.Results.avg_win;
    smooth_win =        p.Results.smooth_win;
    alpha_smooth =      p.Results.alpha_smooth;
    history_delay =     p.Results.history_delay;
    thresh_quantile =   p.Results.thresh_quantile;
    for_color =         p.Results.for_color;
    back_color =        p.Results.back_color;
    background_cmap =   p.Results.background_cmap;
    figsize =           p.Results.figsize;
    filename =          p.Results.filename;
    save_optical_flow = p.Results.save_optical_flow;
    save_ftle =         p.Results.save_ftle;
    img_format =        p.Results.img_format;
    start_frame =       p.Results.start_frame;
    white_background =  p.Results.white_background;
    mask =              p.Results.mask;

    who

    % preprocessing
    dims = size(data);
    mean_data = mean(data,3);
    disp('Computing Optical Flow...');

    % compute Optical Flow using the Horn-Schunck (hs) method
    [x, y, u, v, ~] = smoothScaledHornSchunk(data, history_delay, alpha_hs, ...
                                        max_iter, avg_win, smooth_win, alpha_smooth);

    if save_optical_flow  % save optical flow output if the user specified
        optical_flow.x =     x; 
        optical_flow.y =     y; 
        optical_flow.u =     u; 
        optical_flow.v =     v;
        optical_filename =   strcat(filename,'_optical_flow.mat');
        save(optical_filename, 'optical_flow', '-v7.3');
    end
    clear data; % clear raw data to save space in RAM
    
    disp('Optical Flow Complete. Computing FTLE...');

    % compute Finite Time Lyapunov Exponent from the vector field 
    resolution = [dims(2), dims(1)];
    [ftle] = FTLECompute(start_frame, integration_len, resolution, u, v);

    if save_ftle  % save optical flow output if the user specified
        ftle_filename = strcat(filename,'_ftle.mat');
        save(ftle_filename, 'optical_flow', '-v7.3');
    end

    disp('FTLE Complete. Computing and Saving FLOW portrait...');

    % compute FLOW portrait from the ftle
    [forward_flow, backward_flow] = createFLOWPortrait(mean(ftle.f,3), ...
                                                    mean(ftle.b,3), thresh_quantile);

    % save output image
    flow_filename = strcat(filename, '_flow_portrait', img_format);
    saveFlowPortrait(forward_flow, backward_flow, mean_data, flow_filename, white_background, ...
                        mask, for_color, back_color, background_cmap, figsize);
end