% Nathaniel Linden
% Created: July 5th 2019
% Modified: June 18th 2019
% Brunton Lab

% Function to run Horn-Schunk computation on a data stack and perform time-history
% scaling and temporal filtering operations of the optical flow fields
function [x, y, u, v, sub_data] = smoothScaledHornSchunck(data, history_delay, alpha_hs, max_iter, avg_win, smooth_win, alpha_smooth)
% Syntax:
% [x, y, u, sub_data] = smoothScaledHornSchunk(data, alpha_hs, max_iter, history_delay, avg_win, smooth_win, alpha_smooth)

% Inputs:
% data -          [num_rows, num_cols, time] matrix with data
% history_delay - integer, smoothing timescale for temporal smoothing
% alpha_hs -      double, Horn-Schunck spatial smoothing penalty weight
% max_iter -      int, maximum # of iterations for Horn-Scunck
% avg_win -       int, window length for running avg (should be quite short)
% smooth_win -    int, window length for gaussian smoothing kernel
% alpha_smooth -  double, parameter for gaussian smoothing kernel

% NOTE: Defaults for alpha_hs, max_iter, avg_win, smooth_win, alpha_smooth should be sufficient for FLOW portraits

% Outputs:
% x, y -         meshgrid position matrices with dimensions [num_rows, num_cols]
% u, v -         smoothed and scaled x and y velocity fields, respectively. Same size as data
% sub_data -     shortened data matrix 

    % default params
    if nargin < 3; warning('All default paramters used.'); alpha_hs = 0.75; end
    if nargin < 4; max_iter = 100; end
    if nargin < 5; avg_win = 5; end
    if nargin < 6; smooth_win = 5; end
    if nargin < 7; alpha_smooth = 1.25; end

    % create wait wait_bar
    wait_bar = waitbar(0, 'Intializing HS Method. Please wait.');

    % RUN HS code
    [x, y, u, v wait_bar] = opticalFlowHS(data, alpha_hs, max_iter, wait_bar);

    % Run scaling code
    waitbar(1, wait_bar, 'Runing Scaling Code.')
    weights = timeHistoryScale(data, history_delay, avg_win);

    u_w = u(:,:,history_delay-1:end) .* weights; 
    v_w = v(:,:,history_delay-1:end) .* weights;

    % run temporal smoothing Code
    waitbar(1, wait_bar, 'Runing Smoothing Code.')
    [u, v] = vectorTimeSmoothing(u_w, v_w, smooth_win, alpha_smooth);

    % make sub-data as subset of data correponding to smoothed_scaled VF
    sub_data = data(:,:,history_delay:end);

    % update user and close waitbar
    waitbar(1, wait_bar, 'HS Computation Complete.'); pause(1);
    
    close(wait_bar);
end


% Function to generate a weight matrix based on time history in inpout data
function [weight_matrix] = timeHistoryScale(data, history_delay, window)
    if nargin < 3
        window = 1;
            warning('No smoothing of weights will take place. Window set to 1');
    end

    % take abs difference of current time and time_history in the past
    hist_diff = abs(data(:,:,history_delay:end) - data(:,:,1:end-(history_delay-1)));

    % find the time maximum differencer for every pixel in space
    norm_factor = max(hist_diff,[],3);

    % normalize weights and force between 0 and 1
    weights = hist_diff ./ norm_factor;
    weights(weights > 1) = 1; % all values greater than 1 become 1

    % set  NaN valued weights to 0
    weights(isnan(weights)) = 0;

    % perform windowed average smoothing in time and return the weight matrix for vector scaling
    weight_matrix = movmean(weights, window, 3);
end


% Function to perform gaussian windowed filtering in time on  spatiotemporal
% optical flow field data
function [u_smooth, v_smooth] = vectorTimeSmoothing(u, v, window_length, alpha)
    % reshape vector components to enable convolution in time
    u_s = reshape(u, size(u,1)*size(u,2), size(u,3)); v_s = reshape(v, size(v,1)*size(v,2), size(v,3));

    % creates gaussian filtering kernel
    g = gausswin(window_length, alpha);

    % perform filtering convolution
    u_filt = conv2(u_s, g', 'same'); v_filt = conv2(v_s, g', 'same');

    % reshape and return filtered data
    u_smooth = reshape(u_filt, size(u)); v_smooth = reshape(v_filt, size(v));
end


