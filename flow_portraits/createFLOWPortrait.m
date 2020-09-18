% Nate Linden
% Brunton Lab
% Created: April 13, 2020
% Modified: June 18th, 2020

%% Thresholds and processes forward and backward FTLE to create the FLOW portrait
function varargout = createFLOWPortrait(mean_FTLE_for, mean_FTLE_back, thresh_quantile);
    %% NOTE: THIS FUNCTION REQUIRES THE MATLAB IMAGE PROCESSING TOOLBOX %%
    
    % Syntax
    % [for_FLOW, back_FLOW] = createFLOWPortrait(mean_FTLE_for, mean_FTLE_back);
    % [for_FLOW, back_FLOW] = createFLOWPortrait(mean_FTLE_for, mean_FTLE_back, thresh_quantile);
    % [for_FLOW, back_FLOW, for_skel, back_skel] = createFLOWPortrait(mean_FTLE_for,mean_FTLE_back);
    % [for_FLOW, back_FLOW, for_skel, back_skel] = createFLOWPortrait(mean_FTLE_for,mean_FTLE_back, thresh_quantile);
    
    % Inputs:
    %         mean_FTLE_for - NxM matrix with mean forward FTLE used for FLOW_portrait
    %         mean_FTLE_rev - NxM matrix with mean backward FTLE used for FLOW_portrait
    %         thresh_quantile - OPTIONAL. quantile for thresholding. Default to - 0.95
    % Returns:
    %          for_FLOW - NxM binary image of the forward time FLOW
    %          back_FLOW - NxM binary image of the backward time FLOW
    
        nargoutchk(1,4);
    
        % Defualt thresh_quantile
        if nargin < 3
            thresh_quantile = 0.95;
        end
    
        % perform thresholding of mean FTLE (if < thresh set to 0; else set to 1)
            thresh_for = mean_FTLE_for; thresh_back = mean_FTLE_back;
            % forward
            for_thresh_val = quantile(thresh_for(:),thresh_quantile);
            thresh_for(thresh_for < for_thresh_val) = 0;
            thresh_for(thresh_for >= for_thresh_val) = 1;
            % backward
            back_thresh_val = quantile(thresh_back(:),thresh_quantile);
            thresh_back(thresh_back < back_thresh_val) = 0;
            thresh_back(thresh_back >= back_thresh_val) = 1;
    
        % perform morphological image processing to get skeleton
        % structures of FTLE ridges
        % Operations: close -> thin -> skeletonize
            % forward
            closed = bwmorph(thresh_for, 'close', Inf);
            skel = bwmorph(closed,'thin',Inf);
            opened_for = bwareaopen(skel,4);
            for_skel = opened_for;
    
            % backward
            closed = bwmorph(thresh_back, 'close', Inf);
            skel = bwmorph(closed,'thin',Inf);
            opened_back = bwareaopen(skel,4);
            back_skel = opened_back;
    
        % morphological processing to smooth skeleton structures
        % Operations: diag -> spur -> close
        for_FLOW = bwmorph(bwmorph(bwmorph(opened_for,'diag'),'spur'),'close');
        back_FLOW = bwmorph(bwmorph(bwmorph(opened_back,'diag'),'spur'),'close');
    
        % assign outputs
        switch nargout
            case 2
                varargout{1}=for_FLOW;
                varargout{2}=back_FLOW;
            case 4
                varargout{1}=for_FLOW;
                varargout{2}=back_FLOW;
                varargout{3}=for_skel;
                varargout{4}=back_skel;
            otherwise
                error('Incorrect number of output arguments!');
        end
    end