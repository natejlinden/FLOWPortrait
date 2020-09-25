% Nathaniel Linden
% CREATED: July 2019
% MODIFIED: June 18th, 2020 
% Brunton Lab

% This code provides a function to compute an FTLE field from a vector field input

% The functions to generate a disctete flow map and compute eig_vals of the Cauchy Green
% strain tensor are modified from LCS-Tool-master by Haller et al. (ETH Zurich http://ETHZ.CH)

% DEPENDENCY:
%   The dependent funcitons are in LSC-Tool-master and are as follows:
% SYNTAX:
% [ftle_field] = FTLECompute(start_frame, integration_length, resolution, u, v)
% [ftle_field, cgStrainEigvals, cgStrain] = FTLECompute(start_frame, integration_length, resolution, u, v)

% INPUTS:
% start_frame -        the frame to begin computation
% integration_length - the number of frames used for FTLE computaiton
% resolution -         the number of grid points in the grid used in the format resolution = [num_cols, num_rows]
%                      where num_rows and num_cols specify the number of rows and columns of the data
% u, v -               stacks of velocity data with dimensions [num_rows, num_cols, time]

% OUTPUTS:
% ftle_field - struct with .f for forward stack of FTLE field and .b for backward
% cgStrain -   struct with .f for forward stack of Cauchy-Green strain field and .b for backward 
% det_cg -     struct with .f for forward stack of determinant/eigenvalue product of the cgStrain and .b for backward  


function varargout = FTLECompute(start_frame, integration_length, resolution, u, v, singleInt)
    if nargin < 6; singleInt=false; end
    
    nargoutchk(1,3);
    
     % add LCS-tool
     full_path = pwd; parts = strsplit(full_path, 'FLOWPortrait');
     path_to_LCS = strcat(parts{1}, 'FLOWPortrait/flow_portraits/LCS-tool/');
     addpath(path_to_LCS)

    % define computational domain as the entire frame
    domain  = [1, size(u, 2); 1, size(u, 1)];

    % compute FTLE
    [forward, backward] = FTLEInt(start_frame, u, v, integration_length, resolution, domain, singleInt);
    ftle_field.f = forward.FTLE; ftle_field.b = backward.FTLE;
    cgStrain.f = forward.cgStrain; cgStrain.b = backward.cgStrain;
    det_cg.f = forward.eigen_prod; det_cg.b = backward.eigen_prod;

    switch nargout
        case 1
            varargout{1}=ftle_field;
        case 2
            error('Number of output arguments incorrect.');
        case 3
            varargout{1}=ftle_field;
            varargout{2}=det_cg;
            varargout{3}=cgStrain;
        otherwise
            error('Number of output arguments incorrect.');
    end
end

function [forward, rever] = FTLEInt(start_frame, u_,v_,integration_length, resolution, domain, singleInt)

    % make equal Y resolution from X resolution
    res = resolution;

    % make position grids for griddedInterpolant
    x_ = 1:size(u_, 2);
    y_ = 1:size(u_, 1);
    t_ = 1:size(u_, 3);


    uInterpolant = griddedInterpolant({t_, y_, x_}, permute(u_, [3,1,2]));
    vInterpolant = griddedInterpolant({t_, y_, x_}, permute(v_, [3,1,2]));

    % define derrivative function
    lDeriv = @(t,x,~)derivative(t,x,uInterpolant,vInterpolant); % LCS-tool dependency

    forward.FTLE = []; forward.cgStrain = []; forward.eigen_prod = [];
    rever.FTLE = []; rever.cgStrain = []; rever.eigen_prod = [];

    % create matrices of integration indices to perform
    len = size(u_, 3);
    true_start = start_frame + integration_length-1;
    ints_forward = [[true_start:(len-integration_length)];
                    [(true_start+integration_length):len]]';
    ints_backward = [[true_start:(len-integration_length)];
                    [start_frame:((len+1)-(2*integration_length))]]';

    ints = {ints_forward, ints_backward};
    % Loop through and perform integrations
    for type = 1:2
        if type==2
            disp('Backward')
            ints_backward = ints_backward(1:int_count,:);
            ints = {ints_forward, ints_backward};
        else
            disp('Forward')
        end


        int_count = 0;
        int = cell2mat(ints(type));
        
        if singleInt  % deal with single integrations
            numInts = 1;
        else
            numInts = size(int,1);
        end

        for i = 1:numInts
            time_span = int(i,:);

            % compute cauchy green
            [cg_eig_vect, cg_eig_val, cgStrain] = eig_cgStrain_v2_NJL(lDeriv, domain, res, time_span,'eigenvalueFromMainGrid', true);

            % get maximal eigValue, reshape, and run FTLE.
            cg_eig_val2 = reshape(cg_eig_val(:,2), fliplr(resolution));
            ftle_ = ftle(cg_eig_val2, abs(diff(time_span)));

            % compute det of Cauchy-green (Note this is just the product of the
            % Eigenvalues since it is a symmetric matrix)
            eigen_prod = eig_prod(cg_eig_val);

            % reshape to fit on grid
            eigen_prod = reshape(eigen_prod, fliplr(resolution));
            if type == 1
                forward.FTLE = cat(3, forward.FTLE, ftle_);
                forward.cgStrain = cat(3, forward.cgStrain, cgStrain);
                forward.eigen_prod = cat(3, forward.eigen_prod, eigen_prod);
            else
                rever.FTLE = cat(3, rever.FTLE, ftle_);
                rever.cgStrain = cat(3, rever.cgStrain, cgStrain);
                rever.eigen_prod = cat(3, rever.eigen_prod, eigen_prod);
            end

            int_count = int_count + 1;
            % update user
            update = strcat('Completed int ', num2str(int_count));
            % fprintf(update);
            fprintf('.')
        end
    end
end

% FTLE Calculate Finite-time Lyapunov exponent
function ftle_ = ftle(max_eigenvalue,timespan)
    ftle_ = .5*log(max_eigenvalue)/timespan;
end

% COMPUTE determinant or product of eigenvalues of the cgStrain tensor
function prod_ = eig_prod(eigenvalues,timespan)
    prod_ = prod(eigenvalues,2);
end

