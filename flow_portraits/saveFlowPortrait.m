% Nathaniel Linden
% Bing Brunton Lab - University of Washington
% June 29th, 2020

% Function to save FLOW portrait to an image file
function saveFlowPortrait(for_FLOW, back_FLOW, mean_dFoF, filename, white_background, mask, for_color, back_color, background_cmap, figsize)


% Syntax: saveFlowPortrait(for_FLOW, back_FLOW, mean_dFoF, filename, , white_background, mask, for_color, back_color, background_cmap, figsize)

% Inputs:
% for_FLOW, back_FLOW -   FLOW fields from create_flow_portrait(), dimensions [num_rows, num_cols]
% mean_dFoF -             mean data in time
% filename -              string with filename to save flow portrait, includes file extension (eg. .png for a PNG) - defualts to flow_image.png
% for_color, back_color - RGB vectors [R, G, B] for the forward and backward FLOW colors (Note: RGB vals must be between 0 & 1)
% background_cmap -       colormap for the background image - default to GRAY
% figsize -               position vector for size of figure, format: [0, 0, WIDTH, HEIGHT]
% white_background -      logical, true - use nans to create a white bacground outside field of view
% mask -                  logical matrix, dimensions [num_rows, num_cols], masks field of view
    
    % % defualt parameters
    % if nargin < 5 || nargin < 6; white_background = false; mask = []; end

    % Create Figure and Plot
    figure('Renderer', 'painters', 'Position', figsize,'color','w');

    % adjust range of backgorund image
    scale = [min(min(mean_dFoF)), max(max(mean_dFoF))];
    upper_val = 0.65; 
    if scale(2) < upper_val
        mean_dFoF = imadjust(mean_dFoF, scale, [0 upper_val]);
    else
        mean_dFoF = imadjust(mean_dFoF);
    end
    
    % white_background
    if white_background
        [newImg, imAlpha] = whiteBackground(mean_dFoF, mask);
        mean_dFoF = newImg;
    else
        imAlpha = ones(size(mean_dFoF));
    end

    % create overlay image and plot
    img1 = imoverlay(mean_dFoF, for_FLOW, for_color);
    img2 = imoverlay(img1, back_FLOW, back_color);
    imagesc(img2, 'AlphaData',imAlpha);
   
    % Format and Save Image
    ax = gca;
    ax.CLim = scale;
    axis off; axis tight equal 
    colormap(background_cmap);
    ax.Position = ax.OuterPosition;
    saveas(gcf, filename);

    close all
end
