% Function to make a masked region of an image NaN values
function [newStack, imAlpha] = whiteBackground(imgStack, mask)
    % create new masked stack
    newStack = imgStack;
    mask = mask .* ones(size(newStack)); 
    newStack(~(mask))=NaN;  % set zeroed masked regions to NaN

    % create imAlpha to plot NaN as no color
    imAlpha = ones(size(newStack));
    imAlpha(isnan(newStack))=0;
    imAlpha = imAlpha(:,:,1);
end