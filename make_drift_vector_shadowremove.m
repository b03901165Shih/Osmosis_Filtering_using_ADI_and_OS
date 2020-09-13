% Im: shadowed image
% mask: binary mask for shadow boundary (0 for boundary, 1 for others)
% drift_vec1/ vec2(output): drift vector for shadow removal
% drift_Im(inter): image for calculating drift vector for shadow removal
function [d1, d2] = make_drift_vector_shadowremove(Im, mask)
    mask = double(mask>=0.5);
    
    %masked_image = ones(size(Im,1), size(Im,2), 3);
    %masked_image(:,:,1) = mean2(Im(:,:,1));
    %masked_image(:,:,2) = mean2(Im(:,:,2));
    %masked_image(:,:,3) = mean2(Im(:,:,3));
    %drift_Im = Im.*mask+masked_image.*(1-mask);
    
    pad_i_Im = padarray(Im,[1,0],'symmetric'); % top/bottom
    pad_j_Im = padarray(Im,[0,1],'symmetric'); % left/right
    
    % d1: (h+1)*w
    d1 = (pad_i_Im(2:end,:,:)-pad_i_Im(1:end-1,:,:))./(pad_i_Im(2:end,:,:)+pad_i_Im(1:end-1,:,:)+1e-3);
    % d2: h*(w+1)
    d2 = (pad_j_Im(:,2:end,:)-pad_j_Im(:,1:end-1,:))./(pad_j_Im(:,2:end,:)+pad_j_Im(:,1:end-1,:)+1e-3);
    
    % reverse mask and pad
    pad_i_mask = padarray(1-mask,[1,0],0,'both'); % top/bottom
    pad_j_mask = padarray(1-mask,[0,1],0,'both'); % left/right
    
    mask_i = min(pad_i_mask(2:end,:,:)+pad_i_mask(1:end-1,:,:), 1);
    mask_j = min(pad_j_mask(:,2:end,:)+pad_j_mask(:,1:end-1,:), 1);
    
    d1(mask_i==1) = 0;
    d2(mask_j==1) = 0;
    
end

