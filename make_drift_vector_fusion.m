% Im: shadowed image
% mask: binary mask for shadow boundary (0 for boundary, 1 for others)
% drift_vec1/ vec2(output): drift vector for shadow removal
% drift_Im(inter): image for calculating drift vector for shadow removal
function [d1, d2] = make_drift_vector_fusion(Im1, Im2, mask)
    %masked_image = ones(size(Im,1), size(Im,2), 3);
    %masked_image(:,:,1) = mean2(Im(:,:,1));
    %masked_image(:,:,2) = mean2(Im(:,:,2));
    %masked_image(:,:,3) = mean2(Im(:,:,3));
    %drift_Im = Im.*mask+masked_image.*(1-mask);
    
    pad_i_Im1 = padarray(Im1,[1,0],'symmetric'); % top/bottom
    pad_j_Im1 = padarray(Im1,[0,1],'symmetric'); % left/right
    pad_i_Im2 = padarray(Im2,[1,0],'symmetric'); % top/bottom
    pad_j_Im2 = padarray(Im2,[0,1],'symmetric'); % left/right
    
    % d1: (h+1)*w
    d1_Im1 = (pad_i_Im1(2:end,:,:)-pad_i_Im1(1:end-1,:,:))./(pad_i_Im1(2:end,:,:)+pad_i_Im1(1:end-1,:,:)+1e-3);
    % d2: h*(w+1)
    d2_Im1 = (pad_j_Im1(:,2:end,:)-pad_j_Im1(:,1:end-1,:))./(pad_j_Im1(:,2:end,:)+pad_j_Im1(:,1:end-1,:)+1e-3);
    
    % d1: (h+1)*w
    d1_Im2 = (pad_i_Im2(2:end,:,:)-pad_i_Im2(1:end-1,:,:))./(pad_i_Im2(2:end,:,:)+pad_i_Im2(1:end-1,:,:)+1e-3);
    % d2: h*(w+1)
    d2_Im2 = (pad_j_Im2(:,2:end,:)-pad_j_Im2(:,1:end-1,:))./(pad_j_Im2(:,2:end,:)+pad_j_Im2(:,1:end-1,:)+1e-3);
    
    % reverse mask and pad
    pad_mask = padarray(mask,[1,1],1,'both'); % top/bottom
    pad_mask = pad_mask(2:end, 2:end,:);
    
    d1 = d1_Im1.*pad_mask(1:end,1:end-1,:)+d1_Im2.*(1-pad_mask(1:end,1:end-1,:));
    d2 = d2_Im1.*pad_mask(1:end-1,1:end,:)+d2_Im2.*(1-pad_mask(1:end-1,1:end,:));
    
end

