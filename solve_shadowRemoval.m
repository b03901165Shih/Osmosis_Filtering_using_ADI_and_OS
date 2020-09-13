clear;
addpath(genpath('anisotropic-osmosis-filter'));

testcases = {'44','ball', '017'};
input_folder = 'imgs/shadow/';
output_folder = 'results/shadow/';
anisotrpic_diffusion = true;
local = true;
T_end = {1000, 5000, 150000};
dt = 10;
if anisotrpic_diffusion scheme = 'OS'; else scheme='ADI'; end

if ~exist(output_folder, 'dir')
   mkdir(output_folder)
end

%%
for l = 1:length(testcases)

    Im = im2double(imread([input_folder, testcases{l}, '.png']));
    mask = im2double(imread([input_folder, testcases{l}, '_mask.png']));
    mask = double(mask>=0.5);
    
    if(local)
        %mask_loc(:,:)=0;
        mask_loc = im2double(imread([input_folder, testcases{l}, '_mask_local.png']));
        mask_loc = double(mask_loc(:,:,1)>=0.5);
    else
        mask_loc = zeros(size(Im, 1), size(Im, 2));
    end

    if(anisotrpic_diffusion)
        % b is the minor axis of ellipse b=1 is the circle b=0 is the segment (degenerate ellipse)
        b{1}  = 0.05;
        b{2}  = 1;
        sigma = 0.25;
        rho   = 2;
        [~,TVF] = computefield(Im,mask(:,:,1),sigma,rho);
        [~,WW] = compute_W(TVF.e1,b,mask);
    else %Isotropic
        WW = zeros(size(Im, 1), size(Im, 2), 4);
        WW(:,:,1) = 1;
        WW(:,:,4) = 1;
    end

    %Global


    % Initialize image
    I_init = Im;
    H = size(I_init,1); W = size(I_init, 2);

    t=(0:dt:T_end{l})';
    p = I_init;
    for ch = 1:size(Im, 3)
        fprintf('Processing channel %d...\n', ch);
        [A1, A2, A12]=make_matrix_shadowremove(p, mask, ch, WW, mask_loc);
        [t,p]=solve_advDif(p, A1, A2, A12, t, ch, scheme); %(p, A1, A2, A12, t, ch, scheme)
    end

    figure, imshow([Im, p]), title('Left: original/ Right: shadow removal result');
    imwrite(p, [output_folder, testcases{l}, '_output_aniso.png'])
end
