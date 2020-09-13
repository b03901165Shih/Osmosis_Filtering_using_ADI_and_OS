clear;

% Hyper-parameters
local = false;
T_end = 5000;
dt = 10;
scheme = 'ADI'; %otherwise: 'OS'
input_folder = 'imgs/fusion/';
output_folder = 'results/fusion/';

if ~exist(output_folder, 'dir')
   mkdir(output_folder)
end
%%
Im1 = im2double(imread([input_folder,'p1.png']));
Im2 = im2double(imread([input_folder,'p2.png']));
mask = im2double(imread([input_folder,'alpha.png']));

if(local)
    % Local (setting Dirichlet boundary condition) (fixed value for mask_loc = 1)
    mask_loc = double(mask(:,:,1)>=0.99);
else
    % Global
    mask_loc = zeros(size(mask, 1), size(mask, 2));
end

% Initialize image
I_init = (mask).*Im1+(1-mask).*Im2;
H = size(I_init,1); W = size(I_init, 2);

t=(0:dt:T_end)';
p = I_init;
for ch = 1:3
    fprintf('Processing channel %d...\n', ch);
    [A1, A2]=make_matrix_fusion(p, Im1, Im2, mask, ch, mask_loc);
    [t,p]=solve_advDif(p, A1, A2, [], t, ch, 'ADI'); %(p, A1, A2, A12, t, ch, scheme)
end

figure, imshow([I_init, p(:,:,:,end)]), title('Left: original/ Right: result');
imwrite(I_init, [output_folder, 'fusion_input.png'])
imwrite(p, [output_folder, 'fusion_output.png'])
