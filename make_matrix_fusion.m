function [A1, A2]=make_matrix_fusion(I_init, Im1, Im2, mask, ch, mask_loc)
%[M0,Md,Am,x]=make_matrix(M,N);
% creates matrices to solve advection-diffusion problem
%INPUTS: M specifies the spatial domain x1=(-M1,M1), x2=(-M2,M2)
%N1,N2 are the grid size (size of matrices (N1*N2)x(N*N2) and vectors
%(N1*N2)x1)
%OUTPUT: D1,D2,D11,D22,x1,x2 are used to solve eqn

[d1, d2] = make_drift_vector_fusion(Im1, Im2, mask);
d1 = d1(:,:,ch);  d2 = d2(:,:,ch);
H = size(Im1,1); W = size(Im1, 2);

%----Insert commands for creating D1, D2, and Am here
% A1:
a = []; b = []; c = []; 
a_mat = []; c_mat = []; 
f = [];
for j=1:W
    a_coeff = 1+d1(2:H,j);
    b_coeff = -2-d1(2:H+1, j)+d1(1:H,j);
    c_coeff = 1-d1(2:H, j); 
    % i==0 (boundary)
    b_coeff(1) = -1-d1(2,j); % actually not needed
    c_coeff(1) = 1-d1(2,j); 
    % i==H-1 (boundary)
    b_coeff(H) = -1+d1(H,j); % actually not needed
    a_coeff(H-1) = 1+d1(H,j); 
    
    D = mask_loc(:,j);
    D_neg = 1-D;
    a_neg = ([0; D(1:end-1)]&D_neg).*[0; I_init(1:end-1,j, ch)];
    c_neg = ([D(2:end); 0]&D_neg).*[I_init(2:end,j, ch);0]; 
    
    ac_mask = D_neg(2:end)&D_neg(1:end-1);
    %a_coeff = a_coeff.*ac_mask;
    %c_coeff = c_coeff.*ac_mask;
    b_coeff = b_coeff.*D_neg;
        
    a = [a, 0, (a_coeff.*ac_mask)'];
    c = [c, (c_coeff.*ac_mask)', 0];
    b = [b, b_coeff'];
    
    a_mat = [a_mat, 0, (a_coeff.*D_neg(2:end))'];
    c_mat = [c_mat, (c_coeff.*D_neg(1:end-1))', 0];
    
    f = [f, (a_neg+c_neg)'];
end
A1.a = a; A1.b = b; A1.c = c; A1.f = f;
A1.matrix = make_tridiag(a_mat,b,c_mat);

%%%%%%%%%%%%%%%%%%%%%%%%%%%

% A2:
a = []; b = []; c = []; 
a_mat = []; c_mat = []; 
f = [];
for i=1:H
    a_coeff = 1+d2(i, 2:W);
    b_coeff = -2-d2(i, 2:W+1)+d2(i, 1:W);
    c_coeff = 1-d2(i, 2:W); 
    % i==0 (boundary)
    b_coeff(1) = -1-d2(i,2); % actually not needed
    c_coeff(1) = 1-d2(i,2); 
    % i==H-1 (boundary)
    b_coeff(W) = -1+d2(i, W); % actually not needed
    a_coeff(W-1) = 1+d2(i, W); 
    
    D = mask_loc(i,:);
    D_neg = 1-D;
    a_neg = ([0, D(1:end-1)]&D_neg).*[0, I_init(i,1:end-1,ch)];
    c_neg = ([D(2:end), 0]&D_neg).*[I_init(i,2:end,ch),0]; 
    
    ac_mask = D_neg(2:end)&D_neg(1:end-1);
    %a_coeff = a_coeff.*ac_mask;
    %c_coeff = c_coeff.*ac_mask;
    b_coeff = b_coeff.*D_neg;
    
    a = [a, 0, (a_coeff.*ac_mask)];
    c = [c, (c_coeff.*ac_mask), 0];
    b = [b, b_coeff];    
    
    a_mat = [a_mat, 0, (a_coeff.*D_neg(2:end))];
    c_mat = [c_mat, (c_coeff.*D_neg(1:end-1)), 0];
    
    f = [f, (a_neg+c_neg)];
end
A2.a = a; A2.b = b; A2.c = c; A2.f = f;
A2.matrix = make_tridiag(a_mat,b,c_mat);

%%%%%%%%%%%%%%%%%%%%%%%%%%%