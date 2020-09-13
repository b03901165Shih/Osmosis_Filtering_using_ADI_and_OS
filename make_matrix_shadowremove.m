function [A1, A2, A12]=make_matrix_shadowremove(Im, mask, ch, WW, mask_loc)
%[M0,Md,Am,x]=make_matrix(M,N);
% creates matrices to solve advection-diffusion problem
%INPUTS: M specifies the spatial domain x1=(-M1,M1), x2=(-M2,M2)
%N1,N2 are the grid size (size of matrices (N1*N2)x(N*N2) and vectors
%(N1*N2)x1)
%OUTPUT: D1,D2,D11,D22,x1,x2 are used to solve eqn

[d1, d2] = make_drift_vector_shadowremove(Im, mask);
d1 = d1(:,:,ch);  d2 = d2(:,:,ch);
H = size(Im,1); W = size(Im, 2);

%----Insert commands for creating D1, D2, and Am here
% A1:
a = []; b = []; c = []; 
a_mat = []; c_mat = []; 
f = [];
b_coeff = zeros(H,1);
for j=1:W
    a_coeff = (WW(2:H,j,1)+WW(1:H-1,j,1))/2 +d1(2:H,j);
    b_coeff(2:H-1) = -(WW(3:H,j,1)+2*WW(2:H-1,j,1)+WW(1:H-2,j,1))/2 -d1(3:H, j)+d1(2:H-1,j);
    c_coeff =  (WW(2:H,j,1)+WW(1:H-1,j,1))/2-d1(2:H, j); 
    % i==0 (boundary)
    b_coeff(1) = -(WW(1,j,1)+WW(2,j,1))/2-d1(2,j); 
    c_coeff(1) = (WW(1,j,1)+WW(2,j,1))/2-d1(2,j); 
    % i==H-1 (boundary)
    b_coeff(H) = -(WW(H,j,1)+WW(H-1,j,1))/2+d1(H,j); 
    a_coeff(H-1) = (WW(H,j,1)+WW(H-1,j,1))/2+d1(H,j); 
    
    D = mask_loc(:,j);
    D_neg = 1-D;
    a_neg = ([0; D(1:end-1)]&D_neg).*[0; Im(1:end-1,j, ch)];
    c_neg = ([D(2:end); 0]&D_neg).*[Im(2:end,j, ch);0]; 
    
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
b_coeff = zeros(1, W);
for i=1:H
    a_coeff = (WW(i,2:W,4)+WW(i,1:W-1,4))/2+d2(i, 2:W);
    b_coeff(2:W-1) = -(WW(i,3:W,4)+2*WW(i,2:W-1,4)+WW(i,1:W-2,4))/2 -d2(i, 3:W)+d2(i, 2:W-1);
    c_coeff = (WW(i,2:W,4)+WW(i,1:W-1,4))/2-d2(i, 2:W); 
    % i==0 (boundary)
    b_coeff(1) = -(WW(i,1,4)+WW(i,2,4))/2-d2(i,2); % actually not needed
    c_coeff(1) = (WW(i,1,4)+WW(i,2,4))/2-d2(i,2); 
    % i==H-1 (boundary)
    b_coeff(W) = -(WW(i,W,4)+WW(i,W-1,4))/2+d2(i, W); % actually not needed
    a_coeff(W-1) = (WW(i,W,4)+WW(i,W-1,4))/2+d2(i, W); 
    
    D = mask_loc(i,:);
    D_neg = 1-D;
    a_neg = ([0, D(1:end-1)]&D_neg).*[0, Im(i,1:end-1,ch)];
    c_neg = ([D(2:end), 0]&D_neg).*[Im(i,2:end,ch),0]; 
    
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
%%{
%% First matrix ([u_i0, u_i1, ..., u_iw-1, u_i0,...])
D12_mat = spalloc(H*W, H*W, 4*H*W);
values = zeros(H,W,4);
for i=1:H
    %fprintf('Row index: %d\r',i);
	for j=1:W
        if(mask_loc(i,j)>=0.5) 
            continue;
        end
		dind = (i-1)*W+j;	%diagonal index
        i_decr = (i==1)*i+(i~=1)*(i-1); 
        i_incr = (i==H)*i+(i~=H)*(i+1); sign_i_decr = 1-2*(i==1); sign_i_incr = 1-2*(i==H);
        j_decr = (j==1)*j+(j~=1)*(j-1); 
        j_incr = (j==W)*j+(j~=W)*(j+1); sign_j_decr = 1-2*(j==1); sign_j_incr = 1-2*(j==W);
        coeff = 1/4;%1/(4-(j==W || j==1)*2-(i==H || i==1)*2+((i==1 || i==H) && (j==1 || j==W))*1);
        value1 = (sign_i_decr*WW(i_decr,j,2)+sign_j_decr*WW(i,j_decr,3))*coeff; % i-1, j-1
        value2 = (sign_i_decr*WW(i_decr,j,2)+sign_j_incr*WW(i,j_incr,3))*coeff;%i-1, j+1
        value3 = (sign_i_incr*WW(i_incr,j,2)+sign_j_decr*WW(i,j_decr,3))*coeff;%i+1, j-1
        value4 = (sign_i_incr*WW(i_incr,j,2)+sign_j_incr*WW(i,j_incr,3))*coeff;%i+1, j+1
        values(i,j,:) = [value1, value2, value3, value4];
        if(value1==0 && value2 == 0 && value3==0 && value4 == 0)
            continue;
        end
		if((i==1 || i==H) && (j==1 || j==W))
			i_shift = (i==1)*0+(i==H)*-W;
			j_shift = (j==1)*0+(j==W)*-1;
            sign=-1+2*((i==1&&j==1)||(i==H&&j==W));
			D12_mat(dind, dind+i_shift+j_shift) 	= 1*value1; 
			D12_mat(dind, dind+i_shift+j_shift+1) = -1*value2;
			D12_mat(dind, dind+i_shift+j_shift+W) = -1*value3; 
			D12_mat(dind, dind+i_shift+j_shift+1+W) = 1*value4;
		elseif((i==1 || i==H))
			i_shift = (i==1)*0+(i==H)*-W;
            sign13=-1+2*(i==H); sign24=-1+2*(i==1);
			D12_mat(dind, dind+i_shift-1) = 1*value1; 
			D12_mat(dind, dind+i_shift+1) = -1*value2;
			D12_mat(dind, dind+i_shift+W-1)= -1*value3; 
			D12_mat(dind, dind+i_shift+W+1)= 1*value4;
		elseif((j==1 || j==W))
			j_shift = (j==1)*0+(j==W)*-1;
            sign12=-1+2*(j==W); sign34=-1+2*(j==1);
			D12_mat(dind, dind+j_shift-W) = 1*value1; 
			D12_mat(dind, dind+j_shift-W+1) = -1*value2;
			D12_mat(dind, dind+j_shift+W)= -1*value3; 
			D12_mat(dind, dind+j_shift+W+1)= 1*value4;
        else
			D12_mat(dind, dind-W-1) = value1;   % i-1, j-1
			D12_mat(dind, dind-W+1) = -value2; %i-1, j+1
			D12_mat(dind, dind+W-1)= -value3;  %i+1, j-1
			D12_mat(dind, dind+W+1)= value4;  %i+1, j+1
		end
	end
end
A12.matrix1 = D12_mat;
A12.values1 = values;
%% Second matrix ([u_0j, u_1j, ..., u_h-1j, u_0j,...])
D12_mat = spalloc(H*W, H*W, 4*H*W);
values = zeros(H,W,4);
for j=1:W
    %fprintf('Row index: %d\r',i);
	for i=1:H
        if(mask_loc(i,j)>=0.5) 
            continue;
        end
		dind = (j-1)*H+i;	%diagonal index
        i_decr = (i==1)*i+(i~=1)*(i-1); 
        i_incr = (i==H)*i+(i~=H)*(i+1); sign_i_decr = 1-2*(i==1); sign_i_incr = 1-2*(i==H);
        j_decr = (j==1)*j+(j~=1)*(j-1); 
        j_incr = (j==W)*j+(j~=W)*(j+1); sign_j_decr = 1-2*(j==1); sign_j_incr = 1-2*(j==W);
        coeff = 1/4;%1/(4-(j==W || j==1)*2-(i==H || i==1)*2+((i==1 || i==H) && (j==1 || j==W))*1);
        value1 = (sign_i_decr*WW(i_decr,j,2)+sign_j_decr*WW(i,j_decr,3))*coeff; % i-1, j-1
        value2 = (sign_i_incr*WW(i_incr,j,2)+sign_j_decr*WW(i,j_decr,3))*coeff;%i+1, j-1
        value3 = (sign_i_decr*WW(i_decr,j,2)+sign_j_incr*WW(i,j_incr,3))*coeff;%i-1, j+1
        value4 = (sign_i_incr*WW(i_incr,j,2)+sign_j_incr*WW(i,j_incr,3))*coeff;%i+1, j+1
        values(i,j,:) = [value1, value2, value3, value4];
        if(value1==0 && value2 == 0 && value3==0 && value4 == 0)
            continue;
        end
		if((i==1 || i==H) && (j==1 || j==W))
			i_shift = (i==1)*0+(i==H)*-1;
			j_shift = (j==1)*0+(j==W)*-H;
            sign=-1+2*((i==1&&j==1)||(i==H&&j==W));
			D12_mat(dind, dind+i_shift+j_shift) 	= 1*value1; 
			D12_mat(dind, dind+i_shift+j_shift+1) 	= -1*value2;
			D12_mat(dind, dind+i_shift+j_shift+H) = -1*value3; 
			D12_mat(dind, dind+i_shift+j_shift+1+H) = 1*value4;
		elseif((i==1 || i==H))
			i_shift = (i==1)*0+(i==H)*-1;
            sign12=-1+2*(i==H); sign34=-1+2*(i==1);
			D12_mat(dind, dind+i_shift-H) = 1*value1; 
			D12_mat(dind, dind+i_shift-H+1) = -1*value2;
			D12_mat(dind, dind+i_shift+H)= -1*value3; 
			D12_mat(dind, dind+i_shift+H+1)= 1*value4;
		elseif((j==1 || j==W))
			j_shift = (j==1)*0+(j==W)*-H;
            sign13=-1+2*(j==W); sign24=-1+2*(j==1);
			D12_mat(dind, dind+j_shift-1) =1*value1; 
			D12_mat(dind, dind+j_shift+1) = -1*value2;
			D12_mat(dind, dind+j_shift+H-1)= -1*value3; 
			D12_mat(dind, dind+j_shift+H+1)= 1*value4;
		else
			D12_mat(dind, dind-H-1) = value1; 
			D12_mat(dind, dind-H+1) = -value2;
			D12_mat(dind, dind+H-1)= -value3; 
			D12_mat(dind, dind+H+1)= value4;
		end
	end
end
A12.matrix2 = D12_mat;
A12.values2 = values;