function [t,p]=solve_advDif(p, A1, A2, A12, t, ch, scheme)

H = size(p,1); W = size(p, 2);
dt = t(2)-t(1);

tic() 
p_ch = reshape(p(:,:,ch)', H*W, 1);
%ADI scheme    
assert(strcmp(scheme, 'ADI')||strcmp(scheme, 'OS'));
if(strcmp(scheme, 'ADI')) % Second-order accurate scheme (not suitable for problems with A12(mixed derivatives))
    rhs_first     = sparse(speye(H*W)+.5*dt*(A2.matrix));
    rhs_second = sparse(speye(H*W)+.5*dt*(A1.matrix));

    a_first =-.5*dt*(A1.a);
    b_first =ones(H*W,1)'-.5*dt*(A1.b);
    c_first =-.5*dt*(A1.c);

    a_second =-.5*dt*(A2.a);
    b_second =ones(H*W,1)'-.5*dt*(A2.b);
    c_second =-.5*dt*(A2.c);
    c=0.5;
else %OS scheme: first-order accurate scheme
    rhs_first     = sparse(speye(H*W)+.5*dt*(A12.matrix1));
    rhs_second = sparse(speye(H*W)+.5*dt*(A12.matrix2)); %+.5*dt*(A1.matrix)

    a_first =-dt*(A1.a);
    b_first =ones(H*W,1)'-dt*(A1.b);
    c_first =-dt*(A1.c);

    a_second =-dt*(A2.a);
    b_second =ones(H*W,1)'-dt*(A2.b);
    c_second =-dt*(A2.c);

    c=1;
end

for j=2:length(t)
        if rem(j,100)==0
            fprintf('Number of Iteration %d out of %d finished...\n',j,length(t));
        end
        %----solve for p(:,j)=....(:,j-1)(:,j)
        reshaped_r = reshape(reshape(rhs_first*p_ch, W, H)',H*W,1)-(-c*dt*A1.f');
        p_half = tridiag(a_first,b_first,c_first,reshaped_r); 

        reshaped_r = reshape(reshape(rhs_second*p_half, H, W)',H*W,1)-(-c*dt*A2.f');
        p_ch = tridiag(a_second,b_second,c_second,reshaped_r);  
        %p(:,:,ch, j) = reshape(p_ch,W,H)';
end
p(:,:,ch) = reshape(p_ch,W,H)';    
    
toc()
