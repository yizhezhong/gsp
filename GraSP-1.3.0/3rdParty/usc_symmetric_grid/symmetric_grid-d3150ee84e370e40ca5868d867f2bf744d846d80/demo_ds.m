% This is a demo for block diagonalization implementation in 
% proposed_blk_diagonalization.m
% 
% by: KS Lu
% last update: 20170106

% this is a Laplacian matrix corresponding to a bidiagonal symmetric 
% grid
load L_bs
openfig('gridexample');

[H,R]=proposed_blk_diagonalization(L,'bs');

figure; imagesc(L); colorbar; 
title('Graph Laplacian matrix L, which has a factorization L=HRH^T');

figure; imagesc(R); colorbar; hold on;
NN=size(H,1); N=sqrt(NN);
x=[(NN-2*N)/4, (NN-N)/2, (3*NN-2*N)/4];
for j=1:3
    plot([0,NN+1], [x(j)+0.5, x(j)+0.5], 'w--', 'linewidth', 2);
    plot([x(j)+0.5, x(j)+0.5], [0,NN+1], 'w--', 'linewidth', 2);
end
%set(gca,'visible','off');
title('The block-diagonalized matrix R that is similar to the Laplacian matrix');

figure; imagesc(H); colorbar;
%set(gca,'visible','off');
title('The butterfly-related matrix H for bidiagonal symmetry');