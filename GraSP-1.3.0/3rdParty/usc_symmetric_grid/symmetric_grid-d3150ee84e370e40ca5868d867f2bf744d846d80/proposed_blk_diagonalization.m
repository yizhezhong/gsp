function [H, R, GH, DH] = proposed_blk_diagonalization(L, symtype)
% [H, R, GH, DH]=proposed_blk_diagonalization(L, symtype)
%    return the matrices H and R for the factorization L=HRH^T given 
%    Laplacian matrix L of a symmetric grid with known symmetry type. 
%    H can be factorized as H=GH*DH, as in [1]
%    
% L: Laplacian matrix of a grid, whose size must be N^2xN^2 with some 
%    integer N. It must correpond to a grid (with column-first vertex 
%    ordering) that has the symmetry type specified by the next input, 
%    i.e., symtype. 
% 
% symtype: grid symmetry type, which can be 
%    'ud' for UD-symmetry 
%    'lr' for LR-symmetry
%    'udlr' for UDLR-symmetry 
%    'centro' for centrosymmetry
%    'ds' for diagonal symmetry
%    'as' for anti-diagonal symmetry 
%    'bs' for bidiagonal symmetry 
%    'penta' for pentasymmetry
% 
% H: orthogonal matrix H, as in [1]. It is a scaled version of the 
%    butterfly. (It is the same as the first output of the function hmat.)
% 
% R: butterfly-related block diagonal orthogonal matrix, as in [1]. 
% 
% GH: the butterfly matrix, whose elements are all 0,1,-1. (It is the same 
%    as the second output of the function hmat.)
% 
% DH: diagonal matrix describing the constant factors 
% 
% [1] K.-S. Lu and A. Ortega, "Fast Implementations for Symmetric Non-
%    Separable Transforms Based on Grids," ICASSP 2017.
% 
% by: KS Lu
% last update: 20170105

NN=size(L,1);
N=sqrt(NN);
L_offdiag = L(:); L_offdiag(1:NN+1:end)=0;
if size(L,1)~=size(L,2) || round(N)~=N || any(L_offdiag>1e-8) ...
        || any(diag(L)<-1e-8) || any(sum(L)<-1e-8) || ~issymmetric(L)
    error('L is not a valid Laplacian matrix for a grid');
end

% verify that the Laplacian matrix satisfies the given symmetry
if ~check_grid_sym(L,symtype)
    error('L does not satisfy the given symmetry type.');
end

% perform block diagonalization using the function hmat
[H, GH] = hmat(N,symtype);
fac = sqrt(sum(GH.^2));
DH = diag(fac);
R = H'*L*H;