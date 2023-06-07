function [H, H1]=hmat(gridsize, symtype)
% [H, H1]=hmat(gridsize, symtype)
%    return the butterfly-related matrix H corresponding to a given grid 
%    size and a given grid symmetry type
% 
% gridsize: grid size. For ud/lr/udlr/centro symmetry, 
%    gridsize=[height, width]. For other symmetry types, it is required 
%    that height=width. If this input is a scalar N, it is regarded as a 
%    square grid with N=height=width. In the current version, both width 
%    and height have to be even.
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
%    butterfly
%
% H1: butterfly matrix H, whose elements are all 0,1,-1
%
% [1] K.-S. Lu and A. Ortega, "Fast Implementations for Symmetric Non-
%    Separable Transforms Based on Grids," ICASSP 2017.
% 
% by: KS Lu
% last update: 20170105


% handle grid height and width
if numel(gridsize)==1
    h=gridsize;
    w=gridsize;
else
    h=gridsize(1);
    w=gridsize(2);
end

if any(gridsize<=0) || any(round(gridsize)~=gridsize) || any(mod(gridsize,2)==1)
    error('Invalid grid size');
end

switch symtype
    case {'ds', 'as', 'bs', 'penta'}
        if h~=w
            error('For ds/as/bs/penta symmetry types, we require height=width');
        end
end

% some variables
NN=h*w;
I_NN = eye(NN);
I_NN_half = eye(NN/2);
matI = reshape(1:NN,h,w);
matI_half = reshape(1:NN/2,w,h/2);  % for recovering the order after H_ud
lB=(NN-2*h)/4;  % length of key sub-blocks

% prepare node reordering permutation matrices
switch symtype
    case {'ud', 'udlr'}
        
        % here we use a row-first order to exploit the UD-symmetry
        idx_ud = matI'; idx_ud=idx_ud(:);
        P_ud = I_NN(:,idx_ud);
        
        % This is to recover the node ordering to column-first order.
        % Note that this is not a inverse matrix of P_ud, because the 
        % column-first scanning must be applied to each component 
        % independently.
        idx_ud_recov = matI_half'; idx_ud_recov=idx_ud_recov(:);
        P_ud_recov = kron(eye(2), I_NN_half(:,idx_ud_recov));
        
    case {'bs', 'penta'}
        
        % vertex scanning for bidiagonal symmetry and pentasymmetry cases
        idx_diags=[diag(matI), diag(fliplr(matI))];
        idx_diags(h/2+1:end,:)=idx_diags(h/2+1:end,[2,1]);
        matI_0diag=matI; 
        matI_0diag(1:h+1:end)=0; matI_0diag(h:h-1:end)=0;
        % matI_0diag contains the labels 1 to N^2 in a column-first 
        % manner, but with both diagonals removed
        matI_u=triu(fliplr(triu(fliplr(matI_0diag))));
        idx_u=rot90(matI_u,2);
        idx_u(:,1:2:end)=flipud(idx_u(:,1:2:end));
        idx_u=idx_u(:)'; idx_u(idx_u==0)=[];
        matI_d=tril(fliplr(tril(fliplr(matI_0diag))));
        idx_d=matI_d;
        idx_d(:,1:2:end)=flipud(idx_d(:,1:2:end));
        idx_d=idx_d(:)'; idx_d(idx_d==0)=[];
        matI_l=tril(flipud(tril(flipud(matI_0diag))));
        idx_l=flipud(matI_l');
        idx_l(:,1:2:end)=flipud(idx_l(:,1:2:end));
        idx_l=idx_l(:)'; idx_l(idx_l==0)=[];
        matI_r=triu(flipud(triu(flipud(matI_0diag))));
        idx_r=fliplr(matI_r');
        idx_r(:,1:2:end)=flipud(idx_r(:,1:2:end));
        idx_r=idx_r(:)'; idx_r(idx_r==0)=[];
        idx_bs=[idx_diags(:)', idx_l, idx_d, idx_r, idx_u];
        P_bs = I_NN(:,idx_bs);
        
end

switch symtype
    case 'centro'
        H1=kmat1(NN)';
    case 'ud'
        
        % This is the version in [1]. It works, but the node order will 
        % become row-first. 
        % H1=P_ud*kron(kmat1(h)',eye(w)); 
        
        % This version keeps the column-first order, which simplifies 
        % further decomposition in the next levels
        H1=P_ud*kron(kmat1(h)',eye(w))*P_ud_recov;
        
    case 'lr'
        H1=kron(kmat1(w)',eye(h));
    case 'udlr'
        idx_delta=[matI(1:h/2,:),matI(h/2+1:end,:)]; idx_delta=idx_delta(:);
        P_delta=I_NN(:,idx_delta);
        H1=kron(kmat1(w)',eye(h))*kron(eye(w),kmat1(h)')*P_delta;
        
        % H1=kron(kmat1(w)',eye(h))*P_ud*kron(kmat1(h)',eye(w))*P_delta; %
        % This should work, but currently there is a bug that 
        % should be corrected later
        
    case 'ds'
        idx_ds = zeros(1,NN);
        idx_ds(1:h)=diag(matI,0)';
        curI=h;
        for i=[-h+1:-1, 1:h-1]
            idx_in_zigzag = curI+(1:(h-abs(i)));
            if (i<0 && mod(i,2)==0) || (i>0 && mod(i,2)==1)  % reverse order
                idx_ds(fliplr(idx_in_zigzag)) = diag(matI,i)';
            else
                idx_ds(idx_in_zigzag) = diag(matI,i)';
            end
            curI=curI+h-abs(i);
        end
        P_ds = I_NN(:,idx_ds);
        %P_alpha=I_NN(:,[1:h, h*(h+1)/2+1:NN, h+1:h+h*(h-1)/2]);  % large block first
        P_alpha=I_NN(:,[h+1:h+h*(h-1)/2, 1:h, h*(h+1)/2+1:NN]); % small block first
        H1=P_ds*blkdiag(eye(h),kmat1(NN-h)')*P_alpha;
    case 'as'
        idx_as = zeros(1,NN);
        idx_as(1:h)=diag(fliplr(matI),0)';
        curI=h;
        for i=[-h+1:-1, 1:h-1]
            idx_in_zigzag = curI+(1:(h-abs(i)));
            if (i<0 && mod(i,2)==0) || (i>0 && mod(i,2)==1) 
                idx_as(fliplr(idx_in_zigzag)) = diag(fliplr(matI),i)';
            else
                idx_as(idx_in_zigzag) = diag(fliplr(matI),i)';
            end
            curI=curI+h-abs(i);
        end
        P_as = I_NN(:,idx_as);
        %P_alpha=I_NN(:,[1:h, h*(h+1)/2+1:NN, h+1:h+h*(h-1)/2]);  % large block first
        P_alpha=I_NN(:,[h+1:h+h*(h-1)/2, 1:h, h*(h+1)/2+1:NN]); % small block first
        H1=P_as*blkdiag(eye(h),kmat1(NN-h)')*P_alpha;
    case 'bs'
        P_beta=I_NN(:,[2*h+1:2*h+lB, h/2+1:h, 2*h+lB+1:2*h+2*lB, 1:h/2, ...
            2*h+2*lB+1:2*h+3*lB, h+1:2*h, 2*h+3*lB+1:NN]);
        H1=P_bs*blkdiag(kmat1(2*h)',kmat1(NN-2*h)')...
            *blkdiag(eye(2*h),kron(eye(2),kmat1(NN/2-h))')*P_beta;
    case 'penta'
        P_gamma=I_NN(:,[2*h+1:2*h+lB, h/2+1:h, 2*h+lB+1:2*h+2*lB, 1:h/2, ...
            2*h+2*lB+1:2*h+3*lB, h+1:h+h/2, 2*h+3*lB+1:2*h+3*lB+lB/2, ...
            h+h/2+1:2*h, 2*h+3*lB+lB/2+1:NN]);
        H1=P_bs*blkdiag(kmat1(2*h)',kmat1(NN-2*h)')...
            *blkdiag(eye(2*h),kron(eye(2),kmat1(NN/2-h)'))...
            *blkdiag(eye(h),kmat1(h)',kmat1(NN/4-h/2)',eye(NN/2-h),...
            kmat1(NN/4-h/2)')*P_gamma;
end

H=normc(H1);
