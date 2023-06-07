function issym=check_grid_sym(L,symtype)
% issym=check_grid_sym(L, symtype)
%    check whether the matrix L is a Laplacian matrix of a grid having 
%    input symmetry type. The current version only support square grids.
%    
% L: Laplacian matrix of a grid (with column-first vertex ordering), 
%    whose size must be N^2xN^2 with some integer N. 
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
% by: KS Lu
% last update: 20170105

NN=size(L,1);
N=sqrt(NN);
L_offdiag = L(:); L_offdiag(1:NN+1:end)=0;
if size(L,1)~=size(L,2) || round(N)~=N || any(L_offdiag>1e-4) ...
        || any(diag(L)<-1e-4) || any(sum(L)<-1e-4) || ~issymmetric(L)
    error('L is not a valid Laplacian matrix for a grid');
end

issym=true;
matI=reshape(1:NN,N,N);
I_NN=eye(NN);

% here we check each symmetry type
% note that: although pentasymmetric grids satisfy all 5 symmetries, 
% it suffices to assert pentasymmetry by showing this grid is both 
% UD-symmetric and diagonally symmetric. Thus, the following would be 
% sufficient

% check centrosymmetry
if ismember(symtype, {'centro'})
    issym = issym & isbisym(L);
end

% check UD-symmetry
if ismember(symtype, {'ud', 'udlr', 'penta'})
    % row-first zigzag
    idx_h = matI'; idx_h(:,2:2:end)=idx_h(N:-1:1,2:2:end);
    P_h = I_NN(:,idx_h);
    L_h = P_h'*L*P_h;
    issym = issym & isbisym(L_h);
end

% check LR-symmetry
if ismember(symtype, {'lr', 'udlr'})
    % column-first zigzag
    idx_v = matI; idx_v(:,2:2:end)=idx_v(N:-1:1,2:2:end);
    P_v = I_NN(:,idx_v);
    L_v = P_v'*L*P_v;
    issym = issym & isbisym(L_v);
end

% check diagonal symmetry
if ismember(symtype, {'ds', 'bs', 'penta'})
    idx_ds = zeros(1,NN);
    idx_ds(1:N)=diag(matI,0)';
    curI=N;
    for i=[-N+1:-1, 1:N-1]
        idx_in_zigzag = curI+(1:(N-abs(i)));
        if (i<0 && mod(i,2)==0) || (i>0 && mod(i,2)==1)  % reverse order
            idx_ds(fliplr(idx_in_zigzag)) = diag(matI,i)';
        else
            idx_ds(idx_in_zigzag) = diag(matI,i)';
        end
        curI=curI+N-abs(i);
    end
    P_ds = I_NN(:,idx_ds);
    L_ds = P_ds'*L*P_ds;
    % check B & BJ
    b=L_ds(1:N, N+1:NN); 
    b(:,2:2:end)=b(N:-1:1,2:2:end);
    Nb=numel(b);
    erB=abs(b(1:Nb/2)-b(Nb:-1:Nb/2+1));
    % check C & JCJ and F & JFJ
    cf=L_ds(N+1:NN, N+1:NN); Ncf=numel(cf);
    erCF=abs(cf(1:Ncf/2)-cf(Ncf:-1:Ncf/2+1));
    issym = issym & max(erB)<1e-8 & max(erCF)<1e-8;
end

% check anti-diagonal symmetry
if ismember(symtype, {'as', 'bs'})
    idx_as = zeros(1,NN);
    idx_as(1:N)=diag(fliplr(matI),0)';
    curI=N;
    for i=[-N+1:-1, 1:N-1]
        idx_in_zigzag = curI+(1:(N-abs(i)));
        if (i<0 && mod(i,2)==0) || (i>0 && mod(i,2)==1) 
            idx_as(fliplr(idx_in_zigzag)) = diag(fliplr(matI),i)';
        else
            idx_as(idx_in_zigzag) = diag(fliplr(matI),i)';
        end
        curI=curI+N-abs(i);
    end
    P_as = I_NN(:,idx_as);
    L_as = P_as'*L*P_as;
    % check B & BJ
    b=L_as(1:N, N+1:NN); 
    b(:,2:2:end)=b(N:-1:1,2:2:end);
    Nb=numel(b);
    erB=abs(b(1:Nb/2)-b(Nb:-1:Nb/2+1));
    % check C & JCJ and F & JFJ
    cf=L_as(N+1:NN, N+1:NN); Ncf=numel(cf);
    erCF=abs(cf(1:Ncf/2)-cf(Ncf:-1:Ncf/2+1));
    issym = issym & max(erB)<1e-8 & max(erCF)<1e-8;
end
