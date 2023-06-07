function Xbs=bisym(X)
% make bisymmetric by truncation

Xs = (X+X')/2;
Xbs = (Xs+rot90(Xs,2))/2;