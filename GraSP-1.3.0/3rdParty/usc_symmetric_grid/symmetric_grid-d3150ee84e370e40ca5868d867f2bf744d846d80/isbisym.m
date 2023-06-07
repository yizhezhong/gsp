function b=isbisym(M)
er=abs(M-bisym(M));
b=max(er(:))<1e-8;