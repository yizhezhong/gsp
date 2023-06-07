function K=kmat1(N)
% K matrix with only 1, 0, and -1 (it is not orthogonal, but does not 
% require multiplications
% K=kmat1(N)
% N: size

if mod(N,2)==1
    h=(N-1)/2;
    K=[    eye(h), zeros(h,1),   -jmat(h); 
       zeros(1,h),    sqrt(2), zeros(1,h);
           eye(h), zeros(h,1),    jmat(h)];
else
    K=[eye(N/2), -jmat(N/2); eye(N/2), jmat(N/2)];
end