function K=kmat(N)
% K matrix (a permuted butterfly)
% K=kmat(N)
% N: size, must be even
% K: butterfly matrix

if mod(N,2)==1
    h=(N-1)/2;
    K=1/sqrt(2)*[    eye(h), zeros(h,1),   -jmat(h); 
                 zeros(1,h),    sqrt(2), zeros(1,h);
                     eye(h), zeros(h,1),    jmat(h)];
else
    K=1/sqrt(2)*[eye(N/2), -jmat(N/2); eye(N/2), jmat(N/2)];
end