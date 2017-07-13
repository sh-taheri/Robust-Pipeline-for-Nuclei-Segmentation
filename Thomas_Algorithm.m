%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------------------- Thomas Algorithm ---------------------------%
% All the implementations are vectorized to speed up the program. --------%
% Solving the equation Bu=d , where B is a N*N tridiagonal matrix with the
% diagonal elements alpha(1),...,alpha(N) and the upper diagonal elements
% beta(1),...,beta(N-1) and the lower diagonal elements
% gamma(1),...,gamma(N-1).
%----------- Copyright (c) 2016, Shaghayegh Taheri Hosseinabadi ----------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function u = Thomas_Algorithm(alpha,beta,gamma,d) 

N = size(alpha,1);
m = zeros(N,1);
l = zeros(N-1,1);
y = zeros(N,1);
u = zeros(N,1);
% ---------------------- Step 1 : LR Decomposition ---------------------- %
m(1) = alpha(1);
for i=1:N-1
    l(i) = gamma(i)/m(i);
    m(i+1)= alpha(i+1) - l(i)*beta(i);
 end
% ---------------- Step 2 : Forward Substitution (Ly = d) --------------- %
y(1) = d(1);
 for i=2:N
    y(i) = d(i) - l(i-1)*y(i-1);
 end
% --------------- Step 3 : Backward Substitution (Ru = y) --------------- %
u(N) = y(N)/m(N);
 for i=N-1:-1:1
    u(i) = (y(i) - beta(i)*u(i+1))/m(i);
 end
