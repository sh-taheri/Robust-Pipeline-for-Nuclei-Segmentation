%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----------------- FTCS Solver (Explicit Implementation)------------------%
% Check thesis S. Gao, "A new image segmentation and smoothing method ----%
% based on the mumford-shah variational model", Master of Computer Science%
% Thesis, Concordia University, 2003.-------------------------------------% 
% All the implementations are vectorized to speed up the program. --------%
%----------- Copyright (c) 2016, Shaghayegh Taheri Hosseinabadi ----------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function unew = FTCSsolver(u,F,alpha,dt) % F = feature vector

%%% change the boundary condition such that x-1 = x1 and...!!!! (Later)
gradient_logp = Gradient_LogP(u,F); 
h=1;
epsilon = h;
[nx,ny] = size(u);
unew = u;
for i= 2:nx-1
    for j = 2:ny-1
        
        C1 = 1/sqrt( (u(i+1,j)-u(i,j))^2/(h^2) + (u(i,j+1)-u(i,j-1))^2/(4*h^2) + epsilon);
        C2 = 1/sqrt( (u(i,j)-u(i-1,j))^2/(h^2) + (u(i-1,j+1)-u(i-1,j-1))^2/(4*h^2) + epsilon);
        C3 = 1/sqrt( (u(i+1,j)-u(i-1,j))^2/(4*h^2) + (u(i,j+1)-u(i,j))^2/(h^2) + epsilon);
        C4 = 1/sqrt( (u(i+1,j-1)-u(i-1,j-1))^2/(4*h^2) + (u(i,j)-u(i,j-1))^2/(h^2) + epsilon);        
        D = (1/pi)* epsilon/(epsilon^2 + u(i,j)^2);
%       D = Heaviside_Derivative(u(i,j));
        m = alpha*dt*D/(h^2);        
        C = 1 + m*(C1 + C2 + C3 + C4);        
        unew(i,j) = (1/C)*( u(i,j) + m*(C1*u(i+1,j) + C2*u(i-1,j) + C3*u(i,j+1) + C4*u(i,j-1))...
                    + dt*D*gradient_logp(i,j));
    end
end

