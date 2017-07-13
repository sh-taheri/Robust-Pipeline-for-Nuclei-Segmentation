%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------- ADI Solver (Semi-Implicit Implementation)---------------%
% Check thesis S. Gao, "A new image segmentation and smoothing method ----%
% based on the mumford-shah variational model", Master of Computer Science%
% Thesis, Concordia University, 2003.-------------------------------------% 
% All the implementations are vectorized to speed up the program. --------%
%----------- Copyright (c) 2016, Shaghayegh Taheri Hosseinabadi ----------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [uNew,alphaB,betaB,gammaB,rx1,rx2,ry1,ry2,D] = ADIsolver(u,F,alpha,dt) % F = feature vector
gradient_logp = Gradient_LogP(u,F); 
% gradient_chanvese = Gradient_ChanVese(u,F);
h=1;
epsilon = h;
[nx,ny] = size(u);
uNew = u;
utemp = u;

rx1 = zeros(nx,ny);
rx2 = zeros(nx,ny);
ry1 = zeros(nx,ny);
ry2 = zeros(nx,ny);


% for j = 2:ny-1
%      for i= 2:nx-1
% 
%         
%         C1 = 1/sqrt( (u(i+1,j)-u(i,j))^2/(h^2) + (u(i,j+1)-u(i,j-1))^2/(4*h^2) + epsilon);
%         C2 = 1/sqrt( (u(i,j)-u(i-1,j))^2/(h^2) + (u(i-1,j+1)-u(i-1,j-1))^2/(4*h^2) + epsilon);
%         C3 = 1/sqrt( (u(i+1,j)-u(i-1,j))^2/(4*h^2) + (u(i,j+1)-u(i,j))^2/(h^2) + epsilon);
%         C4 = 1/sqrt( (u(i+1,j-1)-u(i-1,j-1))^2/(4*h^2) + (u(i,j)-u(i,j-1))^2/(h^2) + epsilon);                
% %       D = (1/pi)* epsilon/(epsilon^2 + u(i,j)^2); 
%         D = Heaviside_Derivative(u(i,j));
%         rx1(i,j) = alpha*dt*D*C1/(h^2);
%         rx2(i,j) = alpha*dt*D*C2/(h^2);
%         ry1(i,j) = alpha*dt*D*C3/(h^2);
%         ry2(i,j) = alpha*dt*D*C4/(h^2);  
%         
%     end
% end

i=2:nx-1;
j=2:ny-1;
C1 = 1./sqrt( (u(i+1,j)-u(i,j)).^2/(h^2) + (u(i,j+1)-u(i,j-1)).^2/(4*h^2) + epsilon);
C2 = 1./sqrt( (u(i,j)-u(i-1,j)).^2/(h^2) + (u(i-1,j+1)-u(i-1,j-1)).^2/(4*h^2) + epsilon);
C3 = 1./sqrt( (u(i+1,j)-u(i-1,j)).^2/(4*h^2) + (u(i,j+1)-u(i,j)).^2/(h^2) + epsilon);
C4 = 1./sqrt( (u(i+1,j-1)-u(i-1,j-1)).^2/(4*h^2) + (u(i,j)-u(i,j-1)).^2/(h^2) + epsilon); 
D = Heaviside_Derivative(u(i,j));
rx1(i,j) = alpha*dt*D.*C1/(h^2);
rx2(i,j) = alpha*dt*D.*C2/(h^2);
ry1(i,j) = alpha*dt*D.*C3/(h^2);
ry2(i,j) = alpha*dt*D.*C4/(h^2); 


DT = Gradient(u,gradient_logp,dt);
% DT = Gradient(u,gradient_chanvese,dt); 

for j=1:ny
% j=1:ny;
rhSideB = getrhsB(j,ry1,ry2,DT,u);
[alphaB,betaB,gammaB] = getTridiagB(j,rx1,rx2,nx);
utemp(:,j) = Thomas_Algorithm(alphaB,betaB,gammaB,rhSideB); % utemp : u^(n + 1/2)
end

% for i= 2:nx-1
%     for j = 2:ny-1
%         
%         C1 = 1/sqrt( (utemp(i+1,j)-utemp(i,j))^2/(h^2) + (utemp(i,j+1)-utemp(i,j-1))^2/(4*h^2) + epsilon);
%         C2 = 1/sqrt( (utemp(i,j)-utemp(i-1,j))^2/(h^2) + (utemp(i-1,j+1)-utemp(i-1,j-1))^2/(4*h^2) + epsilon);
%         C3 = 1/sqrt( (utemp(i+1,j)-utemp(i-1,j))^2/(4*h^2) + (utemp(i,j+1)-utemp(i,j))^2/(h^2) + epsilon);
%         C4 = 1/sqrt( (utemp(i+1,j-1)-utemp(i-1,j-1))^2/(4*h^2) + (utemp(i,j)-utemp(i,j-1))^2/(h^2) + epsilon);                
%         D = (1/pi)* epsilon/(epsilon^2 + utemp(i,j)^2);       
%         rx1(i,j) = alpha*dt*D*C1/(h^2);
%         rx2(i,j) = alpha*dt*D*C2/(h^2);
%         ry1(i,j) = alpha*dt*D*C3/(h^2);
%         ry2(i,j) = alpha*dt*D*C4/(h^2);  
%         
%     end
% end
% DT = Gradient(utemp,gradient_logp,dt);

for i=1:nx
rhSideQ = getrhsQ(i,rx1,rx2,DT,utemp);
[alphaQ,betaQ,gammaQ] = getTridiagQ(i,ry1,ry2,ny);
ui = Thomas_Algorithm(alphaQ,betaQ,gammaQ,rhSideQ); % utemp : u^(n + 1/2)
uNew(i,:) = ui';
end

    
function [alphaB,betaB,gammaB] = getTridiagB(j,rx1,rx2,nx)
% Constructing Tridiagonal Matrix elements which are three vectors:
% alpha (n by 1) , beta (n-1 by 1) , gamma (n-1 by 1).

alphaB = zeros(nx,1);
betaB = zeros(nx-1,1);
gammaB = zeros(nx-1,1);
% for i=2:nx-2
%     alphaB(i) = 1 + (rx1(i,j)+rx2(i,j))/2;
%     betaB(i) = -rx1(i,j)/2;
%     gammaB(i) = -rx2(i+1,j)/2;
% end
i=2:nx-2;
alphaB(i) = 1 + (rx1(i,j)+rx2(i,j))/2;
betaB(i) = -rx1(i,j)/2;
gammaB(i) = -rx2(i+1,j)/2;

    alphaB(1) = 1 + (rx1(1,j)+rx2(1,j))/2;
    alphaB(nx-1) = 1 + (rx1(nx-1,j)+rx2(nx-1,j))/2;
    alphaB(nx) = 1 + (rx1(nx,j)+rx2(nx,j))/2;
    betaB(1) = -(rx1(1,j)+rx2(1,j))/2;
    betaB(nx-1) = -rx1(nx-1,j)/2;
    gammaB(1) = -rx2(2,j)/2;
    gammaB(nx-1) = -(rx1(nx,j)+rx2(nx,j))/2;

 


function [alphaQ,betaQ,gammaQ] = getTridiagQ(i,ry1,ry2,ny)
% Constructing Tridiagonal Matrix elements which are three vectors:
% alpha (n by 1) , beta (n-1 by 1) , gamma (n-1 by 1).

alphaQ = zeros(ny,1);
betaQ = zeros(ny-1,1);
gammaQ = zeros(ny-1,1);
% for j=2:ny-2
%     alphaQ(j) = 1 + (ry1(i,j)+ry2(i,j))/2;
%     betaQ(j) = -ry1(i,j)/2;
%     gammaQ(j) = -ry2(i,j+1)/2;
% end
j=2:ny-2;
alphaQ(j) = 1 + (ry1(i,j)+ry2(i,j))/2;
betaQ(j) = -ry1(i,j)/2;
gammaQ(j) = -ry2(i,j+1)/2;

    alphaQ(1) = 1 + (ry1(i,1)+ry2(i,1))/2;
    alphaQ(ny-1) = 1 + (ry1(i,ny-1)+ry2(i,ny-1))/2;
    alphaQ(ny) = 1 + (ry1(i,ny)+ry2(i,ny))/2;
    betaQ(1) = -(ry1(i,1)+ry2(i,1))/2;
    betaQ(ny-1) = -ry1(i,ny-1)/2;
    gammaQ(1) = -ry2(i,2)/2;
    gammaQ(ny-1) = -(ry1(i,ny)+ry2(i,ny))/2;



function rhSideB = getrhsB(j,ry1,ry2,DT,u)
[nx,ny] = size(u);
rhSideB = zeros(nx,1);

    if j==1
%         for i = 1:nx
%         rhSideB(i)=(ry1(i,j)+ry2(i,j))/2*u(i,2) + ...
%         (1 - (ry1(i,j)+ry2(i,j))/2)*u(i,1) + DT(i,j);
%          end
          i = 1:nx;
          rhSideB(i)=(ry1(i,j)+ry2(i,j))/2.*u(i,2) + ...
          (1 - (ry1(i,j)+ry2(i,j))/2).*u(i,1) + DT(i,j);
    else if j== ny
%          for i = 1:nx
%         rhSideB(i)= (ry1(i,j)+ry2(i,j))/2*u(i,ny-1) + ...
%         (1 - (ry1(i,j)+ry2(i,j))/2)*u(i,ny) + DT(i,j);
%          end
        i = 1:nx;
        rhSideB(i)= (ry1(i,j)+ry2(i,j))/2.*u(i,ny-1) + ...
        (1 - (ry1(i,j)+ry2(i,j))/2).*u(i,ny) + DT(i,j);
        else 
%         for i = 1:nx
%         rhSideB(i) = ry1(i,j)/2*u(i,j+1) + (1 - (ry1(i,j)+ry2(i,j))/2)*u(i,j)...
%         + ry2(i,j)/2*u(i,j-1) + DT(i,j);
%         end
        i = 1:nx;
        rhSideB(i) = ry1(i,j)/2.*u(i,j+1) + (1 - (ry1(i,j)+ry2(i,j))/2).*u(i,j)...
        + ry2(i,j)/2.*u(i,j-1) + DT(i,j);
        end
    end

        

function rhSideQ = getrhsQ(i,rx1,rx2,DT,u)
[nx,ny] = size(u);
rhSideQ = zeros(ny,1);

    if i==1
%         for j = 1:ny
%         rhSideQ(j)=(rx1(i,j)+rx2(i,j))/2*u(2,j) + ...
%         (1 - (rx1(i,j)+rx2(i,j))/2)*u(1,j) + DT(i,j);
%         end
          j = 1:ny;
          rhSideQ(j)=(rx1(i,j)+rx2(i,j))/2.*u(2,j) + ...
          (1 - (rx1(i,j)+rx2(i,j))/2).*u(1,j) + DT(i,j);
    else if i== nx
%          for j = 1:ny
%         rhSideQ(j)= (rx1(i,j)+rx2(i,j))/2*u(nx-1,j) + ...
%         (1 - (rx1(i,j)+rx2(i,j))/2)*u(nx,j) + DT(i,j);
%          end
           j = 1:ny;
           rhSideQ(j)= (rx1(i,j)+rx2(i,j))/2.*u(nx-1,j) + ...
           (1 - (rx1(i,j)+rx2(i,j))/2).*u(nx,j) + DT(i,j);
        else 
%         for j = 1:ny
%         rhSideQ(j) = rx1(i,j)/2*u(i+1,j) + (1 - (rx1(i,j)+rx2(i,j))/2)*u(i,j)...
%         + rx2(i,j)/2*u(i-1,j) + DT(i,j);
%         end
        j = 1:ny;
        rhSideQ(j) = rx1(i,j)/2.*u(i+1,j) + (1 - (rx1(i,j)+rx2(i,j))/2).*u(i,j)...
        + rx2(i,j)/2.*u(i-1,j) + DT(i,j);
        end
    end
    
    function DT = Gradient(u,gradient_logp,dt)
%     epsilon = 1;    
%     D = (1/pi)* epsilon./(epsilon^2 + u.^2);
    D = Heaviside_Derivative(u);
    DT = (dt/2).*D.*(gradient_logp);
