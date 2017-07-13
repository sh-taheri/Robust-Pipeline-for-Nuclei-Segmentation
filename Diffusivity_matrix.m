
%----------------------------------------------------------------------%
%--------------------- Non-Linear Diffusion ---------------------------%
%----------------------------------------------------------------------%

function [alpha, beta]= Diffusivity_matrix(u,h,tau,diffusivity_type)

n=size(u,1);
alpha = zeros(n,1);
beta = zeros(n-1,1);

if (strcmp(diffusivity_type,'nonlinear')==1)
    diffusivity = Diffusivity(u,h);
else if (strcmp(diffusivity_type,'linear')==1)
        diffusivity = ones(n,1);
    else
            disp('Wrong parameter');
    return     
    end
end

for i=1:n-1
    beta(i) = -4*tau* ((diffusivity(i) + diffusivity(i+1))/(2*h^2));
    if(i>1)
    alpha(i) = 2 - 4*tau* ((-2*diffusivity(i) - diffusivity(i-1) - ...
    diffusivity(i+1))/(2*h^2)); 
    end
    %comment the coefficient 2 for diffusivity[i] like the paper
end

alpha(n) = 2 - 4*tau* ((-2*diffusivity(n) - diffusivity(n-1))/(2*h^2));
alpha(1) = 2 - 4*tau* ((-2*diffusivity(1) - diffusivity(2))/(2*h^2));
     




function diffusivity = Diffusivity(u,h)

n = size(u,1); % n = size of the vector
m = size(u,2);
diffusivity = zeros(n,1);

for i= 2:n-1        
      sum = 0;
      for j= 1:m
      sum = sum + (u(i-1,j) - u(i+1,j))^2;
      end  
      diffusivity(i) = g( 0.5*sum/((2*h)^2) );
end
sum1 = 0;
sumn = 0;
for j=1:m
    sum1 = sum1 + (u(1,j) - u(2,j))^2;
    sumn = sumn + (u(n,j) - u(n-1,j))^2;
diffusivity(1) =  g(0.5*sum1/((2*h)^2)); % Reflecting boundary
diffusivity(n) =  g(0.5*sumn/((2*h)^2)); % Reflecting boundary
end


function out= g(s)
eps = 0.01;
p = 1.6; % As assumed in the paper (1)
out = 1/((s + eps^2)^(p/2));

% ----------------------------------------------------------------------
% [1] Colour,texture and motion in level set based segmentation and
% tracking


