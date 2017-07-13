%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------------- Gradient Descent Calculation ------------------------%
% All the implementations are vectorized to speed up the program. --------%
%----------- Copyright (c) 2016, Shaghayegh Taheri Hosseinabadi ----------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function gradient_logp = Gradient_Descent(phi,F,alpha) 
 eps = 0.0000001;
 N = size(F,3); % N is the number of feature channels
 Mu1 = zeros(N,1);
 Mu2 = zeros(N,2);
 sigma1 = zeros(N,1);
 sigma2 = zeros(N,2);
  div = Kappa(phi);
 
  for n=1:N    
     [mu1,sigma1,mu2,sigma2] = Gaussian_Parameters(F(:,:,n),phi);
     Mu1(n) = mu1;
     Mu2(n) = mu2;
     sigma1(n) = sigma1;
     sigma2(n) = sigma2;
  end
  
 sum = zeros(size(phi));
 for n=1:N
 sum = sum + log2(normpdf(F(:,:,n),Mu1(n),sigma1(n)+eps)+eps) - log2(normpdf(F(:,:,n),Mu2(n),sigma2(n)+eps)+eps);
 end

gradient_logp = Heaviside_Derivative(phi).*(sum + alpha*div);

 
 
 function [Mu1,sigma1,Mu2,sigma2] = Gaussian_Parameters(F,phi) 
                                                                                                                             
    
Chi1 = Heaviside(phi);
Chi2 = 1 - Heaviside(phi);

num1 = F.*Chi1;
num2 = F.*Chi2; 
Mu1 = sum(num1(:))/sum(Chi1(:));
Mu2 = sum(num2(:))/sum(Chi2(:));

num11 = (F - Mu1).^2.*Chi1;
num22 = (F - Mu2).^2.*Chi2;

sigma1 = sqrt(sum(num11(:))/sum(Chi1(:)));
sigma2 = sqrt(sum(num22(:))/sum(Chi2(:)));
  



%Y = normpdf(X,mu,sigma)

function H = Heaviside(input)

H = 0.5*erf(input) + 0.5;
 
 
 
function H_p = Heaviside_Derivative(input)

H_p = exp(-input.^2)/sqrt(pi);



% function div = Divergence(image)
% % eps = 0.00001;
% % [Ix,Iy] = imgradientxy(image);
% % [Ixx,Ixy] = imgradientxy(Ix);
% % [Iyx,Iyy] = imgradientxy(Iy); 
% % num = Ixx.^2.*(Iy.^2) + Iyy.*(Ix.^2) - 2*Ix.*Iy.*Ixy;
% % den = (Ix.^2 + Iy.^2).^1.5;
% % div = num./(den+eps);
% 
% eps = 0.00001;
% u=image;
% [ux,uy]=gradient(u);
% G=sqrt(ux.^2+uy.^2);
% div =divergence(ux./(G+eps),uy./(G+eps));



 

