



function [gradient_phi] = Gradient_Desecent_NarrowBand(phi,F,alpha,pdf,w)
% F = feature vector 
% w : distance between narrow band and zero level set
% pdf is either 'Gaussian' or 'Parzen'.

 N = size(F,3); % N is the number of feature channels
 div = Kappa(phi);
 Mu1 = zeros(N,1);
 Mu2 = zeros(N,2);
 sigma1 = zeros(N,1);
 sigma2 = zeros(N,2);
 
  for n=1:N    
     [Mu1,sigma1,Mu2,sigma2] = Gaussian_Parameters(F(:,:,n),phi);
     Mu1(n) = Mu1;
     Mu2(n) = Mu2;
     sigma1(n) = sigma1;
     sigma2(n) = sigma2;
  end
 
  P1 = zeros(size(F,1),size(F,2));
  P2 = zeros(size(F,1),size(F,2));
  
sigmah = 4;
 sum = zeros(size(phi));
 for n=1:N
 
if strcmp(pdf,'Gaussian')==1
 P1 = normpdf(F(:,:,n),Mu1(n),sigma1(n)+eps);
 P2 = normpdf(F(:,:,n),Mu2(n),sigma2(n)+eps);
else if strcmp(pdf,'Parzen')==1
[P1,P2] = Parzen(F,phi,sigmah);
    end
end
%sum = sum + log2(normpdf(F(:,:,n),Mu1(n),sigma1(n)+eps)+eps) - log2(normpdf(F(:,:,n),Mu2(n),sigma2(n)+eps)+eps);
sum = sum + log2(P1+eps) - log2(P2+eps); 
 end
% gradient_phi = Heaviside_Derivative(phi).*(sum + alpha*div + beta*(f-1));
gradient_phi = Heaviside_Derivative(phi).*(sum + alpha*div);

 
 
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
  

 function [P1,P2] = Parzen(F,phi,sigmah) % Ps is a column vector of size 256 
     
P1 = zeros(size(F)); 
P2 = zeros(size(F));

Ps1 = zeros(256,1); 
Ps2 = zeros(256,1); 

Chi1 = Heaviside(phi);
Chi2 = 1 - Heaviside(phi);

for s = 0:255
delta_Fs = (floor(F) == s);    
num1 = delta_Fs.*Chi1;
num2 = delta_Fs.*Chi2;

Ps1(s+1) = sum(num1(:))/(sum(Chi1(:))+eps);
Ps2(s+1) = sum(num2(:))/(sum(Chi2(:))+eps);

kernel = Gkernel(sigmah);
Ps1 = conv(Ps1,kernel,'same');
Ps2 = conv(Ps2,kernel,'same');
end

Ps1 = Ps1/sum(Ps1);
Ps2 = Ps2/sum(Ps2);

for x = 1:size(F,1)
    for y = 1:size(F,2)
 P1(x,y) = Ps1(floor(F(x,y))+1); 
 P2(x,y) = Ps2(floor(F(x,y))+1); 
    end
end

function kernel = Gkernel(sigma)
    
size = sigma*2 + 1;    
mu = sigma + 1;
kernel = normpdf(1:size,mu,sigma+eps);
 kernel = kernel/sum(kernel);
   
 

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



 

