%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-- Nonlinear Diffusion Based on AOS Scheme(Semi-Implicit Implementation)-%
% Check paper Weickert, Joachim, BM Ter Haar Romeny, and Max A. Viergever.
% "Efficient and reliable schemes for nonlinear diffusion filtering."
% (1998): 398-410.
% All the implementations are vectorized to speed up the program. --------%
%----------- Copyright (c) 2016, Shaghayegh Taheri Hosseinabadi ----------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function output = Nonlinear_Diffusion(input,iteration,sigma,tau)

number = size(input,3); % number of images for diffusion
output = zeros(size(input));
row = size(input,1);
col = size(input,2);
n = row*col;
h=1;
tau0 = sigma^2/2;

u_x = im2vec_x(input);
u_y = im2vec_y(input);

v_x = zeros(n,number);
v_y = zeros(n,number);

for t=1:iteration 
   
for k = 1:number
         
[alpha_x,beta_x] = Diffusivity_matrix(u_x(:,k),h,tau0,'linear');
[alpha_y,beta_y] = Diffusivity_matrix(u_y(:,k),h,tau0,'linear');
v_x(:,k) = Thomas_Algorithm(alpha_x,beta_x,beta_x,u_x(:,k)); 
v_y(:,k) = Thomas_Algorithm(alpha_y,beta_y,beta_y,u_y(:,k)); 

for i=1:n
    u_x(i,k) =  v_x(i,k)  + v_y( mod(i-1,row)*col + floor((i-1)/row) +1 ,k); 
    u_y(i,k) =  v_y(i,k)  + v_x( mod(i-1,col)*row + floor((i-1)/col) +1 ,k);	
end

[alpha_x,beta_x] = Diffusivity_matrix(u_x,h,tau,'nonlinear');
[alpha_y,beta_y] = Diffusivity_matrix(u_y,h,tau,'nonlinear');
v_x(:,k) = Thomas_Algorithm(alpha_x,beta_x,beta_x,u_x(:,k)); 
v_y(:,k) = Thomas_Algorithm(alpha_y,beta_y,beta_y,u_y(:,k)); 

end

for k=1:number
    for i=1:n
    u_x(i,k) =  v_x(i,k)  + v_y( mod(i-1,row)*col + floor((i-1)/row) +1 ,k); 
    u_y(i,k) =  v_y(i,k)  + v_x( mod(i-1,col)*row + floor((i-1)/col) +1 ,k);	
    end
end

end

for k = 1:number
for i=1:row 
   for j=1:col
       output(i,j,k) = u_y((i-1)*col+j,k);
   end
end
end





function vector_x = im2vec_x(image)

row = size(image,1);
col = size(image,2);
number = size(image,3); % number of images stored in image variable
n=row*col;
vector_x = zeros(n,number);

for k = 1:number
for i=1:row
    for j=1:col
        
        vector_x( (j-1)*row + i ,k) = image(i,j,k);
    end
end
end




function vector_y = im2vec_y(image)

row = size(image,1);
col = size(image,2);
number = size(image,3); % number of images stored in image variable
n=row*col;
vector_y = zeros(n,number);

for k = 1:number
for i=1:row
    for j=1:col
        
        vector_y( (i-1)*col + j,k) = image(i,j,k);
    end
end
end










