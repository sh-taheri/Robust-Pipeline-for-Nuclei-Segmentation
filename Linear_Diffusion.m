%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------------------------ Linear Diffusion -------------------------%
%----------- Copyright (c) 2016, Shaghayegh Taheri Hosseinabadi ----------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function output = Linear_Diffusion(input,iteration,sigma)

output = zeros(size(input));
row = size(input,1);
col = size(input,2);
n = row*col;
h=1;
tau = sigma^2/2;

u_x = im2vec_x(input);
u_y = im2vec_y(input);

for t=1:iteration    
[alpha_x,beta_x] = Diffusivity_matrix(u_x,h,tau,'linear');
[alpha_y,beta_y] = Diffusivity_matrix(u_y,h,tau,'linear');
v_x = Thomas_Algorithm(alpha_x,beta_x,beta_x,u_x); 
v_y = Thomas_Algorithm(alpha_y,beta_y,beta_y,u_y); 

for i=1:n
    u_x(i) =  v_x(i)  + v_y( mod(i-1,row)*col + floor((i-1)/row) +1); 
    u_y(i) =  v_y(i)  + v_x( mod(i-1,col)*row + floor((i-1)/col) +1);	
end
end

for i=1:row 
   for j=1:col
       output(i,j) = u_y((i-1)*col+j);
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




