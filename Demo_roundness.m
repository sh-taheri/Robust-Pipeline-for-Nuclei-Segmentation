
%%% Removing Topology Constraint 

clc;
close all;
clear all;

RGB = imread('cells4.png');
% RGB = imread('test3.png');
% RGB = imread('B.png');
figure(); imshow(RGB); title('Original Image');
% RGB = imresize(RGB,0.3);

[Hematoxylin,Eosin,DAB,RGB_Hematoxylin,RGB_Eosin,...
          RGB_DAB] = Deconvolve(RGB);     
% figure();imshow(RGB);
% figure();imshow(RGB_Hematoxylin);
% figure();imshow(RGB_Eosin);
% figure();imshow(RGB_DAB);

figure();imshow(255-Hematoxylin);
image = double(255-Hematoxylin); 
image = image.^2/256; % Improve the contrast
figure(); imshow(uint8(image)); title('Color Decomposition');


% n = 15:5:35; % cells6.png
% n = 40:5:90; %Shape2.png
n = 10:5:30; % for cells4.png
% n = 20:5:30; %test3.png
% n = 5:2:15; %B.png
% n = 5:2:20; cells5.png ?
% n = 30:5:50

[f] = frst2d(image, n,0.8, 0.5, 'dark');
imagesc(f);
colormap jet;
title('FRST Transform');
% [f] = frst2d(image, n,0.8, 0.5, 'bright');
% [x y z] = localMaximum(f,[9 9],true); %cells4.png
  fmax0 = localMaximum(f,[max(range) min(range)],true,0.001); 
  Index =find(fmax0>fmax);
% [x y z] = localMaximum(f,[5 5],true); cells4.png
count = 1;
for i=1:size(x,1)
    if (f(x(i),y(i))<0.08) % test3.png
    %if (f(x(i),y(i))<0.1) %B,cells4 , ...
        x(i)=NaN;
        y(i)= NaN;
    end
end

x = x(5);
y = y(5);

figure(); imshow(RGB); hold on; title('Detecting Center Points with FRST');
plot(y,x,'g.');


%%% Comment %%%
phi = InitialLevelSet(x,y,10,size(image,2),size(image,1)); %cells4.png
% phi = InitialLevelSet(x,y,5,size(image,2),size(image,1)); %B.png
% phi = InitialLevelSet(x,y,10,size(image,2),size(image,1)); %B.png

phi=double((phi>0).*(bwdist(phi<0))-(phi<0).*(bwdist(phi>0)));

 F(:,:,1)= image;


% dt = 10; for ADIsolver
% alpha = 5; for ADIsolver
x(isnan(x)) = [];
y(isnan(y)) = [];


xc = zeros(size(x));
yc = zeros(size(y));


roundness = zeros(size(phi));
% for i=1:size(phi,1)
%     for j=1:size(phi,2)
%        label = L(i,j);
%        if label~=0
%         roundness(i,j)=1/sqrt((i-xc(label))^2 + (j-yc(label))^2);
%        end
%     end
% end


phi_bw = phi>0;
[L,num] = bwlabel(phi_bw);
% xc = zeros(num,1);
% yc = zeros(num,1);
% for label=1:num
% [r,c] = find(L==label);
%  xc(label) = mean(r);
%  yc(label) = mean(c);
% end
 

figure(); 
imshow(RGB,[]);
hold on;
contour(phi,[0 0],'b','LineWidth',2);
title('Initial Contours');

dt = 1;
% alpha = 0.8; %B.png
% beta = 3; %B.png
% gamma = 0.1; %B.png

% alpha = 0.8; %cells4.png
% beta = 10; %cells4.png
% gamma = 0.3; %cells4.png 

alpha = 0.8; %cells4.png
beta = 0; %cells4.png
gamma = 0.5; %cells4.png 


pdf = 'Gaussian';
figure();
for iteration = 1:150
    
phi = Reinitialize_phi(phi);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phi_bw = phi>0;
% phi_bw = im2bw(phi);
[L,num] = bwlabel(phi_bw,4);
for label=1:num
[r,c] = find(L==label);
 xc(label) = mean(r);
 yc(label) = mean(c);
end

% Label = ordfilt2(L,25,true(5)); % Creating the outer narrow band
Label = ordfilt2(L,25,true(5));
for i=1:size(phi,1)
    for j=1:size(phi,2)
        label = Label(i,j);
        if (label ~=0)
        roundness(i,j)=1/sqrt((i-xc(label))^2 + (j-yc(label))^2 + eps);
       end
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% gradient_phi = Gradient_Descent(phi,F,alpha);
gradient_phi = gamma*Gradient_Descent(phi,F,alpha)+ beta*(roundness).*Heaviside_Derivative(phi);
%gradient_phi = beta*(roundness).*Heaviside_Derivative(phi);


phi = phi + dt*gradient_phi;
% phi_temp = ADIsolver(phi,F,alpha,dt);
imshow(RGB,[]);
hold on;
contour(phi,[0 0],'b','LineWidth',2);
plot(yc,xc,'g.');
title(['iteration ',num2str(iteration)]);
drawnow;
end


BW = im2bw(phi);
Segments = bwlabel(BW,4);
RGBSegments = label2rgb(Segments);
figure(); 
imshow(RGBSegments);








 
