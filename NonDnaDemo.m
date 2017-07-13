%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------------- Demo for H&E images ------------------------------%
% All the implementations are vectorized to speed up the program. --------%
%----------- Copyright (c) 2016, Shaghayegh Taheri Hosseinabadi ----------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;
close all;
clear all;
clc;
 RGB = imread('H&E\\Benign\\Benign5_x40.tif');
[Hematoxylin,Eosin,DAB,RGB_Hematoxylin,RGB_Eosin,...
          RGB_DAB] = Deconvolve(RGB);     
image = 255-Hematoxylin;
image = double(imadjust(image)); %adjusting the contrast.
figure(); imshow(image,[]);

% image = Nonlinear_Diffusion(image,1,1,100);
% figure(); imshow(image,[]);

nrows = size(image,1);
ncols = size(image,2);
fmax = zeros(size(image));
as = zeros(size(image));
bs = zeros(size(image));
thetas = zeros(size(image));
range = 20;
%  range = [25 30];


for a = range
for b= [0.8 0.6]*a
for theta = 0:pi/8:pi
  f = frst2dEllipse(image,a,b,theta, 0.8,0.5,'dark'); 
  fmax0 = localMaximum(f,[min(range) min(range)],true,0.001); 
%   fmax0 = localMaximum(f,[16 16],true,0.001); 
  Index =find(fmax0>fmax);
  as(Index) = a;
  bs(Index)=b;
  thetas(Index)=theta;
  fmax(Index) = fmax0(Index);
end
end
end

fmax2 = localMaximum(fmax,[min(range) min(range)],true,0.001); 
Indexc=find(fmax2~=0);
ac = as(Indexc);
bc = bs(Indexc);
thetac= thetas(Indexc);
[xc,yc] = ind2sub(size(f),Indexc);


% figure(); imshow(RGB); hold on;
% plot(yc,xc,'g.');
% if (~isempty(xc))
% h=ellipse(ac,bc, pi/2 - thetac,yc,xc,'y');
% end


[phi,accpt_index] = InitialLevelSetEllipse(xc,yc,ac,bc,thetac,image,'dark');
figure(); imshow(RGB); hold on;
plot(yc(accpt_index),xc(accpt_index),'g.');
h=ellipse(ac(accpt_index),bc(accpt_index), pi/2 - thetac(accpt_index),yc(accpt_index),xc(accpt_index),'y');

phi = InitialLevelSetEllipse(xc,yc,ac*0.7,bc*0.7,thetac,image,'dark'); % last parameter depends on total image intensity.
phi=double((phi>0).*(bwdist(phi<0))-(phi<0).*(bwdist(phi>0)));
F(:,:,1)= double(image);


dt = 5;
pdf = 'Gaussian';
handle_seg =figure();
B= phi>=0;
W_nb = 3;
load TPLookUpTable.mat Simple %Load the look up table for checking topology preserving constraint

disp('Segmentation');
for iteration = 1:50 %120

alpha = 1;
% alpha = 50;
phi = Reinitialize_phi(phi);

figure(handle_seg); 
imshow(RGB,[]);
hold on;
contour(phi,[0 0],'y','LineWidth',1);
title(['iteration ',num2str(iteration)]);
drawnow;


phi_temp = ADIsolver(phi,F,alpha,dt);
Narrow_Band = (phi>-W_nb) & (phi<W_nb);
changed_sign = (sign(phi_temp)~=sign(phi)).*Narrow_Band;
unchanged_sign = (sign(phi_temp)==sign(phi)).*Narrow_Band;
index = find(changed_sign==1);
Height = size(image,1);
[px,py]= ind2sub(size(phi),index);
phi_new = phi_temp;
X = B.*(phi<W_nb);
Xbar = mod(B+1,2).*(phi>-W_nb);


for i=1:size(index)
 
n = LT(index(i),X,Height); % Look Up table value 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%% Using Look Up Table, instead of direct computation is 100 times
 %%%% faster than the regular method:
 %%%% T4X = TopologicalNumbers(index(i),find(X==1),Height,4);
 %%%% T8Xbar = TopologicalNumbers(index(i),find(Xbar==1),Height,8);
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  finding the coordinates of level set 
%  pixels that are changing but are 
%  topologically simple:
% if(T4X ~=1 || T8Xbar~=1)  % Regular Method  
 if (Simple(n+1)==0)          % Look Up Table
   phi_new(px(i),py(i)) = 0.001*sign(phi(px(i),py(i)));
 else
  B(px(i),py(i))=mod(B(px(i),py(i))+1,2);   %Applying only on simple points with sign change
  X(px(i),py(i))=mod(X(px(i),py(i))+1,2);
  Xbar(px(i),py(i))=mod(Xbar(px(i),py(i))+1,2);
 end
end

phi = phi_new; 
end


%%% Post Processing Operation (Morphology Opening)
phi=double((phi>0).*(bwdist(phi<0))-(phi<0).*(bwdist(phi>0)));
BW = phi>0.5;
BW = bwareaopen(BW,3);
% se = strel('disk',3); 
% BW = imopen(BW,se);
% figure(); imshow(BW);

% BW = im2bw(phi);
% Segments = bwlabel(BW,4);
% RGBSegments = label2rgb(Segments,'jet',[1 1 1],'shuffle');
% figure(); 
% imshow(RGBSegments);


%%% Final Results %%%
figure(); imshow(RGB,[]);
hold on;
contour(phi,[0 0],'g','LineWidth',2);
title('Proposed Method');

Segments = bwlabel(BW,4);
RGBSegments = label2rgb(Segments,'jet',[1 1 1],'shuffle');
imwrite(RGBSegments,'MyMethod.png');
figure(); 
imshow(RGBSegments);

imwrite(BW,'MyMethod2.png');

% CP = imread('CP.png');
% CP = im2bw(CP);
% Segments = bwlabel(CP,4);
% CPSegments = label2rgb(Segments,'jet',[1 1 1],'shuffle');
% imwrite(CPSegments,'CP.png');
% figure(); 
% imshow(RGBSegments);

toc;







