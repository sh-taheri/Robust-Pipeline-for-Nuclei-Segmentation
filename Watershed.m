%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------ Sample Code For Applying Marker Controlled watershed MATLAB ------%
%----------- Copyright (c) 2016, Shaghayegh Taheri Hosseinabadi ----------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;
clear all;
clc;

RGB = imread('Hodgkin_lymphoma_(3)_mixed_cellularity_type.jpg');
image = RGB(:,:,3);
image = double(imadjust(image)); %adjusting the contrast.
figure(); imshow(image,[]);
title('Original image after color deconvolution');



nrows = size(image,1);
ncols = size(image,2);
fmax = zeros(size(image));
as = zeros(size(image));
bs = zeros(size(image));
thetas = zeros(size(image));
range = [8 12];
% range = 20;
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



figure(); imshow(fmax,[]);
fmax2 = localMaximum(fmax,[min(range) min(range)],true,0.001); 
Indexc=find(fmax2~=0);
ac = as(Indexc);
bc = bs(Indexc);
thetac= thetas(Indexc);
[xc,yc] = ind2sub(size(f),Indexc);


figure(); imshow(RGB); hold on;
plot(yc,xc,'g.');
if (~isempty(xc))
h=ellipse(ac,bc, pi/2 - thetac,yc,xc,'y');
end


phi = InitialLevelSetEllipse(xc,yc,ac*0.7,bc*0.7,thetac,image,'dark'); % last parameter depends on total image intensity.
phi1 = InitialLevelSetEllipse(xc,yc,ac*0.3,bc*0.3,thetac,image,'dark'); % last parameter depends on total image intensity.
bw = zeros(size(phi));
bw(find(phi==-1)) = 0;
bw(find(phi==1)) = 1;

bw1 = zeros(size(phi1));
bw1(find(phi1==-1)) = 0;
bw1(find(phi1==1)) = 1;



imgDist=imimposemin(image,bw);
DL = watershed(imgDist);
bgm = DL == 0;
se = strel('disk',1);
bgm = imdilate(bgm,se);
watershed_markers = bgm | (phi1>0) ;
out_red   = RGB(:,:,1);
out_green = RGB(:,:,2);
out_blue  = RGB(:,:,3);
out_red(watershed_markers)   = 0;   % Class of out_red is uint8, so use 255 instead of 1.
out_green(watershed_markers) = 255;
out_blue(watershed_markers)  = 0;
Markers = cat(3, out_red, out_green, out_blue);
figure();imshow(Markers);
imwrite(Markers,'markers.png');
title('Inner and Outer markers');


hy = fspecial('sobel');
hx = hy';
Iy = imfilter(double(image), hy, 'replicate');
Ix = imfilter(double(image), hx, 'replicate');
gradmag = sqrt(Ix.^2 + Iy.^2);
imgDist2 = imimposemin(gradmag, bgm | bw1);
L = watershed(imgDist2);
bgm2 = L==0;

%  se = strel('disk',1);
% bgm2 = imdilate(bgm2,se);
bww = bwdist(bgm2);
figure(); imshow(bww,[])
bgm2(find(bww<=0.5))=1;
Lrgb = label2rgb(L, 'jet', 'w', 'shuffle');
figure ;imshow(Lrgb);
imwrite(Lrgb,'watershed_color.png');
out_red   = RGB(:,:,1);
out_green = RGB(:,:,2);
out_blue  = RGB(:,:,3);
out_red(bgm2)   = 0;   % Class of out_red is uint8, so use 255 instead of 1.
out_green(bgm2) = 255;
out_blue(bgm2)  = 0;
out = cat(3, out_red, out_green, out_blue);
figure();imshow(out);
imwrite(uint8(out),'watershed_outlined.png');



[phi,accpt_index] = InitialLevelSetEllipse(xc,yc,ac,bc,thetac,image,'dark');
figure(); imshow(RGB); hold on;
plot(yc(accpt_index),xc(accpt_index),'g.');
h=ellipse(ac(accpt_index),bc(accpt_index), pi/2 - thetac(accpt_index),yc(accpt_index),xc(accpt_index),'y');








