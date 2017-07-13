%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------------------------------------------------------%
%----------------------- dna-0.png Segmentation --------------------------%
%--------------------------- Result is perfect ---------------------------%
%-------------------------------------------------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;
clear all;
clc;

RGB = imread('DNA\\dna-14.png');
% RGB = imread('dna-0.png'); 
% image = double(rgb2gray((RGB)));
image = rgb2gray(RGB);
image = double(imadjust(image)); %adjusting the contrast.
figure(); imshow(image,[]);


nrows = size(image,1);
ncols = size(image,2);
fmax = zeros(size(image));
as = zeros(size(image));
bs = zeros(size(image));
thetas = zeros(size(image));
range = [50 60];


for a = range
for b=[0.6 0.8]*a
 for theta = 0:pi/8:pi 
  f = frst2dEllipse(image,a,b,theta, 0.8,0.5,'bright'); 
  fmax0 = localMaximum(f,[max(range) min(range)],true,0.001); 
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

mode = 'bright';
phi = InitialLevelSetEllipse(xc,yc,ac*0.9,bc*0.9,thetac,image,mode); % last parameter depends on total image intensity.
phi=double((phi>0).*(bwdist(phi<0))-(phi<0).*(bwdist(phi>0)));
F(:,:,1)= double(image);


dt = 2;
pdf = 'Gaussian';
handle_seg =figure();
B= phi>=0;
W_nb = 3;
load TPLookUpTable.mat Simple %Load the look up table for checking topology preserving constraint

disp('Segmentation');
for iteration = 1:150 %120

alpha = iteration*1.5;

    if(iteration>50 && iteration<140) 
        alpha = 50;
    end
if (iteration>140)
    alpha = 100;
end
  
phi = Reinitialize_phi(phi);
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

%{
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
%}

phi = phi_new; 
figure(handle_seg); 
imshow(RGB,[]);
hold on;
contour(phi,[0 0],'c','LineWidth',1);
title(['iteration ',num2str(iteration)]);
drawnow;
end


%%% Post Processing Operation (Morphology Opening)
phi=double((phi>0).*(bwdist(phi<0))-(phi<0).*(bwdist(phi>0)));
BW = phi>0.5;
BW = bwareaopen(BW,3);
se = strel('disk',3); 
BW = imopen(BW,se);
figure(); imshow(BW);

% BW = im2bw(phi);
% Segments = bwlabel(BW,4);
% RGBSegments = label2rgb(Segments,'jet',[1 1 1],'shuffle');
% figure(); 
% imshow(RGBSegments);

% Adding nonlinear diffusion and thresholding 
% image = imadjust(uint8(image)); %adjusting the contrast.
image = Nonlinear_Diffusion(image,1,1,300);
Threshold_BW = im2bw(uint8(image),0.2);
difference = Threshold_BW - BW;
difference = difference==1; 
[Segments_Threshold,n] = bwlabel(difference,4);
T=5;
A = min(range)^2*pi;
for i=1:n
[x,y] = find(Segments_Threshold==i);
out = CheckObjectRemoval(x,y,size(image,1),size(image,2),T,A);
if out==1
    difference(Segments_Threshold==i)=0;
end
end
se = strel('disk',3); 
difference = imerode(difference,se);
final_BW = (BW==1)|(difference==1);

%%% Final Results %%%
figure(); imshow(RGB,[]);
hold on;
contour(BW,[0 0],'c','LineWidth',1);
title('Proposed Method');

figure(); imshow(RGB,[]);
hold on;
contour(final_BW,[0 0],'c','LineWidth',1);
title('Proposed Method with Nonlinear diffusion and Thresholding');

imwrite(final_BW,'ProposedMethodWithThresholding.png');
imwrite(BW,'ProposedMethod.png');


Segments = bwlabel(final_BW,4);
RGBSegments = label2rgb(Segments,'jet',[1 1 1],'shuffle');
figure(); 
imshow(RGBSegments);
