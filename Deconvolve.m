%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------- Deconvolve Applies Color Deconvolution Algorithm to image -------%
%------ Check paper Ruifrok, Arnout C., and Dennis A. Johnston (2001.-----%
%---- Quantification of histochemical staining by color deconvolution ----%
%----------- Copyright (c) 2016, Shaghayegh Taheri Hosseinabadi. ---------%
%------------------------- Concordia University --------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Hematoxylin,Eosin,BG,RGB_Hematoxylin,RGB_Eosin,...
          RGB_BG] = Deconvolve(image)

 R = double(image(:,:,1));
 G = double(image(:,:,2));
 B = double(image(:,:,3));
 eps = 0.00001;

RGB_Hematoxylin = uint8(zeros(size(image,1),size(image,2),3));
RGB_Eosin = uint8(zeros(size(image,1),size(image,2),3));
RGB_BG = uint8(zeros(size(image,1),size(image,2),3));

 
 %%%---------------------- Hematoxylin, Eosin, DAB ---------------------%%%
  He = [0.65 0.7 0.25]'; 
  Eo = [0.07 0.99 0.11]';
  BG = [0.25 0.57 0.78]';

%%%---------------------------  H&E version 1 ------------------------- %%%
%  He = [0.49015734 0.76897085 0.41040173]';
%  Eo = [0.04615336 0.8420684 0.5373925]';
%  BG = [eps eps eps]'; 

%%%--------------------------- H&E version 2 -------------------------- %%%
%  He = [0.644211 0.716556 0.266844]';
%  Eo = [0.04615336 0.8420684 0.5373925]';
%  BG = [eps eps eps]';

%%%--------------------------- H&E version 3 -------------------------- %%%
%  He = [0.65 0.71 0.31040173]';
%  Eo = [0.04615336 0.7042 0.5373925]';
%  BG = [0.8648 0  0.7842]';
 
 %%%--------------------------- H&E version 4 -------------------------- %%%
%  He = [0.65 0.6 0.31040173]';
%  Eo = [0.07 0.9 0.5373925]';
%  BG = [0.8632 0  0.7842]';

M = [He/norm(He) Eo/norm(Eo) BG/norm(BG)];
index = 1:(size(image,1)*size(image,2));
sai = [R(index);G(index);B(index)];
sai = (sai+1)/256;
OD = -log(sai+eps);
saihat = M\OD; 
Hematoxylin = reshape(saihat(1,:),size(image,1),size(image,2));
Eosin =  reshape(saihat(2,:),size(image,1),size(image,2));
BG =  reshape(saihat(3,:),size(image,1),size(image,2));

% Enhance contrast
Hematoxylin = (Hematoxylin - 0.3) .* (Hematoxylin > 0.3);
Eosin = (Eosin - 0.3) .* (Eosin > 0.3);
BG = (BG - 0.3) .* (BG > 0.3);

RGB_Hematoxylin(:,:,1) = 256*exp(-M(1,1)*Hematoxylin)-1;
RGB_Hematoxylin(:,:,2) = 256*exp(-M(2,1)*Hematoxylin)-1;
RGB_Hematoxylin(:,:,3) = 256*exp(-M(3,1)*Hematoxylin)-1;

RGB_Eosin(:,:,1) = 256*exp(-M(1,2)*Eosin)-1;
RGB_Eosin(:,:,2) = 256*exp(-M(2,2)*Eosin)-1;
RGB_Eosin(:,:,3) = 256*exp(-M(3,2)*Eosin)-1;

RGB_BG(:,:,1) = 256*exp(-M(1,3)*BG)-1;
RGB_BG(:,:,2) = 256*exp(-M(2,3)*BG)-1;
RGB_BG(:,:,3) = 256*exp(-M(3,3)*BG)-1;


Hematoxylin = uint8( 255* mat2gray(Hematoxylin));
Eosin = uint8(255*mat2gray(Eosin));
BG = uint8(255*mat2gray(BG));




