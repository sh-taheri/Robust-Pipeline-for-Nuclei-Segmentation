%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------- Initializing Level Set Based on Initial Input Ellipses --------%
%----------- Copyright (c) 2016, Shaghayegh Taheri Hosseinabadi. ---------%
%------------------------- Concordia University --------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [mask,varargout] = InitialLevelSetEllipse(px,py,a,b,teta,image,mode)

T = mean(mean(image));
% T is a threshold for removing ellipses which mostly contain background
% pixels
height = size(image,1);
width = size(image,2);
 mask = false(height,width);
%  for x = 1:height
%      for y=1:width
%          
%          for i = 1:length(px)
%          term1 = (x-px(i))*cos(teta(i)) + (y-py(i))*sin(teta(i));
%          term2 = (x-px(i))*sin(teta(i)) - (y-py(i))*cos(teta(i));
%          d = term1^2/a(i)^2 + term2^2/b(i)^2;
%          if d<=1
%              mask(x,y) = true;
%          end
%          end
%      end
%  end

indexout = zeros(size(px));
[x,y] = meshgrid(1:height,1:width);
if (strcmp(mode, 'bright'))
 for i = 1:length(px)
 term1 = (x-px(i))*cos(teta(i)) + (y-py(i))*sin(teta(i));
 term2 = (x-px(i))*sin(teta(i)) - (y-py(i))*cos(teta(i));
 d = term1.^2/a(i)^2 + term2.^2/b(i)^2;
 index = find(d<=1);
 index2 = sub2ind(size(mask),x(index),y(index));
 pixels = image(index2);
 if (median(pixels(:))>T)
 mask(index2)= true;
 indexout(i) = 1;
 end
 end
else
  for i = 1:length(px)
 term1 = (x-px(i))*cos(teta(i)) + (y-py(i))*sin(teta(i));
 term2 = (x-px(i))*sin(teta(i)) - (y-py(i))*cos(teta(i));
 d = term1.^2/a(i)^2 + term2.^2/b(i)^2;
 index = find(d<=1);
 index2 = sub2ind(size(mask),x(index),y(index));
 pixels = image(index2);
 if (median(pixels(:))<T)
 mask(index2)= true;
 indexout(i) = 1;
 end
  end
 
end

 mask = double(mask)*2 -1;
 accpt_index = find(indexout==1);
 varargout{1} = accpt_index;

 
 