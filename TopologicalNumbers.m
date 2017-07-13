%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------------ Calculation of Topological Numbers -------------------%
%------ point(px,py) is a simple point if and only if T4 = T8 = 1. -------%
% Check paper Han, Xiao, Chenyang Xu, and Jerry L. Prince. "A topology ---%
% preserving level set method for geometric deformable models." (2003)----%
% All the implementations are vectorized to speed up the program. --------%
%----------- Copyright (c) 2016, Shaghayegh Taheri Hosseinabadi ----------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function T = TopologicalNumbers(p_idx,X,H,n)

if n==4
Nn4k2 = GeodesicNeighbor(4,p_idx,X,H,2); % n=4 , k=2
Nn4k2_window = Window(p_idx,Nn4k2,H);
Nn4k2_CC = bwconncomp(Nn4k2_window,4); % Connected components 
T =  Nn4k2_CC.NumObjects; %T4


else if n==8
Nn8k1 = GeodesicNeighbor(8,p_idx,X,H,1); % n=8 , k=1
Nn8k1_window = Window(p_idx,Nn8k1,H);
Nn8k1_CC = bwconncomp(Nn8k1_window,8); % Connected components 
T = Nn8k1_CC.NumObjects; %T8
    end
end



function window = Window(pindex,Nn,H)
window = zeros(3,3);
if (any(ismember(Nn,pindex)))
    window(2,2)=1;
end
if (any(ismember(Nn,pindex-H-1)))
    window(1,1)=1;
end
if (any(ismember(Nn,pindex-H)))
    window(2,1)=1;
end
if (any(ismember(Nn,pindex-H+1)))
    window(3,1)=1;
end
if (any(ismember(Nn,pindex-1)))
    window(1,2)=1;
end
if (any(ismember(Nn,pindex+1)))
    window(3,2)=1;
end
if (any(ismember(Nn,pindex+H-1)))
    window(1,3)=1;
end
if (any(ismember(Nn,pindex+H)))
    window(2,3)=1;
end
if (any(ismember(Nn,pindex+H+1)))
    window(3,3)=1;
end


function Nn = GeodesicNeighbor(n,p_idx,X,H,k)
% foreground : is a black and white image with 1s as the foreground
% X : indexs of foreground pixels
% p_idx : index of the input point
% n : connectivity type ( 4 or 8 for 2D)
% Nn1 for k=1 , Nn2 for k=2 (k : order of geodesic neighborhood)

Nn_star = nNeighbors(p_idx,n,'star',H);
M = 8; % 2D 
NM_star = nNeighbors(p_idx,M,'star',H);
Nn1 = intersect(X,Nn_star);
Nn2 = [];

if k==1
    Nn = Nn1;
    return;
else if k==2

for i=1:size(Nn1)
     Nn_nbr = nNeighbors(Nn1(i),n,'None',H);
     S = intersect(NM_star,Nn_nbr);
     S = intersect(S,X);
     Nn2 = union(Nn2,S);
end

Nn = Nn2;
    end
end

% function index = Coordinate2Index(px,py,Height)
% index = (py-1)*Height + px;
% 
% function [px,py] = Index2Coordinate(index,Height)
% px = mod(index,Height);
% px = (px==0)*Height + px;
% py = floor(index/Height)+1;


function Nn = nNeighbors(p_idx,n,mode,H)
% This functions get the point coordinates and returns the arrays
% containing n-neighbors of the point. Nn = neighbors including the point,
% Nn_star = neighbors excluding the point.
% p_idx : index of the input point

if (strcmp(mode,'None'))
if (n==8)   
Nn = [p_idx-H-1;p_idx-H;p_idx-H+1;p_idx-1;p_idx;p_idx+1;p_idx+H-1;...
      p_idx+H;p_idx+H+1;];
      
else if (n==4)
Nn = [p_idx-H;p_idx-1;p_idx;p_idx+1;p_idx+H];
    end
end
end

if (strcmp(mode,'star'))
if ( n ==8)    
Nn = [p_idx-H-1;p_idx-H;p_idx-H+1;p_idx-1;p_idx+1;p_idx+H-1;...
      p_idx+H;p_idx+H+1;];
           
else if (n==4)
Nn = [p_idx-H;p_idx-1; p_idx+1;p_idx+H];
    end
end
end


