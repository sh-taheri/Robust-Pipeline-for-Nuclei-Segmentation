%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ------------------------ Geodesic Neighborhood ------------------------ %
% -------- A Topology Preserving Level Set Method for Geometric --------- %
% ----------------------- Deformable Models (TLSM)  --------------------- %
% ------ point(px,py) is a simple point if and only if T4 = T8 = 1. ----- %
% ---------------  To speed up the process, we suggest to Creat --------- %
% --------- Look Up Table for Topology Preserving Constraint ------------ %
% ----------------------------------------------------------------------- %
%  ----------------
% | 256 |  32 |  4 |
% |----------------|
% | 128 |  16 |  2 |
% |----------------|
% | 64  |  8  |  1 |
%  ----------------
% Ex : The 400th element of vector Simple is 110010000 (= dec2bin(400)),
% which leads to the following window:
%  -------------
% | 1 |  0 |  0 |
% |-------------|
% | 1 |  1 |  0 |
% |-------------|
% | 0 |  0 |  0 |
%  -------------
%-------------------------------------------------------------------------%
% --------- Copyright (c) Shaghayegh Taheri Hosseinabadi, 2016 -----------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all;
close all;
clc;
Height =5;
Simple = ones(512,1);

for i=0:511
dec = dec2bin(i,9);
X0 = dec - '0';
X1 = reshape(X0,[3,3]);
X = padarray(X1,[1 1]);
Xbar = 1 - X;
T4X = TopologicalNumbers(13,find(X==1),Height,4);
T8Xbar = TopologicalNumbers(13,find(Xbar==1),Height,8);

 if(T4X ~=1 || T8Xbar~=1)
     Simple(i+1)=0;
 end
end

save 'TPLookUpTable.mat' Simple





