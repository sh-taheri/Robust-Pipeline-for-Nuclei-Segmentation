%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----- Reinitialize level set values to the signed distance function -----%
% All the implementations are vectorized to speed up the program. --------%
%----------- Copyright (c) 2016, Shaghayegh Taheri Hosseinabadi ----------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function out = Reinitialize_phi(phi)
out=double((phi>0).*(bwdist(phi<0)-0.5)-(phi<0).*(bwdist(phi>0)-0.5));
% out = out/max(max(abs(out)));
%dimension = size(phi,1)*size(phi,2);
% phi = phi * (dimension/(10^6));
%out = out * (dimension/(10^6));


