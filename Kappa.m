%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------- Discretization of divergence, Check paper : --------------%
%-------- "Nonlinear total variation based noise removal alg. 1992 -------%
%----------- Copyright (c) 2016, Shaghayegh Taheri Hosseinabadi ----------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function out = Kappa(u) % Kappa(u) = div(grad(u)/|grad(u)|)

out = zeros(size(u));
out1 = zeros(size(u)); 
out2 = zeros(size(u));
u = padarray(u,[1,1],'replicate');
for i=2:size(u,1)-1
    for j=2:size(u,2)-1        
num_x = u(i+1,j)-u(i,j);
den_x = sqrt((u(i+1,j)-u(i,j))^2 + (m(u(i,j+1)-u(i,j) , u(i,j)-u(i,j-1)))^2);
num_y = u(i,j+1)-u(i,j);
den_y = sqrt((u(i,j+1)-u(i,j))^2 + (m(u(i+1,j)-u(i,j) , u(i,j)-u(i-1,j)))^2);
out1(i-1,j-1) = num_x/(den_x+eps);
out2(i-1,j-1) = num_y/(den_y+eps);
    end
end

out1_padded = padarray(out1,[1,1],'replicate');
out2_padded = padarray(out2,[1,1],'replicate');

for i=2:size(u,1)-1
    for j=2:size(u,2)-1  
out(i-1,j-1) = out1_padded(i,j)-out1_padded(i-1,j) + ...
out2_padded(i,j)-out2_padded(i,j-1);   
    end
end


function minmod = m(a,b)
minmod = ((sign(a) + sign(b))/2 )*min(abs(a),abs(b));





%------------------------------ Older codes ------------------------------%
% central difference
% fy = P(3:end,2:n+1)-P(1:m,2:n+1);
% fx = P(2:m+1,3:end)-P(2:m+1,1:n);
% fyy = P(3:end,2:n+1)+P(1:m,2:n+1)-2*I;
% fxx = P(2:m+1,3:end)+P(2:m+1,1:n)-2*I;
% fxy = 0.25.*(P(3:end,3:end)-P(1:m,3:end)+P(3:end,1:n)-P(1:m,1:n));
% G = (fx.^2+fy.^2).^(0.5);
% K = (fxx.*fy.^2-2*fxy.*fx.*fy+fyy.*fx.^2)./((fx.^2+fy.^2+eps).^(1.5));
% KG = K.*G;
% KG(1,:) = eps;
% KG(end,:) = eps;
% KG(:,1) = eps;
% KG(:,end) = eps;
% KG = KG./max(max(abs(KG)));


% eps = 0.00001;
% u=image;
% [ux,uy]=gradient(u);
% G=sqrt(ux.^2+uy.^2);
% div =divergence(ux./(G+eps),uy./(G+eps));
