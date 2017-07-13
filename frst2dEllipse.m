%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------- Generalized Fast Radial Symmetry.-------------------%
%------- Check paper Ni, Jie, Maneesh K. Singh, and Claus Bahlmann.-------% 
% Fast radial symmetry detection under affine transformations. CVPR 2012  %
% --- All the implementations are vectorized to speed up the program. ----%
%------------------- Symmetry under Affine Transformation)----------------%
%----------- Copyright (c) 2016, Shaghayegh Taheri Hosseinabadi ----------%
% (Based on "MATLAB implementation of the Fast Radial Symmetry Transform -% 
% in 2D , Sandro, Mar 2014)-----------------------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [ filtered] = frst2dEllipse( image,a,b,teta, alpha, stdFactor, mode )
    bright = false;
    dark = false;
    if (strcmp(mode, 'bright'))
        bright = true;
    elseif (strcmp(mode, 'dark'))
        dark = true;
    elseif (strcmp(mode, 'both'))
        bright = true;
        dark = true;
    else
        error('invalid mode');
    end
    
    original = double(image);
    [gy, gx] = gradient(original);
    maxRadius = ceil(max(a(:)));
    offset = [maxRadius maxRadius];  
    filtered = zeros(size(original) + 2*offset); 
    S = zeros([numel(a)*numel(b)*numel(teta), size(filtered, 1), size(filtered, 2)]); 
    radiusIndex = 1;

    for ai = a
        for bi = b
            for tetai = teta 
                
        n = ai; 
        
%         O_n = zeros(size(filtered));
%         M_n = zeros(size(filtered));
        
         %%%% Affine Transformation %%%%
        R =[cos(tetai) -sin(tetai); sin(tetai) cos(tetai)];
        Sc = [ai 0; 0 bi];     
        M = [0 1;-1 0];
        G = R*Sc;
        V0 = G*M/G/M;
        VX = V0(1,1)*gx + V0(1,2)*gy;
        VY = V0(2,1)*gx + V0(2,2)*gy;
        

        x0 = 1:size(original,1);
        x = repmat(x0',size(original,2),1);
        y0 = 1:size(original,2);
        y = repmat(y0,size(original,1),1);
        y=y(:);
        gnorm = sqrt((VX(:)).^2 + (VY(:)).^2);
        nonzero = find(gnorm>0);
        g1 = VX(nonzero)./gnorm(nonzero);
        g2 = VY(nonzero)./gnorm(nonzero);
        g = [g1 g2];
        gp = round((g) * n);        
        p = [x(nonzero) y(nonzero)];
        offset = repmat(offset,[size(nonzero),1]);
        

        
         if (bright)
             ppve = p + gp;
             ppve = ppve + offset;
             index = sub2ind(size(filtered),ppve(:,1),ppve(:,2));
             Acc_O_n = accumarray(index,1,[size(filtered,1)*size(filtered,2) 1]);
             Acc_M_n = accumarray(index,gnorm(nonzero),[size(filtered,1)*size(filtered,2) 1]);
         end
          
         
          if (dark)
             pnve = p - gp;
             pnve = pnve + offset;
             index = sub2ind(size(filtered),pnve(:,1),pnve(:,2));
             Acc_O_n = accumarray(index,-1,[size(filtered,1)*size(filtered,2) 1]); 
             Acc_M_n = accumarray(index,-gnorm(nonzero),[size(filtered,1)*size(filtered,2) 1]);
          end
         

          O_n = reshape(Acc_O_n,size(filtered,1),size(filtered,2));
          M_n = reshape(Acc_M_n,size(filtered,1),size(filtered,2));
          
          
          
          
%{      
        for j=1:size(original, 2)
            for i = 1:size(original, 1)
                                            
                    p = [i j];
%                   V = G*M*inv(G)*inv(M)*[gx(i,j);gy(i,j)];
%                     V = G*M/G/M*[gx(i,j);gy(i,j)];
%                     vx = V(1);
%                     vy = V(2);
%                     g = [vx vy];                      
                     g = [VX(i,j) VY(i,j)];
                    gnorm = sqrt( g * g' ) ;
                 if (gnorm > 0)
                      gp = round((g./gnorm) * n);

                   if (bright)
                        ppve = p + gp;
                        ppve = ppve + offset;
                        ppve1 = ppve(1);
                        ppve2 = ppve(2);

                        O_n(ppve1, ppve2) = O_n(ppve1, ppve2) + 1;
                        M_n(ppve1, ppve2) = M_n(ppve1, ppve2) + gnorm;
                    end
                    
                    if (dark)
                        pnve = p - gp;
                        pnve = pnve + offset;
                        pnve1 = pnve(1);
                        pnve2 = pnve(2);

                        O_n(pnve1, pnve2) = O_n(pnve1, pnve2) - 1;
                        M_n(pnve1, pnve2) = M_n(pnve1, pnve2) - gnorm;
                    end
                end
                
            end
        end
        
        %}
      

        O_n = abs(O_n);
        O_n = O_n ./ max(O_n(:));
        
        M_n = abs(M_n);
        M_n = M_n ./ max(M_n(:));
        
        
        S_n = (O_n.^alpha) .* M_n; 
%       S_n = (O_n.^alpha) .* sign(O_n);
         
        
        gaussian = fspecial('gaussian', [ceil(n/2) ceil(n/2)], n*stdFactor);
        S(radiusIndex, :, :) = imfilter(S_n, gaussian);
        
        radiusIndex = radiusIndex + 1;
            end
        end
    end
  
    filtered = squeeze(sum(S, 1));
    filtered = filtered(offset(1)+1:end-offset(2), offset(1)+1:end-offset(2));

    
   
