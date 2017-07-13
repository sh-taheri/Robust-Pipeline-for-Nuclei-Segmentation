
function mask = InitialLevelSet(px,py,r,width,height)

window = false(2*r-1,2*r-1);
mask = false(height,width);
for ix = 1:2*r-1
    for iy = 1:2*r-1
        d = Distance(ix,iy,r,r);
        if (d<r)
            window(ix,iy)=true;
        end   
    end
end

for i = 1:size(px)
    xindex = (px(i) - r + 1) : (px(i) + r - 1);
    yindex = (py(i) - r + 1) : (py(i) + r - 1);
if ((xindex(1)>0) && (xindex(2*r-1)<height) && (yindex(1)>0) && ...
    (yindex(2*r-1)<width))
mask(xindex,yindex) = mask(xindex,yindex)|window;
end
end
mask = double(mask)*2 -1;


                
function d = Distance(px,py,cx,cy)
d = sqrt((px-cx).^2 + (py-cy).^2);
        