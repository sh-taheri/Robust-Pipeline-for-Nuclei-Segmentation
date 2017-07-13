
function out = CheckObjectRemoval(x,y,height,width,T,A) % x : x,y set of all points in object
                                         % threshold value (distance from
                                         % border of image)
                                         % A : (area) size of object.
                                         % height, width: dimension of
                                         % image.
out =1;
xmin = min(x);
xmax = max(x);
ymin = min(y); 
ymax = max(y);

if ((xmin <T) | (xmax > (height-T)) | (ymin<T) | (ymax > (width -T)))  %Checking if object is close to border 
    out = 0;
else if (length(x)>A) %Check if object is large enough.
        out=0;
    end
end

