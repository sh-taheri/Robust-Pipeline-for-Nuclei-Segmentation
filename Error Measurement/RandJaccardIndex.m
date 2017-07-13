%-------------------------------------------------------------------------%
%---------------- Calculating Rand Index and Jaccard Index ---------------%
%-------------------------------------------------------------------------%

function [RI,JI] = RandJaccardIndex(BW_S,BW_R)

%%%% Method 1 : Fast method %%%%
A1 = (BW_S==0)&(BW_R==0);
A2 = (BW_S==0)&(BW_R==1);
A3 = (BW_S==1)&(BW_R==1); 
A4 = (BW_S==1)&(BW_R==0);
n = size(BW_S,1)*size(BW_S,2);
if size(BW_S,1) ~= size(BW_R,1)|| size(BW_S,2) ~= size(BW_R,2)
    disp('Dimension of input images should be the same')
end

nA1 = sum(sum(A1));
nA2 = sum(sum(A2));
nA3 = sum(sum(A3));
nA4 = sum(sum(A4));

num = nchoosek(n,2)- nA1*nA2 - nA1*nA4...
- nA2*nA3 - nA3*nA4;

den1 = nchoosek(n,2);
den2 = nchoosek(n,2) - nchoosek(nA1,2) - nchoosek(nA2,2) - nchoosek(nA3,2) - ...
       nchoosek(nA4,2);
   
   
RI = num/den1;
JI = num/den2;


%%%% Method 2 : Direct method %%%%
% S = BW_S(:);
% R = BW_R(:);
% if length(S) ~= length(R)
%     disp('Dimension of input images should be the same')
% %      exit;
% end
% a = 0;
% b = 0;
% c = 0;
% d = 0;
% L = length(S);
% for i=1:L-1
%     for j=i+1:L
%        if   (R(i)==R(j) && S(i)==S(j))
%            a = a + 1;
%        end
%        if   (R(i)~=R(j) && S(i)==S(j))
%            b = b + 1;
%        end
%        if   (R(i)==R(j) && S(i)~=S(j))
%            c = c + 1; 
%        end
%        if   (R(i)~=R(j) && S(i)~=S(j)) 
%            d = d + 1; 
%        end
%     end
% end
% 
% RI = (a+d)/(a+b+c+d);
% JI = (a+d)/(b+c+d);






