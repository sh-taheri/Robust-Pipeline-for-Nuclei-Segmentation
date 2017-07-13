

function gradient_chanvese = Gradient_ChanVese(phi,F) 

Chi1 = Heaviside(phi);
Chi2 = 1 - Heaviside(phi);

num1 = F.*Chi1;
num2 = F.*Chi2; 
C1 = sum(num1(:))/sum(Chi1(:));
C2= sum(num2(:))/sum(Chi2(:));

gradient_chanvese = -(F-C1).^2 + (F-C2).^2;