%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----- Creating the Look Up table for topology preserving constraint -----%
% All the implementations are vectorized to speed up the program. --------%
%----------- Copyright (c) 2016, Shaghayegh Taheri Hosseinabadi ----------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function n = LT(p_index,X,H) % Look Up table value

n = 0;
Xvector = X(:);
N = size(X,1)*size(X,2);

index = p_index - H - 1;
if index > 0
    n = n + 256* Xvector(index);
end

index = p_index - H ;
if index > 0
    n = n + 128* Xvector(index);
end

index = p_index - H + 1 ;
if index > 0
    n = n + 64* Xvector(index);
end

index = p_index - 1 ;
if index > 0
    n = n + 32 * Xvector(index);
end


n = n + 16 * Xvector(p_index);


index = p_index + 1 ;
if index < N
    n = n + 8 * Xvector(index);
end

index = p_index + H - 1 ;
if index <  N
    n = n + 4 * Xvector(index);
end

index = p_index + H ;
if index < N
    n = n + 2 * Xvector(index);
end

index = p_index + H + 1 ;
if index < N
    n = n + Xvector(index);
end

