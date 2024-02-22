function [s, mask] = vectorize_mask(S, mask)
if nargin == 1
   mask = ~isnan(S(:,:,:,1));
end
if ismatrix(S)
n = size(S, 1);
[x, y, z] = size(mask);
s = NaN([x, y, z, n], 'like', S);
for i = 1:n
    tmp = NaN(x, y, z, 'like', S);
    tmp(mask(:)) = S(i, :);
    s(:,:,:,i) = tmp;
end
else
nl = length( mask(mask(:)) );  % new RickCL  
s  = zeros(size(S, 4), nl , 'like', S); % new RickCL    
for i = 1:size(S, 4)
   Si = S(:,:,:,i);
   s(i, :) = Si(mask(:));
end
end
end
   
