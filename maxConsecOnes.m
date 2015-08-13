function Nmax = maxConsecOnes(td)
% This function determines if the number of consecutive ones is >= than
% N
%
% Input: 
%
% td - logical array
%
% Output: 
%
% Nmax - maximum number of consecutive ones 

Nmax = max( diff( [0 (find( ~ (td > 0) ) ) numel(td) + 1] ) - 1);

end