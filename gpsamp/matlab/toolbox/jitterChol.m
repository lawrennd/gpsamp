function [L er] = jitterChol(K)
% function [L er] = jitterChol(K)
%
% Description:  Computing Choleski decomposition by adding jitter  
%              when the matrix is semipositive definite  
%

jitter = 1e-7;
m = size(K,1); 
[L er] = chol(K);
if er > 0 % add jitter
   K = K + jitter*mean(diag(K))*eye(m);
   [L er] = chol(K);
end