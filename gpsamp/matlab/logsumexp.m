function y = logsumexp(x, dim),

if nargin < 2,
  dim = 1;
end

m = max(x, [], dim);

y = m + log(sum(exp(bsxfun(@minus, x, m)), dim));
