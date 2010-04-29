function y = logsumexp(x, dim),

if nargin < 2,
  dim = 1;
end

m = max(x, [], dim);
mm = repmat(m, size(x)./size(m));

y = m + log(sum(exp(x - mm), dim));
