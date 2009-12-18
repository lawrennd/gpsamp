function [res, mlass] = drosSparseRegression(drosexp, inputs, output, priorindices),

spikepriors = 1 - [0.0751 0.1163 0.1729  0.0378  0.2387];

spikepriors = spikepriors(priorindices);

F = drosexp.fitmean(inputs, :)';
y = drosexp.fitmean(output, :)';
yvar = drosexp.fitvar(output, :)';

F = [F ./ repmat(sqrt(var(F)), [size(F, 1), 1]), ones(length(y), 1)];
yscale = sqrt(var(y));
y = y ./ yscale;
yvar = yvar ./ (yscale .^ 2);

A = F' * diag(1./yvar) * F;
Ainv = inv(A);
b = sum(F .* repmat(y ./ yvar, [1, length(inputs) + 1]))';
c = .5 * sum(y.^2 ./ yvar);

lls = zeros(1, 2^length(inputs));

for k=0:(2^length(inputs)-1),
  bitmask = bitget(uint32(k), length(inputs):-1:1);
  prior = double(bitmask) .* (1-spikepriors) + (double(1-bitmask) .* spikepriors);
  sigma_w = [double(bitmask).*1 + double((1-bitmask)).* 0.0001, 1];
  sigmainv = diag(1 ./ sigma_w) + A;
  lls(k+1) = -.5 * log(det(sigmainv)) - .5 * sum(log(sigma_w)) ...
      +.5 * b'*(sigmainv\b) + sum(log(prior));
end

res = exp(lls - max(lls));
res = res ./ sum(res);

[foo, k] = max(res);
mlass0 = bitget(uint32(k-1), length(inputs):-1:1);
mlass = zeros(1, 5);
mlass(priorindices) = mlass0;
