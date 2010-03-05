function results = gpmtfSummariseResults(testGenes, method),

if nargin < 2,
  method='zscore';
end

if strcmp(method, 'loglike'),
  results = zeros(1, length(testGenes));
  for k=1:length(testGenes),
    results(k) = mean(testGenes{k}.LogL);
  end
  return;
elseif strcmp(method, 'harmmeanlike'),
  results = zeros(1, length(testGenes));
  for k=1:length(testGenes),
    results(k) = log(1./mean(1 ./ exp(testGenes{k}.LogL)));
  end
  return;
elseif strcmp(method, 'meansens'),
  results = zeros(1, length(testGenes));
  for k=1:length(testGenes),
    results(k) = mean(testGenes{k}.kinetics(2, :));
  end
  return;
elseif strcmp(method, 'meansigma'),
  if isfield(testGenes{1}, 'sigmas'),
    results = zeros(1, length(testGenes));
    for k=1:length(testGenes),
      results(k) = mean(testGenes{k}.sigmas);
    end
  else
    results = [];
  end
  return;
elseif strcmp(method, 'sigma2f'),
  if isfield(testGenes{1}, 'sigma2f'),
    results = zeros(1, length(testGenes));
    for k=1:length(testGenes),
      results(k) = mean(testGenes{k}.sigma2f);
    end
  else
    results = [];
  end
  return;
elseif strcmp(method, 'lengthScale'),
  if isfield(testGenes{1}, 'lengthScale'),
    results = zeros(1, length(testGenes));
    for k=1:length(testGenes),
      results(k) = mean(testGenes{k}.lengthScale);
    end
  else
    results = [];
  end
  return;
end

if isfield(testGenes{1}, 'Weights'),
  results = zeros(size(testGenes{1}.Weights, 1), length(testGenes));
else
  results = zeros(size(testGenes{1}.W, 1), length(testGenes));
end

for k=1:length(testGenes),
  if isfield(testGenes{k}, 'Weights'),
    W = testGenes{k}.Weights;
  else
    W = testGenes{k}.W;
  end
  d = size(W, 1);
  switch method,
   case 'meansensweight',
    results(:, k) = mean(W .* repmat(testGenes{k}.kinetics(2, :), [d, 1]), 2);
   case 'meanweight',
    results(:, k) = mean(W, 2);
   case 'zscore',
    results(:, k) = mean(W, 2) ./ std(W, 0, 2);
   case 'pplr2',
    results(:, k) = max(mean(W < -.02, 2), mean(W > .02, 2));
  end
end
