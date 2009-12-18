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
    results(k) = 1./mean(1 ./ testGenes{k}.LogL);
  end
  return;
end

results = zeros(size(testGenes{1}.Weights, 1), length(testGenes));

for k=1:length(testGenes),
  W = testGenes{k}.Weights;
  switch method,
   case 'zscore',
    results(:, k) = mean(W, 2) ./ std(W, 0, 2);
   case 'pplr2',
    results(:, k) = max(mean(W < -.02, 2), mean(W > .02, 2));
  end
end
