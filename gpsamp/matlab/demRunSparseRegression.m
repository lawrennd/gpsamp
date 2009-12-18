tfs = {drosTF.names, {'mef2', 'twi', 'bin'}, {'mef2', 'twi', 'tin'}};
priorinds = {1:5, [5, 3, 2], [5, 3, 1]};

for l=1:length(tfs),
  fbgns = {};
  for k=1:length(tfs{l}),
    fbgns{k} = drosTF.fbgns.(tfs{l}{k});
  end
  inputs = drosFindGeneinds(drosexp, fbgns);
  
  regression_res{l} = zeros(length(testset.indices), 2^length(inputs));
  regression_ass{l} = zeros(length(testset.indices), 5);

  for k=1:length(testset.indices),
    [regression_res{l}(k, :), regression_ass{l}(k, :)] = drosSparseRegression(drosexp, inputs, testset.indices(k), priorinds{l});
  end
end
