function results = gpmtfSummariseResults(testGenes, method, Genes, GeneVars, TFs, models, tGrid)

if nargin < 2,
  method='zscore';
end

if strcmp(method, 'thermodyn')
  Dtn = tGrid(2:end)-tGrid(1:end-1);
  for k=1:size(testGenes,1) 
    fprintf('Running gene %d/%d...\n', k, size(testGenes, 1));
    for c=1:size(testGenes,2)
       tmpResults = zeros( size(tGrid) );   
       for f=1:size(testGenes,3)
          tmpResults(f) = mean(testGenes{k,c,f}.LogL);
       end
       results(k,c) = 0.5*sum(Dtn.*(tmpResults(1:end-1) + tmpResults(2:end)));
    end
  end
  return;
elseif strcmp(method, 'margLogLik2')
  % Chib's approximatino to the marginal likelihood  
  results = zeros( size(testGenes) );
  % compute all marginal likelihoods for all trained models
  for k=1:size(testGenes,1) 
    fprintf('Running gene %d/%d...\n', k, size(testGenes, 1));
  for c=1:size(testGenes,2) 
  for f=1:size(testGenes,3)    
    %
    if size(testGenes,3) > 1    
        modelTest = models{c,f};
    else
        modelTest = models{c};
    end
    lnpriorKin = ['ln',modelTest.prior.kinetics.type,'pdf'];
    TrspaceKin = modelTest.prior.kinetics.priorSpace;
    lnpriorNoiseWHITE = ['ln',modelTest.prior.sigma2.type,'pdf'];
    TrspaceNoiseWHITE = modelTest.prior.sigma2.priorSpace; 
    
    % approximate the entropy using Gaussian density estimation 
    if isfield(testGenes{k,c,f}, 'W')
    %  
       lnpriorW0 = ['ln',modelTest.prior.W0.type,'pdf'];
       TrspaceW0 = modelTest.prior.W0.priorSpace;  
       lnpriorW = ['ln',modelTest.prior.W.type,'pdf'];
       TrspaceW = modelTest.prior.W.priorSpace;
       
       kin = feval(TrspaceKin, testGenes{k,c,f}.kinetics);
       W = feval(TrspaceW, testGenes{k,c,f}.W);
       W0 = feval(TrspaceW0, testGenes{k,c,f}.W0);
       sigma2 = feval(TrspaceNoiseWHITE, testGenes{k,c,f}.sigma2);
       
       X = [kin; W; W0; sigma2]';
       % Store also a representative paramrters for the Chib;s marginal
       % likelihood (currently is the max in term of fititng the data)
       [mm mind] = max(testGenes{k,c,f}.LogL);
       parmu = [kin(:,mind); W(:,mind); W0(1,mind); sigma2(1,mind)]';
       %X = [testGenes{k,c}.kinetics; testGenes{k,c}.W; testGenes{k,c}.W0]';
    else
       
       kin = feval(TrspaceKin, testGenes{k,c,f}.kinetics);
       sigma2 = feval(TrspaceNoiseWHITE, testGenes{k,c,f}.sigma2);
        
       X = [kin; sigma2]'; 
       % Store also a representative paramrters for the Chib;s marginal
       % likelihood (currently is the max in term of fititng the data)
       [mm mind] = max(testGenes{k,c,f}.LogL);
       parmu = [kin(:,mind); sigma2(1,mind)]';
       %X = [testGenes{k,c}.kinetics]';
    end
    
    mu = mean(X,1);
    parmu = mu;
    jit = 1e-8;
    Sigma = cov(X) + jit*eye(size(X,2));
    [N D] = size(X);                              
    logPosteriorOfSelectedPoint = - (0.5*D)*log(2*pi)  - 0.5*log(det(Sigma)) ...
                                  - 0.5*((parmu-mu)*(Sigma\(parmu-mu)'));
    
    % Chib's approximation to the marginal likelihood 
    % use values with likelhiood under the model and evaluate the Chib's 
    % method there 
    modelTest.Likelihood.Genes = Genes(k,:,:);
    if isfield(modelTest.Likelihood.noiseModel, 'pumaSigma2') == 1
        modelTest.Likelihood.noiseModel.pumaSigma2(1,:,:) = GeneVars(k,:,:);
    end
    if modelTest.Likelihood.numTFs == 0
    %    
       if strcmp(TrspaceKin, 'log')
           modelTest.Likelihood.kinetics(1) = exp(parmu(1));
           modelTest.Likelihood.kinetics(2) = exp(parmu(2));
           modelTest.Likelihood.kinetics(3) = exp(parmu(3));
       elseif strcmp(TrspaceKin, 'lin')
           modelTest.Likelihood.kinetics(1) = parmu(1);
           modelTest.Likelihood.kinetics(2) = parmu(2);
           modelTest.Likelihood.kinetics(3) = parmu(3);
       end
       if strcmp(TrspaceNoiseWHITE, 'log')
           modelTest.Likelihood.noiseModel.sigma2 = exp(parmu(4));
       elseif strcmp(TrspaceNoiseWHITE, 'lin')
           modelTest.Likelihood.noiseModel.sigma2 = parmu(4); 
       end
       if ~isfield(modelTest.Likelihood, 'crValMask') 
           modelTest.Likelihood.crValMask = 1:modelTest.Likelihood.numTimes;
       end
       % compute the maeginal likelhiod (GPs fucntino marginalized out)
       PredGenes = modelTest.Likelihood.kinetics(1)/modelTest.Likelihood.kinetics(2) +...
                  (modelTest.Likelihood.kinetics(3) - ...
                  modelTest.Likelihood.kinetics(1)/modelTest.Likelihood.kinetics(2))*exp(-modelTest.Likelihood.TimesG*modelTest.Likelihood.kinetics(2));
       L = 0;
       %     
       for r=1:modelTest.Likelihood.numReplicas
            sigmas = zeros(1, modelTest.Likelihood.numTimes);
            if modelTest.Likelihood.noiseModel.active(1) == 1
                sigmas = modelTest.Likelihood.noiseModel.pumaSigma2(1, : , r);
            end
            if modelTest.Likelihood.noiseModel.active(2) == 1
                sigmas = sigmas + repmat(modelTest.Likelihood.noiseModel.sigma2, 1, modelTest.Likelihood.numTimes ); 
            end
            L = L - 0.5*sum(log(2*pi*sigmas(modelTest.Likelihood.crValMask)),2)....
                  - 0.5*sum(((modelTest.Likelihood.Genes(1, modelTest.Likelihood.crValMask, r)...
                  - PredGenes(modelTest.Likelihood.crValMask)).^2)./sigmas(modelTest.Likelihood.crValMask),2);
       end
       LogLik = ones(1,size(TFs,2))*L;
    else
      if strcmp(TrspaceKin, 'log')
           modelTest.Likelihood.kinetics(1,:) = exp(parmu(1:4));
      elseif strcmp(TrspaceKin, 'lin')
           modelTest.Likelihood.kinetics(1,:) = parmu(1:4);
      end
      if strcmp(TrspaceW, 'log')
           modelTest.Likelihood.W(1,:) = exp(parmu(5:4+modelTest.Likelihood.numTFs));
      elseif strcmp(TrspaceW, 'lin')
           modelTest.Likelihood.W(1,:) = parmu(5:4+modelTest.Likelihood.numTFs);
      end
      if strcmp(TrspaceW0, 'log')
           modelTest.Likelihood.W0 = exp(parmu(5+modelTest.Likelihood.numTFs));
      elseif strcmp(TrspaceW0, 'lin')
           modelTest.Likelihood.W0 = parmu(5+modelTest.Likelihood.numTFs);
      end
      if strcmp(TrspaceNoiseWHITE, 'log')
           modelTest.Likelihood.noiseModel.sigma2 = exp(parmu(end));
      elseif strcmp(TrspaceNoiseWHITE, 'lin')
           modelTest.Likelihood.noiseModel.sigma2 = parmu(end);
      end
      TFset = find(modelTest.Likelihood.TFcomb==1);
      LogLik = zeros(1,size(TFs,2));
      %     
      for t=1:size(TFs,2)
         modelTest.Likelihood.TF = TFs{t}(TFset,:,:);
         for r=1:modelTest.Likelihood.numReplicas          
            LogLik(t) = LogLik(t) + gpmtfLogLikelihoodGene(modelTest.Likelihood, zeros(size(modelTest.Likelihood.TF(:,:,1))), r, 1);
         end   
      end
    %  
    end
    
    % evaluate the priors on the mu vector
    Likkin = feval(TrspaceKin, modelTest.Likelihood.kinetics);
    LogPriorKin = feval(lnpriorKin, Likkin, modelTest.prior.kinetics);
    
    % white GP noise model is used
    temp = feval(TrspaceNoiseWHITE, modelTest.Likelihood.noiseModel.sigma2);
    LogPriorNoiseWHITE = feval(lnpriorNoiseWHITE, temp, modelTest.prior.sigma2);
    LogPrior = sum(LogPriorKin) + sum(LogPriorNoiseWHITE);
    
    if modelTest.Likelihood.numTFs > 0
    %    
       % evaluation of the prior for the interaction bias
       LikW0 = feval(TrspaceW0, modelTest.Likelihood.W0);
       LogPriorW0 = feval(lnpriorW0, LikW0, modelTest.prior.W0);
       % evaluation of the prior for the interaction weights
       LikW = feval(TrspaceW, modelTest.Likelihood.W);
       LogPriorW = feval(lnpriorW, LikW, modelTest.prior.W);
       LogPrior = LogPrior + sum(LogPriorW0) + sum(LogPriorW);
    %
    end
    LogMargL(k,c,f) = logsumexp(LogLik, 2) - log(size(TFs,2))  +  LogPrior - logPosteriorOfSelectedPoint;
    %
  end
  end
  end
  results = real(LogMargL);
  return;
  %
elseif strcmp(method, 'margLogLik1')
  % 
  results = zeros( size(testGenes) );
  % compute all marginal likelihoods for all trained models
  for k=1:size(testGenes,1) 
  for c=1:size(testGenes,2) 
  for f=1:size(testGenes,3)    
    %
    if size(testGenes,3) > 1    
        modelTest = models{c,f};
    else
        modelTest = models{c};
    end
    lnpriorKin = ['ln',modelTest.prior.kinetics.type,'pdf'];
    TrspaceKin = modelTest.prior.kinetics.priorSpace;
    lnpriorNoiseWHITE = ['ln',modelTest.prior.sigma2.type,'pdf'];
    TrspaceNoiseWHITE = modelTest.prior.sigma2.priorSpace; 
    
    % approximate the entropy using Gaussian density estimation 
    if isfield(testGenes{k,c,f}, 'W')  
       lnpriorW0 = ['ln',modelTest.prior.W0.type,'pdf'];
       TrspaceW0 = modelTest.prior.W0.priorSpace;  
       lnpriorW = ['ln',modelTest.prior.W.type,'pdf'];
       TrspaceW = modelTest.prior.W.priorSpace;
       kin = feval(TrspaceKin, testGenes{k,c,f}.kinetics);
       W = feval(TrspaceW, testGenes{k,c,f}.W);
       W0 = feval(TrspaceW0, testGenes{k,c,f}.W0);
       sigma2 = feval(TrspaceNoiseWHITE, testGenes{k,c,f}.sigma2);
       X = [kin; W; W0; sigma2]';
       %X = [testGenes{k,c}.kinetics; testGenes{k,c}.W; testGenes{k,c}.W0; testGenes{k,c,f}.sigma2]';
    else
       kin = feval(TrspaceKin, testGenes{k,c,f}.kinetics);
       sigma2 = feval(TrspaceNoiseWHITE, testGenes{k,c,f}.sigma2);
       X = [kin; sigma2]'; 
       %X = [testGenes{k,c,f}.kinetics; testGenes{k,c,f}.sigma2]'; 
       %X = [testGenes{k,c}.kinetics]';
    end
    mu = mean(X,1);
    Sigma = cov(X);
    [N D] = size(X);
    jit = 1e-10;
    entropy = (0.5*D)*log(2*pi) + (0.5*D) + 0.5*log(det(Sigma + jit*eye(size(Sigma,1))));
    LogMargL(k,c,f) = mean(testGenes{k,c,f}.LogL) + mean(testGenes{k,c,f}.LogPrior) + entropy;
    %
  end
  end
  end
  results = real(LogMargL);
  return;
  %
elseif strcmp(method, 'loglike'),
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

