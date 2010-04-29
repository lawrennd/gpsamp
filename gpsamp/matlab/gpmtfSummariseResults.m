function results = gpmtfSummariseResults(testGenes, method, Genes, GeneVars, TFs, models)

if nargin < 2,
  method='zscore';
end


if strcmp(method, 'margLogLik2')
  % Chib's approximatino to the marginal likelihood  
  results = zeros( size(testGenes) );
  % compute all marginal likelihoods for all trained models
  for k=1:size(testGenes,1) 
    fprintf('Running gene %d/%d...\n', k, size(testGenes, 1));
  for c=1:size(testGenes,2) 
  for f=1:size(testGenes,3)    
    %
    % approximate the entropy using Gaussian density estimation 
    if isfield(testGenes{k,c,f}, 'W')
       X = [testGenes{k,c,f}.kinetics; testGenes{k,c,f}.W; testGenes{k,c,f}.W0; testGenes{k,c,f}.sigma2]';
       % Store also a representative paramrters for the Chib;s marginal
       % likelihood (currently is the max in term of fititng the data)
       [mm mind] = max(testGenes{k,c,f}.LogL);
       parmu = [testGenes{k,c,f}.kinetics(:,mind); testGenes{k,c,f}.W(:,mind);...
               testGenes{k,c,f}.W0(1,mind); testGenes{k,c,f}.sigma2(1,mind)]';
       %X = [testGenes{k,c}.kinetics; testGenes{k,c}.W; testGenes{k,c}.W0]';
    else
       X = [testGenes{k,c,f}.kinetics; testGenes{k,c,f}.sigma2]'; 
       % Store also a representative paramrters for the Chib;s marginal
       % likelihood (currently is the max in term of fititng the data)
       [mm mind] = max(testGenes{k,c,f}.LogL);
       parmu = [testGenes{k,c,f}.kinetics(:,mind); testGenes{k,c,f}.sigma2(1,mind)]';
       %X = [testGenes{k,c}.kinetics]';
    end
    
    mu = mean(X,1);
    jit = 1e-8;
    Sigma = cov(X) + jit*eye(size(X,2));
    [N D] = size(X);                              
    logPosteriorOfSelectedPoint = - (0.5*D)*log(2*pi)  - 0.5*log(det(Sigma)) ...
                                  - 0.5*((parmu-mu)*(Sigma\(parmu-mu)'));
    
    % Chib's approximation to the marginal likelihood 
    %
    % use values with likelhiood under the model and evaluate the Chib's 
    % method there 
    if size(testGenes,3) > 1    
        modelTest = models{c,f};
    else
        modelTest = models{c};
    end
    modelTest.Likelihood.Genes = Genes(k,:,:);
    if isfield(modelTest.Likelihood.noiseModel, 'pumaSigma2') == 1
        modelTest.Likelihood.noiseModel.pumaSigma2(1,:,:) = GeneVars(k,:,:);
    end
    if modelTest.Likelihood.numTFs == 0
    %    
       modelTest.Likelihood.kinetics(1) = parmu(1);
       modelTest.Likelihood.kinetics(2) = parmu(2);
       modelTest.Likelihood.kinetics(3) = parmu(3);
       modelTest.Likelihood.noiseModel.sigma2 = parmu(4);
       if ~isfield(modelTest.Likelihood, 'crValMask') 
       modelTest.Likelihood.crValMask = 1:modelTest.Likelihood.numTimes;
       end
       % compute the maeginal likelhiod (GPs fucntino marginalized out)
       PredGenes = parmu(1)/parmu(2)  + (parmu(3) - parmu(1)/parmu(2))*exp(-modelTest.Likelihood.TimesG*parmu(2));
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
      modelTest.Likelihood.kinetics(1,:) = parmu(1:4);
      modelTest.Likelihood.W(1,:) = parmu(5:4+modelTest.Likelihood.numTFs);
      modelTest.Likelihood.W0 = parmu(5+modelTest.Likelihood.numTFs);
      modelTest.Likelihood.noiseModel.sigma2 = parmu(end);
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
    % evaluate the priors on the mu vectro
    lnpriorKin = ['ln',modelTest.prior.kinetics.type,'pdf'];
    TrspaceKin = modelTest.prior.kinetics.priorSpace; 
    Likkin = feval(TrspaceKin, modelTest.Likelihood.kinetics);
    LogPriorKin = feval(lnpriorKin, Likkin, modelTest.prior.kinetics);
      
    % white GP noise model is used
    lnpriorNoiseWHITE = ['ln',modelTest.prior.sigma2.type,'pdf'];
    TrspaceNoiseWHITE = modelTest.prior.sigma2.priorSpace; 
    temp = feval(TrspaceNoiseWHITE, modelTest.Likelihood.noiseModel.sigma2);
    LogPriorNoiseWHITE = feval(lnpriorNoiseWHITE, temp, modelTest.prior.sigma2);
    
    LogPrior = sum(LogPriorKin + LogPriorNoiseWHITE);
    
    if modelTest.Likelihood.numTFs > 0
    % evaluation of the prior for the interaction bias
    lnpriorW0 = ['ln',modelTest.prior.W0.type,'pdf'];
    TrspaceW0 = modelTest.prior.W0.priorSpace;  
    LikW0 = feval(TrspaceW0, modelTest.Likelihood.W0);
    LogPriorW0 = feval(lnpriorW0, LikW0, modelTest.prior.W0);
    % evaluation of the prior for the interaction weights
    lnpriorW = ['ln',modelTest.prior.W.type,'pdf'];
    TrspaceW = modelTest.prior.W.priorSpace; 
    LikW = feval(TrspaceW, modelTest.Likelihood.W);
    LogPriorW = feval(lnpriorW, LikW, modelTest.prior.W);
    %
    LogPrior = LogPrior + sum(LogPriorW0 + LogPriorW);
    %
    end
    %mm
    %logsumexp(LogLik)
    %max(LogLik)
    %modelTest.Likelihood.numTFs
    %pause
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
    % approximate the entropy using Gaussian density estimation 
    if isfield(testGenes{k,c,f}, 'W')
       X = [testGenes{k,c,f}.kinetics; testGenes{k,c,f}.W; testGenes{k,c,f}.W0; testGenes{k,c,f}.sigma2]';
       %X = [testGenes{k,c}.kinetics; testGenes{k,c}.W; testGenes{k,c}.W0]';
    else
       X = [testGenes{k,c,f}.kinetics; testGenes{k,c,f}.sigma2]'; 
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

