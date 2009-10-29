% create artificial regression datasets for several input dimensions and apply 
% the control-points algorithm. Keep the hyperparameters fixed and compute the KL
% divergence between exact the GP posterior Gaussian and the one obtained 
% from sampling 

% dimensions 
Dim = 10; 
Ntrain = 200;
% simple exponential kernel 
sigmaf = 1; % kernel variance   
sigma2 = 0.3^2; % nosie variance
ell2 = 0.01; % lengthscale
iter = 1; % how many times to perform the experiment  

% model options
options = gpsampOptions('regression');
mcmcoptions = mcmcOptions('controlPnts');
mcmcoptions.train.StoreEvery = 10;
mcmcoptions.train.T = 30000; 

for it = 1:iter 
%    
  for D=1:Dim 
    % 
    %set the the kernel and noise parameters in the log scale
    kernLogtheta = zeros(1,D+1);
    kernLogtheta(1:D) = log(ell2); 
    kernLogtheta(D+1) = log(sigmaf);
    likLogtheta(1) = log(sigma2);
 
    % randomly generate input data in the unit hypercube 
    Xtrain = rand(Ntrain,D); 
    % if one-dimensional, then sort the input data 
    Xtrain = sort(Xtrain); 
  
    % compute the covariance matrix 
    gptmp.logtheta = kernLogtheta;
    Knn = kernCompute(gptmp, Xtrain, Xtrain);
   
    % average correlation coefficient
    Correl(it,D) = (sum(Knn(:)) - size(Knn,1))/(Ntrain*(Ntrain-1)); 
   
    % generate a sample from the Gaussian process
    mu = zeros(1,Ntrain);
    F = gaussianSample(1, mu, Knn);
    % produce the output data by adding nosie to noise to F 
    Ytrain = F(:) + sqrt(sigma2)*randn(Ntrain,1);
   
    if D == 1
       plot(Xtrain,Ytrain,'+'); 
       hold on; 
       plot(Xtrain,F,'r');
       pause(1);
    end
  
    % set up the model to run MCMC 
    % create the model
    model = gpsampCreate(Ytrain, Xtrain, options);
    % fix the hyperparameter to ground-truth 
    model.GP.logtheta = kernLogtheta;
    model.Likelihood.logtheta = likLogtheta; 
    model.constraints.kernHyper = 'fixed';
    model.constraints.likHyper = 'fixed'; 
    mcmcoptions.adapt.incrNumContrBy = floor(D/3) + 1;
    [model PropDist samples accRates] = gpsampControlAdapt(model, mcmcoptions.adapt);
    [model PropDist samples accRates] = gpsampControlTrain(model, PropDist, mcmcoptions.train);
    samplesControl{D} = samples; 
    
    % compute KL divergence between full GP posteriot and MCMC Gaussian
    KLControl(it,D) = computeKLs(model, samples);
   
    data{D}.Xtrain = Xtrain;
    data{D}.Ytrain = Ytrain;
    numControl(it,D) = size(model.Xu,1);
    latentF{it}(D,:) = F(:)';
  %
end
end

