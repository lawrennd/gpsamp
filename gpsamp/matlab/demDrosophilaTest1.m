

dataName = 'drosophila_dataTest';
expNo = 1;
storeRes = 0;
printPlot = 0;


% This should go inside the loadDatasets function later
%%%%%%%%%%%%%%  Load the test genes data  %%%%%%%%%%%%%%%% 
%load datasets/drosophila_data.mat;
% load traning samples 
load drosTrainTotal;
load trainGenes;

TestGenes = Genes(1,:,:);
TestGenesVar = GenesVar(1,:,:);
TimesG = 0:11;
numTFs = 5;
%%%%%%%%%%%%%%  Load data  %%%%%%%%%%%%%%%% 

% model options
options = gpmtfOptions(TestGenes,numTFs); 

options.tauMax = 0; % no delays
% define the dense discretized grid in the time axis for the TF latent functions 
[options, TimesF] = gpmtfDiscretize(TimesG, options); 

% CREATE the model
modelTest = gpmtfCreate(TestGenes, TestGenesVar, [], [], TimesG, TimesF, options);

mcmcoptions = mcmcOptions('controlPnts'); 

% change the mcmc options
% MCMC OPTIONS (you can change these options if you like)   
mcmcoptions.train.StoreEvery = 1; % store samples after burn in  every StoreEvery iterations
mcmcoptions.train.Burnin = 1000;  % burn in time
mcmcoptions.train.T = 1; % keep only one sample after burn in  
        
% options for the adaptive phase in MCMC 
mcmcoptions.adapt.T = 40;          
mcmcoptions.adapt.Burnin = 10;
mcmcoptions.adapt.disp = 0;

% for each traininng sample for the TFs run a Markov chain to get an
% independent sample for the parameters of all test genes
% adaption phase
for cnt =1:size(samples.LogL,2)
    %
     % place the training sample in the model structure 
     modelTest.Likelihood.kineticsTF = samples.kineticsTF(:,:,cnt); 
     modelTest.F = samples.F{cnt};
     
     % precompute the TFs
     for r=1:modelTest.Likelihood.numReplicas
         modelTest.Likelihood.TF(:,:,r) = gpmtfComputeTF(modelTest.Likelihood,  modelTest.F(:,:,r), 1:numTFs);
     end
     
     tic;
     [modelTest PropDist samplesTest accRates] = gpmtfTestGenesAdapt(modelTest, mcmcoptions.adapt);
     % training/sampling phase
     [modelTest PropDist samplesTest accRates] = gpmtfTestGenesSample(modelTest, PropDist, mcmcoptions.train);
     toc;
     
     % store the collected single sample (obtained after burn-in) in a local variable
     ss.predGenes(:,:,cnt) = squeeze(samplesTest.predGenes{1})';
     ss.kinetics(:,cnt) = squeeze(samplesTest.kinetics); 
     ss.Weights(:,cnt) = squeeze(samplesTest.Weights); 
     ss.Weights0(:,cnt) = samplesTest.Weights0;
     ss.LogL(cnt) = samplesTest.LogL;
    %
end
