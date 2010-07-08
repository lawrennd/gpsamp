
dataName = 'toy_dataTwoConds';
expNo = 1;
storeRes = 1;
printPlot = 0;
normGenes = 0;
noiseM = {'white'};

numGenes = 30;
numTFs = 3; 
load datasets/toy4TFs28_June_10.mat; 
Genes = Genes(1:30,:,:); 
GenesTF = GenesTF(1:3,:,:); 

% model options
options = gpmtfOptions(Genes,numTFs); 
options.constraints.X = Net(1:30,1:3); 
options.noiseModel = noiseM;
%options.constraints.Ft0 = ones(1,numTFs);
options.tauMax = 0; % no delays

% define the dense discretized grid in the time axis for the TF latent functions 
[options, TimesF] = gpmtfDiscretize(TimesG, options); 

% CREATE the model
model = gpmtfCreate(Genes, [], GenesTF, [], TimesG, TimesF, options);

mcmcoptions = mcmcOptions('controlPnts'); 
mcmcoptions.adapt.T = 100;
mcmcoptions.train.StoreEvery = 100;
mcmcoptions.train.Burnin = 5000;
mcmcoptions.train.T = 100000;
% adaption phase
[model PropDist samples accRates] = gpmtfAdaptImdata(model, mcmcoptions.adapt);
% training/sampling phase
[model PropDist samples accRates] = gpmtfSampleImdata(model, PropDist, mcmcoptions.train);

if storeRes == 1
    d = date; 
    save(['dem' dataName num2str(expNo) d '.mat'], 'model','samples','accRates');
end

% Plot/Print the results 
%gpmtfPlot(model, samples, dataName, printPlot);
