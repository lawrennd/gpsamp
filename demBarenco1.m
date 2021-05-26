
dataName = 'p53Barenco5Genes'; %'toy3TfsDelays';
expNo = 1;
storeRes = 0;
printPlot = 0;

% load gene expression data for all replicas as well as the time points 
% (the expressions of each replica are assumed to have been obtained at
% common time points)
[Genes, TimesG, GenesVar, GenesTF, GenesTFVar, GroundTruth] = loadDataset(dataName);

% number of TFs 
numTFs = 1; 

% model options
options = gpmtfOptions(Genes,numTFs); 

options.tauMax = 0;
% define the dense discretized grid in the time axis for the TF latent functions 
[options, TimesF] = gpmtfDiscretize(TimesG, options); 

% CREATE the model
model = gpmtfCreate(Genes, GenesVar, GenesTF, GenesTFVar, TimesG, TimesF, options);

% store groundTruth parameters if available
if ~isempty(GroundTruth)
    model.groundtr = GroundTruth; 
end

mcmcoptions = mcmcOptions('controlPnts'); 
% adaption phase
[model PropDist samples accRates] = gpmtfAdapt(model, mcmcoptions.adapt);
% training/sampling phase
[model PropDist samples accRates] = gpmtfSample(model, PropDist, mcmcoptions.train);

if storeRes == 1
    d = date; 
    save(['dem' dataName num2str(experimentNo) d '.mat'], 'model','samples','accRates');
end

% Plot/Print the results 
gpmtfPlot(model, samples, dataName, printPlot);


