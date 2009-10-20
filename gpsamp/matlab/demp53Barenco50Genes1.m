
dataName = 'p53Barenco50Genes';
expNo = 1;
storeRes = 0;
printPlot = 1; 
% load gene expression data for all replicas as well as the time points 
% (the expressions of each replica are assumed to have been obtained at the same
%  time points)
[Genes, TimesG, GeneVars, TFGenes, TFGeneVars, GroundTruth] = loadDataset(dataName);

% number of TFs 
numTFs = 1; 

% model options
options = gpsampOptions(Genes,numTFs); 

options.tauMax = 4;
% define the dense discretized grid in the time axis for the TF latent functions 
[options, TimesF] = TFDiscretize(TimesG, options); 

% CREATE the model
if ~isempty(GeneVars) 
model = modelCreate(options, Genes, TimesG, TimesF, GeneVars);
else
model = modelCreate(options, Genes, TimesG, TimesF);
end

% store groundTruth parameters if available
if ~isempty(GroundTruth)
    model.groundtr = GroundTruth; 
end

mcmcoptions = mcmcOptions('controlPnts'); 
% adaption phase
[model PropDist samples accRates] = gpTFAdapt(model, mcmcoptions.adapt);
% training/sampling phase
[model PropDist samples accRates] = gpTFSample(model, PropDist, mcmcoptions.train);

if storeRes == 1
    d = date; 
    save(['dem' dataName num2str(experimentNo) d '.mat'], 'model','samples','accRates');
end

% Plot/Print the results 
plotCreate(model, samples, dataName, printPlot);
