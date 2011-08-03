
dataName = 'ToyData';
expNo = 1;
storeRes = 0;
printPlot = 0;

load datasets/toy5TFs25_March_10Final.mat;

TFset = [3 5];
chip = Net(1:20,TFset);
numTFs = size(TFset,2); 

testset.indices = 1:20;
ok = Genes(testset.indices, :, :); 
Genes = ok(:,:,1); 
Genes = [Genes, ok(:,:,2)]; 
Genes = [Genes, ok(:,:,3)];
ok = GenesVar(testset.indices, :, :); 
GenesVar = ok(:,:,1); 
GenesVar = [GenesVar, ok(:,:,2)]; 
GenesVar = [GenesVar, ok(:,:,3)];

numGenes = size(Genes,1);


% normalize separely the individucal genes
sc = 10./max(Genes, [], 2);
Genes = Genes.*repmat(sc, 1, size(Genes,2));
Genes = reshape(Genes, numGenes, 12, 3);
GenesVar = GenesVar.*repmat(sc.^2, 1, size(GenesVar,2));
GenesVar = reshape(GenesVar,numGenes, 12, 3);

sc = 10./max(GenesTF, [], 2);
GenesTF = GenesTF.*repmat(sc, 1, size(GenesTF,2));
GenesTF = reshape(GenesTF, 5, 12, 3);
GenesTFVar = GenesTFVar.*repmat(sc.^2, 1, size(GenesTFVar,2));
GenesTFVar = reshape(GenesTFVar, 5, 12, 3);

GenesTF = GenesTF(TFset, :, :);
GenesTFVar = GenesTFVar(TFset, :, :);
TimesG = 0:11;

%%%%%%%%%%%%%%  Load data  %%%%%%%%%%%%%%%% 

% model options
options = gpmtfOptions(Genes,numTFs); 
options.constraints.X = chip; 
options.noiseModel = {'white'};

options.tauMax = 0; % no delays
% define the dense discretized grid in the time axis for the TF latent functions 
[options, TimesF] = gpmtfDiscretize(TimesG, options); 

% CREATE the model
model = gpmtfCreate(Genes, GenesVar, GenesTF, GenesTFVar, TimesG, TimesF, options);

mcmcoptions = mcmcOptions('controlPnts'); 
mcmcoptions.adapt.T = 100;
mcmcoptions.train.StoreEvery = 50;
mcmcoptions.train.Burnin = 5000;
mcmcoptions.train.T = 50000;
% adaption phase
[model PropDist samples accRates] = gpmtfAdapt(model, mcmcoptions.adapt);
% training/sampling phase
[model PropDist samples accRates] = gpmtfSample(model, PropDist, mcmcoptions.train);

if storeRes == 1
    d = date; 
    save(['dem' dataName num2str(expNo) d '.mat'], 'model','samples','accRates');
end


