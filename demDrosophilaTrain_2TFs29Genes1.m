

dataName = 'drosophila_data';
expNo = 1;
storeRes = 0;
printPlot = 0;
normGenes = 1;
noiseM = {'pumaWhite', 'white'};

% This should go inside the loadDatasets function later
%%%%%%%%%%%%%%  Load data  %%%%%%%%%%%%%%%% 
load datasets/drosophila_data.mat;
genesAndChip = importdata('datasets/eileen_nature_training_set.txt'); 
TFset = [3 5]; %[1 2 3 4 5];   %% user spesificed quantity  
if ~isfield(drosexp, 'fbgns'),
  drosexp.fbgns = drosexp.genes;
end

fbgns = genesAndChip.textdata(2:end,1);
Genes = [];
GenesVar = [];
for i=1:size(fbgns,1) 
%
   I = find(strcmp(fbgns(i), drosexp.fbgns));
   prScore = mean(drosexp.fitmean(I, :) ./ sqrt(drosexp.fitvar(I, :)), 2);
   % select one probe for the gene 
   [val j] = max(prScore);  
   Genes = [Genes; drosexp.fitmean(I(j), :)];
   GenesVar = [GenesVar; drosexp.fitvar(I(j), :)];
% 
end

% collect the gens for the TFs
if isstruct(drosTF.fbgns),
  fbgnsTF = struct2cell(drosTF.fbgns);
else
  fbgnsTF = drosTF.fbgns; 
end
GenesTF = [];
GenesTFVar = [];
for i=1:size(fbgnsTF(:),1) 
%   
   I = find(strcmp(fbgnsTF(i), drosexp.fbgns));
   prScore = mean(drosexp.fitmean(I, :) ./ sqrt(drosexp.fitvar(I, :)), 2);  
   % select one probe for the geneTF 
   [val j] = max(prScore);
   GenesTF = [GenesTF; drosexp.fitmean(I(j), :)];
   GenesTFVar = [GenesTFVar; drosexp.fitvar(I(j), :)];
%   
end

% choose only the genes that are regulated only by one TF 
genesAndChip.data(genesAndChip.data~=0)=1;
ok = find(sum(genesAndChip.data(:, [1 2 4]),2)==0); 
Genes = Genes(ok,:); 
GenesVar = GenesVar(ok,:); 
numGenes = size(Genes,1);
chip = genesAndChip.data(ok,:); 
chip = chip(:, TFset);
numTFs = size(TFset,2); 
numGenes = size(Genes,1);


% separate are normalized seprately 
if normGenes == 1
%    
   % scale the genes expression to roughly be  
   % in range [0 10]
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
   
% gene are normalized globally
else
%    
   sc = 0.1*(max(Genes(:)) - min(Genes(:)));
   Genes = Genes/sc;
   Genes = reshape(Genes,numGenes,12,3);
   GenesVar = GenesVar/(sc.^2);
   GenesVar = reshape(GenesVar,numGenes,12,3);
   GenesTF = GenesTF/sc;
   GenesTF = reshape(GenesTF,5,12,3);
   GenesTFVar = GenesTFVar/(sc.^2);
   GenesTFVar = reshape(GenesTFVar,5,12,3);
   GenesTF = GenesTF(TFset, :, :);
   GenesTFVar = GenesTFVar(TFset, :, :);
%
end

TimesG = 0:11;

% model options
options = gpmtfOptions(Genes,numTFs); 
options.constraints.X = chip; 
options.noiseModel = noiseM; 

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
