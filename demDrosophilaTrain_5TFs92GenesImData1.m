% (This demo used in the paper to produce the initial training with the 92 genes)
% DEMO that has created the file drosTrainTotal_28NOV2010.mat
dataName = 'drosophila_data';
expNo = 1;
storeRes = 0;
printPlot = 0;
normGenes = 3;
noiseM = {'pumaWhite' 'white'};

% This should go inside the loadDatasets function later
%%%%%%%%%%%%%%  Load data  %%%%%%%%%%%%%%%% 
load datasets/drosophila_data.mat;
genesAndChip = importdata('datasets/eileen_nature_training_set.txt'); 

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

% collect the genes for the TFs
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

numGenes = 92;
numTFs = 5; 
glbSc = 1;

% normalization based on Antti suggestion
if normGenes == 3
   for j=1:size(Genes,1)
       G = Genes(j,:);
       G = G(:);
       sc = mean(G.^2, 1); 
       G =Genes(j,:)/sqrt(sc); 
       Genes(j,:) = G;
       GenesVar(j,:) = GenesVar(j,:)/sc;
   end
   Genes = reshape(Genes, numGenes, 12, 3);
   GenesVar = reshape(GenesVar,numGenes, 12, 3);
   for j=1:size(GenesTF,1)
       G = GenesTF(j,:);
       G = G(:);
       sc = mean(G.^2, 1);
       G =GenesTF(j,:)/sqrt(sc); 
       GenesTF(j,:) = G;
       GenesTFVar(j,:) = GenesTFVar(j,:)/sc;
   end
   GenesTF = reshape(GenesTF, numTFs, 12, 3);
   GenesTFVar = reshape(GenesTFVar,numTFs, 12, 3);
% normalize so that each gene has unit variance 
elseif normGenes == 2
   for j=1:size(Genes,1)
       G = Genes(j,:);
       G = G(:); 
       sc = var(G);
       G =Genes(j,:)/sqrt(sc); 
       Genes(j,:) = G;
       GenesVar(j,:) = GenesVar(j,:)/sc;
   end
   Genes = reshape(Genes, numGenes, 12, 3);
   GenesVar = reshape(GenesVar,numGenes, 12, 3);
   for j=1:size(GenesTF,1)
       G = GenesTF(j,:);
       G = G(:);
       sc = var(G);
       G =GenesTF(j,:)/sqrt(sc); 
       GenesTF(j,:) = G;
       GenesTFVar(j,:) = GenesTFVar(j,:)/sc;
   end
   GenesTF = reshape(GenesTF, numTFs, 12, 3);
   GenesTFVar = reshape(GenesTFVar,numTFs, 12, 3);
% separate are normalized separately based on the maximum value
elseif normGenes == 1
%    
   sc = glbSc./max(Genes, [], 2);
   Genes = Genes.*repmat(sc, 1, size(Genes,2));
   Genes = reshape(Genes, numGenes, 12, 3);
   GenesVar = GenesVar.*repmat(sc.^2, 1, size(GenesVar,2));
   GenesVar = reshape(GenesVar,numGenes, 12, 3);
   sc = glbSc./max(GenesTF, [], 2);
   GenesTF = GenesTF.*repmat(sc, 1, size(GenesTF,2));
   GenesTF = reshape(GenesTF, numTFs, 12, 3);
   GenesTFVar = GenesTFVar.*repmat(sc.^2, 1, size(GenesTFVar,2));   
   GenesTFVar = reshape(GenesTFVar, numTFs, 12, 3);
%
% gene are normalized globally
else
%    
   sc = 0.1*(max(Genes(:)) - min(Genes(:)));
   Genes = Genes/sc;
   Genes = reshape(Genes,numGenes,12,3);
   GenesVar = GenesVar/(sc.^2);
   GenesVar = reshape(GenesVar,numGenes,12,3);
   GenesTF = GenesTF/sc;
   GenesTF = reshape(GenesTF,numTFs,12,3);
   GenesTFVar = GenesTFVar/(sc.^2);
   GenesTFVar = reshape(GenesTFVar,numTFs,12,3);
%
end

TimesG = 0:11;

% model options
options = gpmtfOptions(Genes,numTFs); 
genesAndChip.data(genesAndChip.data~=0)=1;
%
%selSubset = sum(genesAndChip.data,2)<=3;
%Genes = Genes(selSubset,:,:); 
%GenesVar = GenesVar(selSubset,:,:); 
%genesAndChip.data = genesAndChip.data(selSubset,:);
%fbgns = fbgns(selSubset);
%
options.constraints.X = genesAndChip.data; 
options.noiseModel = noiseM;
options.tauMax = 0; % no delays
options.lengthScalePrior = 'invGamma';
% robust Gaussian means: a two-component mixture with 
% a Gaussian  and a uniform distribution
%options.likelihoodModel = 'robustGaussian';
% define the dense discretized grid in the time axis for the TF latent functions 
[options, TimesF] = gpmtfDiscretize(TimesG, options); 

% Fix seeds
randn('seed', 1e6);
rand('seed', 1e6);
% CREATE the model
%options.constraints.spaceW = 'positive';
model = gpmtfCreate(Genes, GenesVar, GenesTF, GenesTFVar, TimesG, TimesF, options);
model.prior.lengthScale.constraint(1) = 2;
for j=1:model.Likelihood.numTFs
    model.GP{j}.lengthScale = model.prior.lengthScale.constraint(1) + 0.3;
end
model.Likelihood.kineticsTF(:,1) = 0.1;
model.constraints.Ft0(2) = 1;
%model.constraints.TFinitCond(2) = 1;

mcmcoptions = mcmcOptions('controlPnts'); 
mcmcoptions.adapt.T = 150;
mcmcoptions.adapt.Burnin = 50;
mcmcoptions.train.StoreEvery = 50;
mcmcoptions.train.Burnin = 5000;
mcmcoptions.train.T = 50000;
% adaption phase
[model PropDist samples accRates] = gpmtfAdaptImdata(model, mcmcoptions.adapt);
% training/sampling phase
[model PropDist samples accRates] = gpmtfSampleImdata(model, PropDist, mcmcoptions.train);

if storeRes == 1
    d = date; 
    save(['dem' dataName num2str(expNo) d '.mat'], 'model','samples','accRates');
end
