

dataName = 'drosophila_data';
expNo = 1;
storeRes = 0;
printPlot = 0;

% [Genes, TimesG, GenesVar, GenesTF, GenesTFVar, GroundTruth] = loadDataset(dataName); 

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
 
   %% Visualization: plot all the probes 
   %plot(drosexp.fitmean(I, :)','b'); 
   %hold on; 
   %% make red the selected one
   %plot(drosexp.fitmean(I(j), :)','r');
   %pause(0.1) 
   %hold off;
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
   
   I = find(strcmp(fbgnsTF(i), drosexp.fbgns));
  
   prScore = mean(drosexp.fitmean(I, :) ./ sqrt(drosexp.fitvar(I, :)), 2);
   
   % select one probe for the geneTF 
   [val j] = max(prScore);
   
   GenesTF = [GenesTF; drosexp.fitmean(I(j), :)];
   GenesTFVar = [GenesTFVar; drosexp.fitvar(I(j), :)];
  
   %% Visualization: plot all the probes 
   %plot(drosexp.fitmean(I, :)','b'); 
   %hold on; 
   %% make red the selected one
   %plot(drosexp.fitmean(I(j), :)','r');
   %pause 
   %hold off;
end

% scale the genes expression to roughly be
% in range [0 10]
numGenes = 92;
numTFs = 5; 


sc = 10./max(Genes, [], 2);
Genes = Genes.*repmat(sc, 1, size(Genes,2));
Genes = reshape(Genes, numGenes, 12, 3);
GenesVar = GenesVar.*repmat(sc.^2, 1, size(GenesVar,2));
GenesVar = reshape(GenesVar,numGenes, 12, 3);


%sc = 0.1*(max(Genes(:)) - min(Genes(:)));
%Genes = Genes/sc;
%Genes = reshape(Genes,numGenes,12,3);
%GenesVar = GenesVar/(sc.^2);
%GenesVar = reshape(GenesVar,numGenes,12,3);
%

sc = 10./max(GenesTF, [], 2);
GenesTF = GenesTF.*repmat(sc, 1, size(GenesTF,2));
GenesTF = reshape(GenesTF, numTFs, 12, 3);
GenesTFVar = GenesTFVar.*repmat(sc.^2, 1, size(GenesTFVar,2));
GenesTFVar = reshape(GenesTFVar, numTFs, 12, 3);

%GenesTF = GenesTF/sc;
%GenesTF = reshape(GenesTF,numTFs,12,3);
%GenesTFVar = GenesTFVar/(sc.^2);
%GenesTFVar = reshape(GenesTFVar,numTFs,12,3);

TimesG = 0:11;
%%%%%%%%%%%%%%  Load data  %%%%%%%%%%%%%%%% 

% model options
options = gpmtfOptions(Genes,numTFs); 
genesAndChip.data(genesAndChip.data~=0)=1;
options.constraints.X = genesAndChip.data; 
options.noiseModel = {'pumaWhite' 'white'};
%options.constraints.replicas = 'coupled'; 
%options.constraints.Ft0 = ones(1,numTFs);

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

% Plot/Print the results 
%gpmtfPlot(model, samples, dataName, printPlot);
