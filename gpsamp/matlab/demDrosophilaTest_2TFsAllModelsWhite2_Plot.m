% demDrosophilaTest_2TFsAllModelsWhite2_Plot: plots trained 
% models obtain by running  demDrosophilaTest_2TFsAllModelsWhite2_Plot
%demDrosophilaTest_2TFsAllModelsWhite2_Plot(modulus, remainder),

addpath ~/mlprojects/gpsamp/matlab
addpath ~/mlprojects/gpsamp/matlab/activFuncts
addpath ~/mlprojects/gpsamp/matlab/toolbox

modulus = 1;
remainder = 1;
%outdir = '~/mlprojects/gpsamp/matlab/results';
outdir = '/usr/local/michalis/mlprojects/gpsamp/matlab/results';
Date = '2010-03-09';
outfile = sprintf('%s/multitf3_%s_m%d_r%d.mat', outdir, Date, modulus, remainder);
load(outfile); 

dataName = 'Drosophila';
expNo = 1;
storeRes = 0;
printPlot = 0;

%%%%%%%%%%%%%%  Load the test genes data  %%%%%%%%%%%%%%%% 
%load ~/mlprojects/gpsamp/matlab/datasets/drosophila_data;
%load ~/mlprojects/gpsamp/matlab/datasets/trainScaleDros;
%load ~/mlprojects/gpsamp/matlab/datasets/testset;
%load ~/mlprojects/gpsamp/matlab/drosTrainTotal;
load datasets/drosophila_data;
load datasets/trainScaleDros;
load datasets/testset;
load drosTrainTotal;


if 1
load topranked10GenesMef2Twi;   
numGenes = 20;
indices = G;
mygenes = drosexp.genes(indices);
Genes = drosexp.fitmean(indices, :);
GenesVar = drosexp.fitvar(indices, :);
else
testindices = remainder:modulus:length(testset.indices);
indices = testset.indices(testindices);
numGenes = length(indices);
mygenes = drosexp.genes(indices);
Genes = drosexp.fitmean(indices, :);
GenesVar = drosexp.fitvar(indices, :);
end

% normalize separely the individucal genes
sc = 10./max(Genes, [], 2);
Genes = Genes.*repmat(sc, 1, size(Genes,2));
%Genes = Genes/sc;
Genes = reshape(Genes, numGenes, 12, 3);

GenesVar = GenesVar.*repmat(sc.^2, 1, size(GenesVar,2));
%GenesVar = GenesVar/(sc.^2);
GenesVar = reshape(GenesVar,numGenes, 12, 3);

TimesG = 0:11;
% All possible models os twi and mef2
comb = [0 0; 1 0; 0 1; 1 1];
%

TFset = [3 5];

if 1 %exist('models') == 0
%
   for c=1:size(comb,1)  
   %
       numTFs = sum(comb(c,:)); 
       % model options 
       options = gpmtfOptions(ones(1,12,3), numTFs); 
       options.jointAct = 'sigmoid';   
       %options.spikePriorW = 'yes';
       options.noiseModel = {'pumaWhite' 'white'};
       options.constraints.spaceW = 'positive'; 
       options.tauMax = 0; % no delays
       % define the dense discretized grid in the time axis for the TF latent functions 
       [options, TimesF] = gpmtfDiscretize(TimesG, options); 
       modelTest = gpmtfCreate(ones(1,12,3), ones(1,12,3), [], [], TimesG, TimesF, options);
       modelTest.Likelihood.TFcomb = comb(c,:);
       models{c} = modelTest;
    %
   end
%
end


% precompute the TFs
for cnt=1:size(samples.F,2)
%
    %load drosTrainTotal;
    numTFs = sum(comb(end,:)); 
    % model options 
    options = gpmtfOptions(ones(1,12,3), numTFs); 
    options.jointAct = 'sigmoid';   
    %options.spikePriorW = 'yes';
    options.noiseModel = {'pumaWhite' 'white'};
    options.constraints.spaceW = 'positive'; 
    options.tauMax = 0; % no delays
    % define the dense discretized grid in the time axis for the TF latent functions 
    [options, TimesF] = gpmtfDiscretize(TimesG, options); 
    modelTest = gpmtfCreate(ones(1,12,3), ones(1,12,3), [], [], TimesG, TimesF, options);
    modelTest.Likelihood.kineticsTF = samples.kineticsTF(TFset,:,cnt);
    for r=1:modelTest.Likelihood.numReplicas
         TFs{cnt}(:,:,r) = gpmtfComputeTFODE(modelTest.Likelihood, samples.F{cnt}(TFset,:,r), 1:numTFs);
    end    
    TFs{cnt} = log(TFs{cnt} + 1e-100);
%
end
%
  
%
for n=1:size(Genes,1) 
    gpmtfTestPlot(models, testGene(n,:), TFs, Genes(n,:,:), GenesVar(n,:, :), char(mygenes(n)), dataName, 1); close all;
end