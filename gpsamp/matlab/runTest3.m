function runTest3(modulus, reminder),

addpath ~/mlprojects/gpsamp/matlab
addpath ~/mlprojects/gpsamp/matlab/activFuncts
addpath ~/mlprojects/gpsamp/matlab/toolbox

outdir = '~/mlprojects/gpsamp/matlab/results';

outfile = sprintf('%s/multitf3_%s_m%d_r%d.mat', outdir, datestr(now, 29), modulus, reminder);

dataName = 'drosophila_dataTest';
expNo = 1;
storeRes = 0;
printPlot = 0;

%%%%%%%%%%%%%%  Load the test genes data  %%%%%%%%%%%%%%%% 
load ~/mlprojects/gpsamp/matlab/datasets/drosophila_data;
load ~/mlprojects/gpsamp/matlab/datasets/trainScaleDros;
load ~/mlprojects/gpsamp/matlab/datasets/testset;
load ~/mlprojects/gpsamp/matlab/drosTrainTotal;

testindices = reminder:modulus:length(testset.indices);
indices = testset.indices(testindices);
numGenes = length(indices);
%numGenes = 100;
mygenes = drosexp.genes(indices);
Genes = drosexp.fitmean(indices, :);
GenesVar = drosexp.fitvar(indices, :);

Genes = Genes/sc;
Genes = reshape(Genes,numGenes,12,3);
GenesVar = GenesVar/(sc.^2);
GenesVar = reshape(GenesVar,numGenes,12,3);
%
TimesG = 0:11;
numTFs = 3;
%%%%%%%%%%%%%%  Load data  %%%%%%%%%%%%%%%% 

mcmcoptions = mcmcOptions('controlPnts'); 
mcmcoptions.adapt.T = 40; 
mcmcoptions.adapt.Burnin = 40;
mcmcoptions.train.StoreEvery = 10;
mcmcoptions.train.T = 30000;
mcmcoptions.train.Burnin = 1000;

% model options 
options = gpmtfOptions(ones(1,12,3),numTFs); 
options.jointAct = 'sigmoid';
options.spikePriorW = 'yes';
options.constraints.spaceW = 'positive';
% prior probablity for each interaction weight to be around zero 
% for each TF
TFpis = [0.1163 0.1729  0.2387];
options.spikepriors = 1 - [0.1163 0.1729  0.2387];


options.tauMax = 0; % no delays
% define the dense discretized grid in the time axis for the TF latent functions 
[options, TimesF] = gpmtfDiscretize(TimesG, options); 
modelTest = gpmtfCreate(ones(1,12,3), ones(1,12,3), [], [], TimesG, TimesF, options);



% select the 3 most interesting TFs 
TFset = [2 3 5];
% precompute the TFs
for cnt=1:size(samples.F,2)
    modelTest.Likelihood.kineticsTF = samples.kineticsTF(TFset,:,cnt);
    for r=1:modelTest.Likelihood.numReplicas
        TFs{cnt}(:,:,r) = gpmtfComputeTFODE(modelTest.Likelihood, samples.F{cnt}(TFset,:,r), 1:numTFs);
    end
    samples.F{cnt} = samples.F{cnt}(TFset,:,r);
    
    TFs{cnt} = log(TFs{cnt} + 1e-100);
end
samples.kinetics = samples.kinetics(TFset, :);
samples.kineticsTF = samples.kineticsTF(TFset, :);

testGene = {};
if exist(outfile, 'file'),
  load(outfile);
end

%
for n=1:size(Genes,1)
  if length(testGene) >= n,
    continue;
  end
  
    TestGenes = Genes(n,:,:);
    % !!! NO PUMA variancwes will be used. The model will samples the
    % variances !!!!!
    TestGenesVar = [];
    % CREATE the model
    modelTest = gpmtfCreate(TestGenes, TestGenesVar, [], [], TimesG, TimesF, options);
   
    tic;
    [modelTest PropDist samplesTest accRates] = gpmtfTestGenesAdapt2(modelTest, TFs, [], mcmcoptions.adapt); 
    % training/sampling phase
    [modelTest PropDist samplesTest accRates] = gpmtfTestGenesSample2(modelTest, TFs, [], PropDist, mcmcoptions.train);
    toc;
    %
    if isfield(samplesTest, 'predGenes'),
      samplesTest = rmfield(samplesTest, 'predGenes');
    end
    testGene{n} = samplesTest;
    testaccRates{n} = accRates;
    save(outfile, 'testGene', 'testaccRates', 'mygenes');
    % if you want to plot genes
    %gpmtfTestPlot(modelTest, Genes(n,:,:), GenesVar(n,:,:),  'n', testGene{n}, TFs, 'Dros', 0);

end
