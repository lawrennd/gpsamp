
dataName = 'drosophila_dataTest';
expNo = 1;
storeRes = 0;
printPlot = 0;

%%%%%%%%%%%%%%  Load the test genes data  %%%%%%%%%%%%%%%% 
load datasets/drosophila_data;
load datasets/trainScaleDros;
%load datasets/testset;
load drosTrainTotal;

%testset.indices = [12197];
testset.indices = [12197, 2428, 1021, 983, 13538, 10898, 10339, 1958, 11320, 12182];
%numGenes = 500;
numGenes = length(testset.indices);
Genes = drosexp.fitmean(testset.indices(1:numGenes), :);
GenesVar = drosexp.fitvar(testset.indices(1:numGenes), :);

%sc = 10./max(Genes, [], 2);
%Genes = Genes.*repmat(sc, 1, size(Genes,2));
%Genes = reshape(Genes, numGenes, 12, 3);
%GenesVar = GenesVar.*repmat(sc.^2, 1, size(GenesVar,2));
%GenesVar = reshape(GenesVar,numGenes, 12, 3);

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
options.noiseModel = {'pumaWhite'};
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


%
for n=1:size(Genes,1)

    TestGenes = Genes(n,:,:);
    % !!! NO PUMA variancwes will be used. The model will samples the
    % variances !!!!!
    TestGenesVar = GenesVar(n,:,:); % [];
    % CREATE the model
    modelTest = gpmtfCreate(TestGenes, TestGenesVar, [], [], TimesG, TimesF, options);
   
    tic;
    [modelTest PropDist samplesTest accRates] = gpmtfTestGenesAdapt2(modelTest, TFs, mcmcoptions.adapt); 
    % training/sampling phase
    [modelTest PropDist samplesTest accRates] = gpmtfTestGenesSample2(modelTest, TFs, PropDist, mcmcoptions.train);
    toc;
    %
    testGene{n} = samplesTest;
    testaccRates{n} = accRates;
    % if you want to plot genes
    %gpmtfTestPlot(modelTest, Genes(n,:,:), GenesVar(n,:,:),  'n', testGene{n}, TFs, 'Dros', 0);

end
