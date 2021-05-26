% demDrosophilaTest_2TFsAllModelsWhite1 runs the multi-TF ffor screening
% using all possible combinations of two TFS(twist and mef2) 

dataName = 'drosophila_dataTest';
expNo = 1;
storeRes = 0;
printPlot = 0;

%%%%%%%%%%%%%%  Load the test genes data  %%%%%%%%%%%%%%%% 
load datasets/drosophila_data;
load datasets/trainScaleDros;
load datasets/testset;
load drosTrainTotal;

numGenes = 500;
Genes = drosexp.fitmean(testset.indices(1:numGenes), :);
GenesVar = drosexp.fitvar(testset.indices(1:numGenes), :);

sc = 10./max(Genes, [], 2);
Genes = Genes.*repmat(sc, 1, size(Genes,2));
%Genes = Genes/sc;
Genes = reshape(Genes, numGenes, 12, 3);

GenesVar = GenesVar.*repmat(sc.^2, 1, size(GenesVar,2));
%GenesVar = GenesVar/(sc.^2);
GenesVar = reshape(GenesVar,numGenes, 12, 3);

mcmcoptions = mcmcOptions('controlPnts'); 
mcmcoptions.adapt.T = 40; 
mcmcoptions.adapt.Burnin = 40;
mcmcoptions.train.StoreEvery = 10;
mcmcoptions.train.T = 30000;
mcmcoptions.train.Burnin = 1000;

%
TimesG = 0:11;

% All possible models: 1) twist only, 2) mef2 only 3) twist and mef2
comb = [1 0; 0 1; 1 1]; 

for c=1:size(comb,1)
%

   numTFs = sum(comb(c,:)); 
 
   % model options 
   options = gpmtfOptions(ones(1,12,3),numTFs); 
   options.jointAct = 'sigmoid';   
   %options.spikePriorW = 'yes';
   options.noiseModel = {'pumaWhite' 'white'};
   %options.constraints.spaceW = 'positive';
   %
   % prior probablity for each interaction weight to be around zero 
   % for each TF
   %TFpis = [0.1163 0.1729 0.2387];
   %options.spikepriors = 1 - [0.1163 0.1729  0.2387];

   options.tauMax = 0; % no delays
   % define the dense discretized grid in the time axis for the TF latent functions 
   [options, TimesF] = gpmtfDiscretize(TimesG, options); 
   modelTest = gpmtfCreate(ones(1,12,3), ones(1,12,3), [], [], TimesG, TimesF, options);

   if comb(c,1) == 1
      TFset = 3; 
      if comb(c,2) == 1
         TFset = [3 5];
      end
   else
      TFset = 5; 
   end

   TFs = [];
   % precompute the TFs
   for cnt=1:size(samples.F,2)
       %
       modelTest.Likelihood.kineticsTF = samples.kineticsTF(TFset,:,cnt);
       for r=1:modelTest.Likelihood.numReplicas
           TFs{cnt}(:,:,r) = gpmtfComputeTFODE(modelTest.Likelihood, samples.F{cnt}(TFset,:,r), 1:numTFs);
       end    
       TFs{cnt} = log(TFs{cnt} + 1e-100);
       %
   end

   %
   for n=1:size(Genes,1)
    %
       TestGenes = Genes(n,:,:);
       % cheat the puma variances
       TestGenesVar = GenesVar(n,:,:); 
       % CREATE the model
       modelTest = gpmtfCreate(TestGenes, TestGenesVar, [], [], TimesG, TimesF, options);
   
       tic;
       [modelTest PropDist samplesTest accRates] = gpmtfTestGenesAdapt2(modelTest, TFs, mcmcoptions.adapt); 
       % training/sampling phase
       [modelTest PropDist samplesTest accRates] = gpmtfTestGenesSample2(modelTest, TFs, PropDist, mcmcoptions.train);
       toc;
       %
       testGene{c,n} = samplesTest;
       testaccRates{c,n} = accRates;
       % if you want to plot genes
       %gpmtfTestPlot(modelTest, Genes(n,:,:), GenesVar(n,:,:),  'n', testGene{n}, TFs, 'Dros', 0);
    %
   end
   %
   %
end
