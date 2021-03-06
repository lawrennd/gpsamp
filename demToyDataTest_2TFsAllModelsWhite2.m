% demToyDataTest_2TFsAllModelsWhite2 runs the multi-TF for screening
% using all possible combinations of two TFs in toye -Drosophila data (twist and mef2)
% the model has been previously trained using 
% demToyDataTrain_2TFsWhite2
function demToyDataTest_2TFsAllModelsWhite2(modulus, remainder)

addpath ~/mlprojects/gpsamp/matlab
addpath ~/mlprojects/gpsamp/matlab/activFuncts
addpath ~/mlprojects/gpsamp/matlab/toolbox

%outdir = '~/mlprojects/gpsamp/matlab/results';
outdir = '/usr/local/michalis/mlprojects/gpsamp/matlab/results';

outfile = sprintf('%s/multitf3_%s_m%d_r%d.mat', outdir, datestr(now, 29), modulus, remainder);

dataName = 'Toy_dataTest';
expNo = 1;
storeRes = 0;
printPlot = 0;

%%%%%%%%%%%%%%  Load the test genes data  %%%%%%%%%%%%%%%% 
%load ~/mlprojects/gpsamp/matlab/datasets/drosophila_data;
%load ~/mlprojects/gpsamp/matlab/datasets/trainScaleDros;
%load ~/mlprojects/gpsamp/matlab/datasets/testset;
%load ~/mlprojects/gpsamp/matlab/drosTrainTotal; 
load demToyDataTrainedSamples27-Mar-2010.mat;
load datasets/toy5TFs25_March_10Final.mat;
testset.indices = 21:1000;
ok = Genes(testset.indices, :, :); 
Genes = ok(:,:,1); 
Genes = [Genes, ok(:,:,2)]; 
Genes = [Genes, ok(:,:,3)];
ok = GenesVar(testset.indices, :, :); 
GenesVar = ok(:,:,1); 
GenesVar = [GenesVar, ok(:,:,2)]; 
GenesVar = [GenesVar, ok(:,:,3)];


testindices = remainder:modulus:length(testset.indices);
indices = testset.indices(testindices);
numGenes = length(indices);
mygenes = indices;


% normalize separely the individucal genes
sc = 10./max(Genes, [], 2);
Genes = Genes.*repmat(sc, 1, size(Genes,2));
Genes = reshape(Genes, numGenes, 12, 3);
GenesVar = GenesVar.*repmat(sc.^2, 1, size(GenesVar,2));
GenesVar = reshape(GenesVar,numGenes, 12, 3);

mcmcoptions = mcmcOptions('controlPnts'); 
mcmcoptions.adapt.T = 40; 
mcmcoptions.adapt.Burnin = 40;
mcmcoptions.train.StoreEvery = 10;
mcmcoptions.train.T = 30000;
mcmcoptions.train.Burnin = 1000;

TimesG = 0:11;

% All possible models of twi and mef2
comb = [0 0; 1 0; 0 1; 1 1];
%

models = {};
testGene = {};
if exist(outfile, 'file'),
     load(outfile);
end 

%
for n=1:size(Genes,1)
%    
    if size(testGene,1) > n,
        continue;
    end
    
    TestGenes = Genes(n,:,:);
    TestGenesVar = []; % GenesVar(n,:,:); 
    %
    for c=1:size(comb,1)  
    %
        numTFs = sum(comb(c,:)); 
        % model options 
        options = gpmtfOptions(ones(1,12,3), numTFs); 
        options.jointAct = 'sigmoid';   
        %options.spikePriorW = 'yes';
        options.noiseModel = {'white'};
        %options.constraints.spaceW = 'positive'; 
        options.tauMax = 0; % no delays
        % define the dense discretized grid in the time axis for the TF latent functions 
        [options, TimesF] = gpmtfDiscretize(TimesG, options); 
        modelTest = gpmtfCreate(ones(1,12,3), ones(1,12,3), [], [], TimesG, TimesF, options);

        TFs = [];
        if numTFs > 0
        %
           if comb(c,1) == 1
              TFset = 1; 
              if comb(c,2) == 1
                 TFset = [1 2];
              end
           else
              TFset = 2; 
           end
           
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
        end
        
        % CREATE the model
        modelTest = gpmtfCreate(TestGenes, TestGenesVar, [], [], TimesG, TimesF, options);
   
        if numTFs > 0
           [modelTest PropDist samplesTest accRates] = gpmtfTestGenesAdapt2(modelTest, TFs, mcmcoptions.adapt); 
           % training/sampling phase
           [modelTest PropDist samplesTest accRates] = gpmtfTestGenesSample2(modelTest, TFs, PropDist, mcmcoptions.train);
        else
           [modelTest PropDist samplesTest accRates] = gpmtfOnlyDecayModel(modelTest, mcmcoptions);
        end
        %
        testGene{n,c} = samplesTest;
        testaccRates{n,c} = accRates; 
        models{c} = modelTest;
        models{c}.Likelihood.TFcomb = comb(c,:);
        save(outfile, 'testGene', 'testaccRates', 'mygenes', 'models');
        %
    end
   %
   %
end
save(outfile, 'testGene', 'testaccRates', 'mygenes', 'models');
