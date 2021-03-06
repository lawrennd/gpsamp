% demDrosophilaTest_2TFsAllModelsCrossVal1 runs the multi-TF for screening
% using all possible combinations of two TFs(twist and mef2)
% It also used cros validation validatino by excluding each time 
% point at a time
<<<<<<< .mine
function [Genes, GenesVar, TFs, models, mygenes] = demDrosophilaTest_2TFsAllModelsCrossVal1(modulus, remainder, identifier, flag)
=======
function [Genes, GenesVar, TFs, models, mygenes] = demDrosophilaTest_2TFsAllModelsCrossVal1(modulus, remainder, identifier, flag),
>>>>>>> .r885

addpath ~/mlprojects/ndlutil/matlab
addpath ~/mlprojects/gpsamp/matlab
addpath ~/mlprojects/gpsamp/matlab/activFuncts
addpath ~/mlprojects/gpsamp/matlab/toolbox

if nargin < 3,
  identifier = datestr(now, 29);
end

if nargin < 4,
  flag = 1;
end

outdir = '~/mlprojects/gpsamp/matlab/results';
%outdir = '/usr/local/michalis/mlprojects/gpsamp/matlab/results';
%outdir = 'results/';
outfile = sprintf('%s/multitf5c_%s_m%d_r%d.mat', outdir, identifier, modulus, remainder);

dataName = 'drosophila_dataTest';
expNo = 1;
storeRes = 0;
printPlot = 0;
noiseM = {'pumaWhite' 'white'};

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
mygenes = drosexp.genes(G);
Genes = drosexp.fitmean(G, :);
GenesVar = drosexp.fitvar(G, :);
else
testindices = remainder:modulus:length(testset.indices);
indices = testset.indices(testindices);
numGenes = length(indices);
mygenes = drosexp.genes(indices);
Genes = drosexp.fitmean(indices, :);
GenesVar = drosexp.fitvar(indices, :);
end

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

if flag==0,
    TestGenes = Genes(1,:,:);
    TestGenesVar = GenesVar(1,:,:); 
    for m=1:size(comb,1)  
    %
      for cv=1:size(TimesG,2)
      %    
        numTFs = sum(comb(m,:)); 
        % model options 
        options = gpmtfOptions(ones(1,12,3), numTFs); 
        options.jointAct = 'sigmoid';   
        %options.spikePriorW = 'yes';
        options.noiseModel = noiseM;
        options.constraints.spaceW = 'positive'; 
        options.tauMax = 0; % no delays
        % define the dense discretized grid in the time axis for the TF latent functions
        [options, TimesF] = gpmtfDiscretize(TimesG, options); 
        modelTest = gpmtfCreate(ones(1,12,3), ones(1,12,3), [], [], TimesG, TimesF, options);  
        TFs = [];
        if numTFs > 0
        %
           if comb(m,1) == 1
              TFset = 3; 
              if comb(m,2) == 1
                 TFset = [3 5];
              end
           else
              TFset = 5; 
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
        modelTest.Likelihood.TFcomb = comb(m,:);
        modelTest.Likelihood.crValMask = [1:cv-1, cv+1:modelTest.Likelihood.numTimes];
        models{m,cv} = modelTest;
      end
    end

else
  models = {};
  testGene = {};
  if exist(outfile, 'file'),
    fprintf('Loading existing results from %s...\n', outfile);
    load(outfile);
    fprintf('Loaded results for %d genes.\n', size(testGene, 1));
  end 

  %
  for n=1:size(Genes,1)
  %    
    if size(testGene,1) > n,
        continue;
    end    
    fprintf('Running gene %d/%d...\n', n, size(Genes, 1));
    TestGenes = Genes(n,:,:);
    TestGenesVar = GenesVar(n,:,:); 
    for m=1:size(comb,1)  
    %
      for cv=1:size(TimesG,2)
      %    
        numTFs = sum(comb(m,:)); 
        % model options 
        options = gpmtfOptions(ones(1,12,3), numTFs); 
        options.jointAct = 'sigmoid';   
        %options.spikePriorW = 'yes';
        options.noiseModel = noiseM;
        options.constraints.spaceW = 'positive'; 
        options.tauMax = 0; % no delays
        % define the dense discretized grid in the time axis for the TF latent functions
        [options, TimesF] = gpmtfDiscretize(TimesG, options); 
        modelTest = gpmtfCreate(ones(1,12,3), ones(1,12,3), [], [], TimesG, TimesF, options);  
        TFs = [];
        if numTFs > 0
        %
           if comb(m,1) == 1
              TFset = 3; 
              if comb(m,2) == 1
                 TFset = [3 5];
              end
           else
              TFset = 5; 
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
        % This mask will zero the log ikelihood of the cvth data point 
        % which implies that the contribution of the cvth time point is p(y_cv|parameters)=1
        % and the data point is essentially excluded from the model
        modelTest.Likelihood.crValMask = [1:cv-1, cv+1:modelTest.Likelihood.numTimes];
        if numTFs > 0
           [modelTest PropDist samplesTest accRates] = gpmtfTestGenesAdapt2(modelTest, TFs, mcmcoptions.adapt); 
           % training/sampling phase
           [modelTest PropDist samplesTest accRates] = gpmtfTestGenesSample2(modelTest, TFs, PropDist, mcmcoptions.train);
        else
           [modelTest PropDist samplesTest accRates] = gpmtfOnlyDecayModel(modelTest, mcmcoptions);
        end
        %
        testGene{n, m, cv} = samplesTest;
        testaccRates{n, m, cv} = accRates; 
        modelTest.Likelihood.TFcomb = comb(c,:);
        models{m, cv} = modelTest;
        %
      end
      %
    end
    % Only save after each gene is completed
    safeSave(outfile, 'testGene', 'testaccRates', 'mygenes', 'models');
    %
  end
  %safeSave(outfile, 'testGene', 'testaccRates', 'mygenes', 'models');
  fprintf('Completed.\n');
end
