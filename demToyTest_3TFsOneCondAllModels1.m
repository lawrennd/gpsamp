% demToyTest_3TFsOneCondAllModels1.m runs the multi-TF for screening
% using all possible combinations of 3 TFs in toy data. One experimental
% condition is used 
function [Genes, GenesVar, TFs, models, mygenes] = demToyTest_3TFsOneCondAllModels1(modulus, remainder, identifier, flag)

if nargin < 4,
  flag=1;
end

% All possible models
comb = [0 0 0; 1 0 0; 0 1 0; 0 0 1;
         1 1 0; 1 0 1; 0 1 1;
         1 1 1];
        
% if flag =0 , then just retutn the precomputations 
if flag == 0
%    
    addpath ~/mlprojects/ndlutil/matlab
    addpath ~/mlprojects/gpsamp/matlab
    addpath ~/mlprojects/gpsamp/matlab/activFuncts
    addpath ~/mlprojects/gpsamp/matlab/toolbox
    load datasets/toy4TFs28_June_10.mat;
    load demtoy_dataOneCond108-Jul-2010.mat
    noiseM = {'white'};
    testsetIndices = (30+remainder):modulus:1030;
    Genes = Genes(testsetIndices, :, 1); 
    GenesVar = GenesVar(testsetIndices, :, 1);     
    numGenes = size(Genes,1);
    mygenes = testsetIndices;
    
    %
    TestGenes = Genes(1,:,:);
    TestGenesVar = GenesVar(1,:,:);
    numTFs = size(samples.F{1},1);
    % model options
    options = gpmtfOptions(ones(1,size(TimesG,2),3), numTFs);
    options.jointAct = 'sigmoid';
    %options.spikePriorW = 'yes';
    options.noiseModel = noiseM;
    %options.constraints.spaceW = 'positive';
    options.tauMax = 0; % no delays
    % define the dense discretized grid in the time axis for the TF latent functions
    [options, TimesF] = gpmtfDiscretize(TimesG, options);
    modelTest = gpmtfCreate(ones(size(TestGenes)), ones(size(TestGenesVar)), [], [], TimesG, TimesF, options);
    % precompute the TFs
    TFs =[];
    for cnt=1:size(samples.F,2)
        %
        modelTest.Likelihood.kineticsTF = samples.kineticsTF(:,:,cnt);
        for r=1:modelTest.Likelihood.numReplicas
            TFs{cnt}(:,:,r) = gpmtfComputeTFODE(modelTest.Likelihood, samples.F{cnt}(:,:,r), 1:numTFs);
        end
        TFs{cnt} = log(TFs{cnt} + 1e-100);
        %
    end
    for c=1:size(comb,1)
        %
        TFset = find(comb(c,:));
        numTFs = sum(comb(c,:));
        % model options
        options = gpmtfOptions(ones(1,12,3), numTFs);
        options.jointAct = 'sigmoid';
        %options.spikePriorW = 'yes';
        options.noiseModel = noiseM;
        %options.constraints.spaceW = 'positive';
        options.tauMax = 0; % no delays
        % define the dense discretized grid in the time axis for the TF latent functions
        [options, TimesF] = gpmtfDiscretize(TimesG, options);
        modelTest = gpmtfCreate(ones(1,12,3), ones(1,12,3), [], [], TimesG, TimesF, options);

        % CREATE the model
        modelTest = gpmtfCreate(TestGenes, TestGenesVar, [], [], TimesG, TimesF, options);
        modelTest.Likelihood.TFcomb = comb(c,:);
        models{c} = modelTest;
        %
    end   
%   
else % otherwise run the demo     
%  
    addpath ~/mlprojects/ndlutil/matlab
    addpath ~/mlprojects/gpsamp/matlab
    addpath ~/mlprojects/gpsamp/matlab/activFuncts
    addpath ~/mlprojects/gpsamp/matlab/toolbox

    if nargin < 3,
        identifier = datestr(now, 29);
    end

    outdir = '~/mlprojects/gpsamp/matlab/results';
    %outdir = '/usr/local/michalis/mlprojects/gpsamp/matlab/results';
    outfile = sprintf('%s/multitfToyOneCond6a_%s_m%d_r%d.mat', outdir, identifier, modulus, remainder);

    dataName = 'toy_dataOneCond';
    expNo = 1;
    storeRes = 0;
    printPlot = 0;
    noiseM = {'white'};

    %%%%%%%%%%%%%%  Load the test genes data  %%%%%%%%%%%%%%%%
    %load ~/mlprojects/gpsamp/matlab/datasets/drosophila_data;
    %load ~/mlprojects/gpsamp/matlab/datasets/trainScaleDros;
    %load ~/mlprojects/gpsamp/matlab/datasets/testset;
    %load ~/mlprojects/gpsamp/matlab/drosTrainTotal;
    load datasets/toy4TFs28_June_10.mat;
    load demtoy_dataOneCond108-Jul-2010.mat
    noiseM = {'white'};
    testsetIndices = (30+remainder):modulus:1030;
    Genes = Genes(testsetIndices, :, 1); 
    GenesVar = GenesVar(testsetIndices, :, 1);     
    numGenes = size(Genes,1);
    mygenes = testsetIndices;

    mcmcoptions = mcmcOptions('controlPnts');
    mcmcoptions.adapt.T = 40;
    mcmcoptions.adapt.Burnin = 40;
    mcmcoptions.adapt.disp = 0;
    mcmcoptions.train.StoreEvery = 10;
    mcmcoptions.train.T = 30000;
    mcmcoptions.train.Burnin = 1000;

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
        for c=1:size(comb,1)
            %
            numTFs = sum(comb(c,:));
            % model options
            options = gpmtfOptions(ones(size(TestGenes)), numTFs);
            options.jointAct = 'sigmoid';
            %options.spikePriorW = 'yes';
            options.noiseModel = noiseM;
            %options.constraints.spaceW = 'positive';
            options.tauMax = 0; % no delays
            % define the dense discretized grid in the time axis for the TF latent functions
            [options, TimesF] = gpmtfDiscretize(TimesG, options);
            modelTest = gpmtfCreate(ones(size(TestGenes)), ones(size(TestGenesVar)), [], [], TimesG, TimesF, options);

            TFs = [];
            TFset = find(comb(c,:));
            if numTFs > 0
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
            modelTest.Likelihood.TFcomb = comb(c,:);
            models{c} = modelTest;
            %
        end
        % Only save after each gene is completed 
        %save(outfile, 'testGene', 'testaccRates', 'mygenes', 'models');
        safeSave(outfile, 'testGene', 'testaccRates', 'mygenes', 'models');
        %
        %
    end
    %safeSave(outfile, 'testGene', 'testaccRates', 'mygenes', 'models');
    fprintf('Completed.\n');
%
end
