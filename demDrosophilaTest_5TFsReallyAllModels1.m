% demDrosophilaTest_5TFsReallyAllModels1 runs the multi-TF for screening
% using all possible combinations of 5 TFs 
function [Genes, GenesVar, TFs, models, mygenes] = demDrosophilaTest_5TFsReallyAllModels1(modulus, remainder, identifier, flag)

if nargin < 4,
  flag=1;
end

% All possible models
comb = [0 0 0 0 0;
	1 0 0 0 0; 0 1 0 0 0; 0 0 1 0 0; 0 0 0 1 0; 0 0 0 0 1;
	1 1 0 0 0; 1 0 1 0 0; 1 0 0 1 0; 1 0 0 0 1;
	0 1 1 0 0; 0 1 0 1 0; 0 1 0 0 1;
	0 0 1 1 0; 0 0 1 0 1;
	0 0 0 1 1;
	1 1 1 0 0; 1 1 0 1 0; 1 1 0 0 1;
	1 0 1 1 0; 1 0 1 0 1;
	1 0 0 1 1;
	0 1 1 1 0; 0 1 1 0 1;
	0 1 0 1 1;
	0 0 1 1 1;
	1 1 1 1 0; 1 1 1 0 1; 1 1 0 1 1; 1 0 1 1 1; 0 1 1 1 1;
	1 1 1 1 1];

% if flag =0 , then just return the precomputations 
if flag == 0
%    
    addpath ~/mlprojects/ndlutil/matlab
    addpath ~/mlprojects/gpsamp/matlab
    addpath ~/mlprojects/gpsamp/matlab/activFuncts
    addpath ~/mlprojects/gpsamp/matlab/toolbox
    load datasets/drosophila_data;
    load drosTrainTotal;
    load datasets/testset;
    noiseM = {'pumaWhite' 'white'};

    if 0
        load topranked10GenesMef2Twi;
        numGenes = 20;
        mygenes = drosexp.genes(G);
        Genes = drosexp.fitmean(G, :);
        GenesVar = drosexp.fitvar(G, :);
    else
        testindices = remainder:modulus:length(testset.indices);
        indices = testset.indices(testindices);
        %indices = [2282 10486];
        %indices = 2282;
        %indices = 10486;
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
    TimesG = 0:11;

    %
    TestGenes = Genes(1,:,:);
    TestGenesVar = GenesVar(1,:,:);
    numTFs = size(samples.F{1},1);
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
    % precompute the TFs
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
        options.constraints.spaceW = 'positive';
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
    outfile = sprintf('%s/multitf6a_%s_m%d_r%d.mat', outdir, identifier, modulus, remainder);

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

    if 0
        load topranked10GenesMef2Twi;
        numGenes = 20;
        mygenes = drosexp.genes(G);
        Genes = drosexp.fitmean(G, :);
        GenesVar = drosexp.fitvar(G, :);
    else
        testindices = remainder:modulus:length(testset.indices);
        indices = testset.indices(testindices);
        %indices = [2282 10486];
        %indices = 2282;
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
    mcmcoptions.adapt.T = 80;
    mcmcoptions.adapt.Burnin = 1;
    mcmcoptions.adapt.disp = 0;
    mcmcoptions.train.StoreEvery = 12;
    mcmcoptions.train.T = 36000;
    mcmcoptions.train.Burnin = 1000;

    TimesG = 0:11;

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
            % Fix seeds
            randn('seed', 1e6);
            rand('seed', 1e6);
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
