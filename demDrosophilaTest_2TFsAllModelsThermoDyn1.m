% demDrosophilaTest_2TFsAllModels1 runs the multi-TF for screening
% using all possible combinations of two TFs(twist and mef2) 
function [Genes, GenesVar, TFs, models, mygenes, tGrid] = demDrosophilaTest_2TFsAllModels1(modulus, remainder, identifier, flag)

if nargin < 4,
  flag=1;
end
% set the grid of inverse temperatures 
tN = 30;
% power
p = 5; 

% if flag =0 , then jsut retutn the precomputations 
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

    % All possible models os twi and mef2
    comb = [0 0; 1 0; 0 1; 1 1];
    %
    TestGenes = Genes(1,:,:);
    TestGenesVar = GenesVar(1,:,:);
    numTFs = sum(comb(end,:));
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
    TFset = [3 5];
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

        tGrid = 0:tN;
        tGrid = (tGrid/tN).^p;
        % CREATE the model
        for ther=1:size(tGrid,2)
            modelTest = gpmtfCreate(TestGenes, TestGenesVar, [], [], TimesG, TimesF, options);
            modelTest.Likelihood.TFcomb = comb(c,:);
            modelTest.Likelihood.invT = tGrid(ther); 
            models{c,ther} = modelTest;
        end
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

    %outdir = '~/mlprojects/gpsamp/matlab/results';
    outdir = '/usr/local/michalis/mlprojects/gpsamp/matlab/results';
    outfile = sprintf('%s/multitf5a_%s_m%d_r%d.mat', outdir, identifier, modulus, remainder);

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

    % All possible models os twi and mef2
    comb = [0 0; 1 0; 0 1; 1 1];
    %

    models = {};
    testGene = {};
    if exist(outfile, 'file'),
        fprintf('Loading existing results from %s...\n', outfile);
        load(outfile);
        fprintf('Loaded results for %d genes.\n', size(testGene, 1));
    end

    %
    for n=2:size(Genes,1)
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
            if numTFs > 0
                %
                if comb(c,1) == 1
                    TFset = 3;
                    if comb(c,2) == 1
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

       
            tGrid =0:tN;
            tGrid = (tGrid/tN).^p;
            % run the thermodynamic integration loop 
            for ther=1:size(tGrid,2)
                modelTest = gpmtfCreate(TestGenes, TestGenesVar, [], [], TimesG, TimesF, options);
                % set inverse temperature
                modelTest.Likelihood.invT = tGrid(ther); 
                if numTFs > 0
                    [modelTest PropDist samplesTest accRates] = gpmtfTestGenesAdapt2(modelTest, TFs, mcmcoptions.adapt);
                    % training/sampling phase
                    [modelTest PropDist samplesTest accRates] = gpmtfTestGenesSample2(modelTest, TFs, PropDist, mcmcoptions.train);
                else
                    [modelTest PropDist samplesTest accRates] = gpmtfOnlyDecayModel(modelTest, mcmcoptions);
                end
                % 
                testGene{n, c, ther} = samplesTest;
                testaccRates{n, c, ther} = accRates; 
                modelTest.Likelihood.TFcomb = comb(c,:);
                models{c, ther} = modelTest;
            end
            %
        end
        % Only save after each gene is completed 
        save(outfile, 'testGene', 'testaccRates', 'mygenes', 'models');
        %safeSave(outfile, 'testGene', 'testaccRates', 'mygenes', 'models');
        %
        %
    end
    %safeSave(outfile, 'testGene', 'testaccRates', 'mygenes', 'models');
    fprintf('Completed.\n');
%
end