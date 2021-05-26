


%%%%%%%%%%%%%%  Load the test genes data  %%%%%%%%%%%%%%%% 
load /usr/local/michalis/mlprojects/gpsamp/matlab/datasets/drosophila_data;
load /usr/local/michalis/mlprojects/gpsamp/matlab/datasets/trainScaleDros;
load /usr/local/michalis/mlprojects/gpsamp/matlab/datasets/testset;
load /usr/local/michalis/mlprojects/gpsamp/matlab/drosTrainTotal;



% generate randomly 100 genes
numGenes = 100;
numTFs = 5;
TimesG = 0:11;
% model options 
options = gpmtfOptions(ones(1,12,3),numTFs);
options.constraints.initialConds = 0;
options.constraints.W0 = zeros(1,numTFs);

options.jointAct = 'sigmoid';
options.spikePriorW = 'yes';
% prior probablity for each interaction weight to be around zero 
% for each TF
TFpis = [0.0751 0.1163 0.1729  0.0378  0.2387];
options.spikepriors = 1 - [0.0751 0.1163 0.1729  0.0378  0.2387];

options.tauMax = 0; % no delays
% define the dense discretized grid in the time axis for the TF latent functions 
[options, TimesF] = gpmtfDiscretize(TimesG, options); 
modelTest = gpmtfCreate(ones(1,12), ones(1,12), [], [], TimesG, TimesF, options);

sigma2 = 0.01;

% precompute the TFs
for cnt=1:size(samples.F,2)
    modelTest.Likelihood.kineticsTF = samples.kineticsTF(:,:,cnt);
    for r=1:modelTest.Likelihood.numReplicas
        TFs{cnt}(:,:,r) = gpmtfComputeTF(modelTest.Likelihood, samples.F{cnt}(:,:,r), 1:numTFs);
    end
    TFs{cnt} = log(TFs{cnt});
end

if 1
for n=1:numGenes
    %
    % choose which TFs will be present 
    ch = rand(1,numTFs); 
    ch(ch<=TFpis)=1;
    ch(ch~=1)=0;
    if sum(ch(:)) == 0
        % at least TF should regulate
        ok = randperm(numTFs);
        ch(ok(1)) = 1; 
    end
    
    Net(n,:) = ch;
    
    % randomly select TFs profiles, interaction weigths and kinetics
    ch = randperm(size(samples.F,2));
    
    % only the first replica 
    TF = TFs{ch(1)}(:,:,1);
    TFindex(n) = ch(1); 
    
    LikParams = modelTest.Likelihood; 
    LikParams.TF = TF;
   
    Kinetics = 2*rand(1,4)-1;
    Kinetics = exp(Kinetics);%exp(rand(1,4));
   
    Kinetics(4) = Kinetics(1)/Kinetics(2);  
    LikParams.kinetics = Kinetics;
    
    W(n,:) = randn(1,numTFs).*Net(n,:);  
    W0(n) = 0;%randn;
    
    LikParams.W = W(n,:);
    LikParams.W0 = W0(n);
    
    Gns = gpmtfComputeGeneODE(LikParams, zeros(5,111), 1, 1);

    % add the some noise 
    Gns = Gns + sqrt(sigma2)*randn(size(Gns));
    
    Genes(n,:) = Gns(:, LikParams.comInds);;
    GenesVar(n,:) = sigma2*ones(1,12);
    GrTruth{n}.kinetics = LikParams.kinetics;
    GrTruth{n}.W = LikParams.W;
    GrTruth{n}.W0 = LikParams.W0;
    GrTruth{n}.TFindex = TFindex(n);
    %
end
end

mcmcoptions = mcmcOptions('controlPnts'); 
mcmcoptions.adapt.T = 40; 
mcmcoptions.adapt.Burnin = 40;
mcmcoptions.train.StoreEvery = 10;
mcmcoptions.train.T = 30000;
mcmcoptions.train.Burnin = 1000;


for n=1:size(Genes,1)

    TestGenes = Genes(n,:);
    TestGenesVar = GenesVar(n,:);
    % CREATE the model
    modelTest = gpmtfCreate(TestGenes, TestGenesVar, [], [], TimesG, TimesF, options);
    
    tic;
    [modelTest PropDist samplesTest accRates] = gpmtfTestGenesAdapt2(modelTest, TFs, [], mcmcoptions.adapt); 
    % training/sampling phase
    [modelTest PropDist samplesTest accRates] = gpmtfTestGenesSample2(modelTest, TFs, [], PropDist, mcmcoptions.train);
    toc;
    %
    testGene{n} = samplesTest;
    testaccRates{n} = accRates;
end

