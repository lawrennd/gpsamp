


%%%%%%%%%%%%%%  Load the test genes data  %%%%%%%%%%%%%%%% 
load datasets/drosophila_data;
load datasets/trainScaleDros;
load datasets/testset;
load drosTrainTotal;



% generate randomly 100 genes
numGenes = 500;
numTFs = 3;
TimesG = 0:11;
% model options 
options = gpmtfOptions(ones(1,12,3),numTFs);
%options.constraints.initialConds = 0;
%options.constraints.W0 = zeros(1,numTFs);

options.jointAct = 'sigmoid';
options.spikePriorW = 'yes';
% prior probablity for each interaction weight to be around zero 
% for each TF
TFpis = [0.1163 0.1729  0.2387];
options.spikepriors = 1 - [0.1163 0.1729  0.2387];

options.tauMax = 0; % no delays
% define the dense discretized grid in the time axis for the TF latent functions 
[options, TimesF] = gpmtfDiscretize(TimesG, options); 
modelTest = gpmtfCreate(ones(1,12), ones(1,12), [], [], TimesG, TimesF, options);

sigma2 = 0.001;
TFset = [2 3 5];


% precompute the TFs
for cnt=1:size(samples.F,2)
    modelTest.Likelihood.kineticsTF = samples.kineticsTF(TFset,:,cnt);
    for r=1:modelTest.Likelihood.numReplicas
        TFs{cnt}(:,:,r) = gpmtfComputeTF(modelTest.Likelihood, samples.F{cnt}(TFset,:,r), 1:numTFs);
    end
    samples.F{cnt} = samples.F{cnt}(TFset,:,r);
    
    TFs{cnt} = log(TFs{cnt});
end
samples.kinetics = samples.kinetics(TFset, :);
samples.kineticsTF = samples.kineticsTF(TFset, :);

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
    W0(n) = randn;
    
    LikParams.W = W(n,:);
    LikParams.W0 = W0(n);
    
    Gns = gpmtfComputeGeneODE(LikParams, zeros(numTFs,111), 1, 1);

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

ss = 2*sqrt(modelTest.prior.weights.spikeSigma2);

for n=1:size(testGene,2)
    %
    for j=1:numTFs
        ok = testGene{n}.Weights(j,:);
        ok(ok>=-ss & ok<=ss) = 0; 
        ok(ok~=0)=1;
        prob(n,j) = sum(ok)/3000;
    end
    %
end

disp('Frequencies of TF (first 3 are the ground truth, last 3  are the learned ones)')
sum([Net, round(prob)])/numGenes
disp('Number of incorrectly predicted network connections out of numGenes x numTFs in total')
sum(sum(abs([Net-round(prob)])))

% visualization
for n=1:numGenes
    % 
    gpmtfTestPlot(modelTest, Genes(n,:,:), GenesVar(n,:,:),  'n', testGene{n}, TFs, 'toy', 0, GrTruth{n}); disp(n); 
    disp('close the windows and press any key to continue')    
    pause; 
    close all;
    %
end