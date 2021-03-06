demBarenco1.m                                                                                       0000700 0003466 0000024 00000002342 11300762661 013062  0                                                                                                    ustar   mtitsias                        games                                                                                                                                                                                                                  
dataName = 'p53Barenco5Genes'; %'toy3TfsDelays';
expNo = 1;
storeRes = 0;
printPlot = 0;

% load gene expression data for all replicas as well as the time points 
% (the expressions of each replica are assumed to have been obtained at
% common time points)
[Genes, TimesG, GenesVar, GenesTF, GenesTFVar, GroundTruth] = loadDataset(dataName);

% number of TFs 
numTFs = 1; 

% model options
options = gpmtfOptions(Genes,numTFs); 

options.tauMax = 0;
% define the dense discretized grid in the time axis for the TF latent functions 
[options, TimesF] = gpmtfDiscretize(TimesG, options); 

% CREATE the model
model = gpmtfCreate(Genes, GenesVar, GenesTF, GenesTFVar, TimesG, TimesF, options);

% store groundTruth parameters if available
if ~isempty(GroundTruth)
    model.groundtr = GroundTruth; 
end

mcmcoptions = mcmcOptions('controlPnts'); 
% adaption phase
[model PropDist samples accRates] = gpmtfAdapt(model, mcmcoptions.adapt);
% training/sampling phase
[model PropDist samples accRates] = gpmtfSample(model, PropDist, mcmcoptions.train);

if storeRes == 1
    d = date; 
    save(['dem' dataName num2str(experimentNo) d '.mat'], 'model','samples','accRates');
end

% Plot/Print the results 
gpmtfPlot(model, samples, dataName, printPlot);


                                                                                                                                                                                                                                                                                              demBinaryClass1.m                                                                                   0000700 0003466 0000024 00000002076 11275534044 013732  0                                                                                                    ustar   mtitsias                        games                                                                                                                                                                                                                  
% illustration in 1-d classification

dataName = 'binaryclassToy';
expNo = 1;
storeRes = 0;


% CREATE DATA 
N = 200;
% simple exponential kernel 
sigmaf = 1; % kernel variance   
sigma2 = 0.03^2; % noise variance
ell2 = 0.01; % lengthscale
D = 1;
% randomly generate input data 
X = rand(N,D); 
X = sort(X);  
gptmp.logtheta = [log(ell2) log(sigmaf)];
Knn = kernCompute(gptmp, X, X);
F = gaussianSample(1, zeros(1,N), Knn);
Y = F(:) + sqrt(sigma2)*randn(N,1);
Y(Y<0)=-1;
Y(Y>=0)=1;

% model options
options = gpsampOptions('classification'); 

% create the model
model = gpsampCreate(Y, X, options);

mcmcoptions = mcmcOptions('controlPnts');
mcmcoptions.adapt.incrNumContrBy = 1;
% adaption phase
[model PropDist samples accRates] = gpsampControlAdapt(model, mcmcoptions.adapt);
% training/sampling phase
[model PropDist samples accRates] = gpsampControlTrain(model, PropDist, mcmcoptions.train);

if storeRes == 1
    d = date; 
    save(['dem' dataName num2str(experimentNo) d '.mat'], 'model', 'PropDist','samples','accRates');
end

% plot the samples
gpsampPlot(model, samples);
                                                                                                                                                                                                                                                                                                                                                                                                                                                                  demDrosophila1.m                                                                                    0000600 0003466 0000024 00000006203 11312241627 013611  0                                                                                                    ustar   mtitsias                        games                                                                                                                                                                                                                  

dataName = 'drosophila_data';
expNo = 1;
storeRes = 0;
printPlot = 0;

% [Genes, TimesG, GenesVar, GenesTF, GenesTFVar, GroundTruth] = loadDataset(dataName); 

% This should go inside the loadDatasets function later
%%%%%%%%%%%%%%  Load data  %%%%%%%%%%%%%%%% 
load datasets/drosophila_data.mat;
genesAndChip = importdata('datasets/eileen_nature_training_set.txt'); 

if ~isfield(drosexp, 'fbgns'),
  drosexp.fbgns = drosexp.genes;
end

fbgns = genesAndChip.textdata(2:end,1);
Genes = [];
GenesVar = [];
for i=1:size(fbgns,1) 
   %
   I = find(strcmp(fbgns(i), drosexp.fbgns));
   
   prScore = mean(drosexp.fitmean(I, :) ./ sqrt(drosexp.fitvar(I, :)), 2);
   
   % select one probe for the gene 
   [val j] = max(prScore);
   
   Genes = [Genes; drosexp.fitmean(I(j), :)];
   GenesVar = [GenesVar; drosexp.fitvar(I(j), :)];
 
   %% Visualization: plot all the probes 
   %plot(drosexp.fitmean(I, :)','b'); 
   %hold on; 
   %% make red the selected one
   %plot(drosexp.fitmean(I(j), :)','r');
   %pause(0.1) 
   %hold off;
end

% collect the gens for the TFs
if isstruct(drosTF.fbgns),
  fbgnsTF = struct2cell(drosTF.fbgns);
else
  fbgnsTF = drosTF.fbgns; 
end
GenesTF = [];
GenesTFVar = [];
for i=1:size(fbgnsTF(:),1) 
   
   I = find(strcmp(fbgnsTF(i), drosexp.fbgns));
  
   prScore = mean(drosexp.fitmean(I, :) ./ sqrt(drosexp.fitvar(I, :)), 2);
   
   % select one probe for the geneTF 
   [val j] = max(prScore);
   
   GenesTF = [GenesTF; drosexp.fitmean(I(j), :)];
   GenesTFVar = [GenesTFVar; drosexp.fitvar(I(j), :)];
  
   %% Visualization: plot all the probes 
   %plot(drosexp.fitmean(I, :)','b'); 
   %hold on; 
   %% make red the selected one
   %plot(drosexp.fitmean(I(j), :)','r');
   %pause 
   %hold off;
end

% scale the genes expression to roughly be
% in range [0 10]
numGenes = 92;
numTFs = 5; 
sc = 0.1*(max(Genes(:)) - min(Genes(:)));
Genes = Genes/sc;
Genes = reshape(Genes,numGenes,12,3);
GenesVar = GenesVar/(sc.^2);
GenesVar = reshape(GenesVar,numGenes,12,3);
%
GenesTF = GenesTF/sc;
GenesTF = reshape(GenesTF,numTFs,12,3);
GenesTFVar = GenesTFVar/(sc.^2);
GenesTFVar = reshape(GenesTFVar,numTFs,12,3);
TimesG = 0:11;
%%%%%%%%%%%%%%  Load data  %%%%%%%%%%%%%%%% 

% model options
options = gpmtfOptions(Genes,numTFs); 
genesAndChip.data(genesAndChip.data~=0)=1;
options.constraints.X = genesAndChip.data; 
%options.constraints.replicas = 'coupled'; 
%options.constraints.Ft0 = ones(1,numTFs);

options.tauMax = 0; % no delays
% define the dense discretized grid in the time axis for the TF latent functions 
[options, TimesF] = gpmtfDiscretize(TimesG, options); 

% CREATE the model
model = gpmtfCreate(Genes, GenesVar, GenesTF, GenesTFVar, TimesG, TimesF, options);

mcmcoptions = mcmcOptions('controlPnts'); 
mcmcoptions.adapt.T = 100;
mcmcoptions.train.StoreEvery = 10;
% adaption phase
[model PropDist samples accRates] = gpmtfAdapt(model, mcmcoptions.adapt);
% training/sampling phase
[model PropDist samples accRates] = gpmtfSample(model, PropDist, mcmcoptions.train);

if storeRes == 1
    d = date; 
    save(['dem' dataName num2str(expNo) d '.mat'], 'model','samples','accRates');
end

% Plot/Print the results 
%gpmtfPlot(model, samples, dataName, printPlot);
                                                                                                                                                                                                                                                                                                                                                                                             demDrosophilaTest1.m                                                                                0000600 0003466 0000024 00000004561 11311763014 014454  0                                                                                                    ustar   mtitsias                        games                                                                                                                                                                                                                  

dataName = 'drosophila_dataTest';
expNo = 1;
storeRes = 0;
printPlot = 0;


% This should go inside the loadDatasets function later
%%%%%%%%%%%%%%  Load the test genes data  %%%%%%%%%%%%%%%% 
%load datasets/drosophila_data.mat;
% load traning samples 
load drosTrainTotal;
load trainGenes;

TestGenes = Genes(1,:,:);
TestGenesVar = GenesVar(1,:,:);
TimesG = 0:11;
numTFs = 5;
%%%%%%%%%%%%%%  Load data  %%%%%%%%%%%%%%%% 

% model options
options = gpmtfOptions(TestGenes,numTFs); 

options.tauMax = 0; % no delays
% define the dense discretized grid in the time axis for the TF latent functions 
[options, TimesF] = gpmtfDiscretize(TimesG, options); 

% CREATE the model
modelTest = gpmtfCreate(TestGenes, TestGenesVar, [], [], TimesG, TimesF, options);

mcmcoptions = mcmcOptions('controlPnts'); 

% change the mcmc options
% MCMC OPTIONS (you can change these options if you like)   
mcmcoptions.train.StoreEvery = 1; % store samples after burn in  every StoreEvery iterations
mcmcoptions.train.Burnin = 1000;  % burn in time
mcmcoptions.train.T = 1; % keep only one sample after burn in  
        
% options for the adaptive phase in MCMC 
mcmcoptions.adapt.T = 40;          
mcmcoptions.adapt.Burnin = 10;
mcmcoptions.adapt.disp = 0;

% for each traininng sample for the TFs run a Markov chain to get an
% independent sample for the parameters of all test genes
% adaption phase
for cnt =1:size(samples.LogL,2)
    %
     % place the training sample in the model structure 
     modelTest.Likelihood.kineticsTF = samples.kineticsTF(:,:,cnt); 
     modelTest.F = samples.F{cnt};
     
     % precompute the TFs
     for r=1:modelTest.Likelihood.numReplicas
         modelTest.Likelihood.TF(:,:,r) = gpmtfComputeTF(modelTest.Likelihood,  modelTest.F(:,:,r), 1:numTFs);
     end
     
     tic;
     [modelTest PropDist samplesTest accRates] = gpmtfTestGenesAdapt(modelTest, mcmcoptions.adapt);
     % training/sampling phase
     [modelTest PropDist samplesTest accRates] = gpmtfTestGenesSample(modelTest, PropDist, mcmcoptions.train);
     toc;
     
     % store the collected single sample (obtained after burn-in) in a local variable
     ss.predGenes(:,:,cnt) = squeeze(samplesTest.predGenes{1})';
     ss.kinetics(:,cnt) = squeeze(samplesTest.kinetics); 
     ss.Weights(:,cnt) = squeeze(samplesTest.Weights); 
     ss.Weights0(:,cnt) = samplesTest.Weights0;
     ss.LogL(cnt) = samplesTest.LogL;
    %
end
                                                                                                                                               demDrosophilaTest2.m                                                                                0000600 0003466 0000024 00000010266 11321374142 014455  0                                                                                                    ustar   mtitsias                        games                                                                                                                                                                                                                  
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

Genes = Genes/sc;
Genes = reshape(Genes,numGenes,12,3);
GenesVar = GenesVar/(sc.^2);
GenesVar = reshape(GenesVar,numGenes,12,3);

%%!!!!
%numGenes = 92;
%load datasets/dros92proprocessedTraininggenes;
%%!!!!

%
TimesG = 0:11;
numTFs = 5;
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
<<<<<<< .mine
options.constraints.spaceW = 'positive';
=======
options.constraints.spaceW = '';
>>>>>>> .r676
% prior probablity for each interaction weight to be around zero 
% for each TF
options.spikepriors = 1 - [0.0751 0.1163 0.1729  0.0378  0.2387];
%options.constraints.replicas = 'coupled'; 

options.tauMax = 0; % no delays
% define the dense discretized grid in the time axis for the TF latent functions 
[options, TimesF] = gpmtfDiscretize(TimesG, options); 
modelTest = gpmtfCreate(ones(1,12,3), ones(1,12,3), [], [], TimesG, TimesF, options);

% precompute the TFs
if strcmp(modelTest.constraints.replicas,'free')
    for cnt=1:size(samples.F,2)
        modelTest.Likelihood.kineticsTF = samples.kineticsTF(:,:,cnt);
        for r=1:modelTest.Likelihood.numReplicas
            TFs{cnt}(:,:,r) = gpmtfComputeTFODE(modelTest.Likelihood, samples.F{cnt}(:,:,r), 1:numTFs);
        end
        TFs{cnt} = log(TFs{cnt} + 1e-100);
    end
else 
    % replicas are coupled 
    for cnt=1:size(samples.F,2)
        modelTest.Likelihood.kineticsTF = samples.kineticsTF(:,:,cnt);
        TFs{cnt} = gpmtfComputeTFODE(modelTest.Likelihood, samples.F{cnt}, 1:numTFs);
        TFs{cnt} = log(TFs{cnt} + 1e-100);
    end
    %
end


% process the test genes one at a time 
%simMat = [];
%
for n=1:size(Genes,1)

    TestGenes = Genes(n,:,:);
    TestGenesVar = GenesVar(n,:,:);
    
    %options.constraints.X = Chip(n,:); 
    
    % CREATE the model
    modelTest = gpmtfCreate(TestGenes, TestGenesVar, [], [], TimesG, TimesF, options);
   
    % % compute all distances between samples (useful for sampling using the
    % % geometric distribution)
    % if strcmp(modelTest.constraints.replicas,'free')
    %     for j=1:modelTest.Likelihood.numTFs
    %     for r=1:modelTest.Likelihood.numReplicas
    %         % collect all TFs  of the r replica in matrix
    %         MatTF = zeros(size(samples.F,2), size(modelTest.Likelihood.TimesF,2));
    %         for cnt=1:size(samples.F,2)
    %             MatTF(cnt,:)  = TFs{cnt}(j,:,r);
    %         end
    %
    %         % compute the similarity matrix
    %         dd = dist2(MatTF,MatTF);
    %         for cnt=1:size(samples.F,2)
    %            [ok Ind] = sort(dd(cnt,:));
    %            simMat{j,r}(cnt,:) = uint16(Ind(2:end));
    %         end
    %     end
    %     end
    %     %
    % else
    %     %
    %     for j=1:modelTest.Likelihood.numTFs
    %     MatTF = zeros(size(samples.F,2), size(modelTest.Likelihood.TimesF,2));
    %     for cnt=1:size(samples.F,2)
    %         MatTF(cnt,:)  = TFs{cnt}(j,:);
    %     end
    %
    %     % compute the similarity matrix
    %     dd = dist2(MatTF,MatTF);
    %     for cnt=1:size(samples.F,2)
    %        [ok Ind] = sort(dd(cnt,:));
    %        simMat{j}(cnt,:) = uint16(Ind(2:end));
    %     end
    %     end
    %     %
    % end
    %
    %% parameter in the geometric distribution
    %modelTest.geo = 0.1;

    tic;
    [modelTest PropDist samplesTest accRates] = gpmtfTestGenesAdapt2(modelTest, TFs, [], mcmcoptions.adapt); 
    % training/sampling phase
    [modelTest PropDist samplesTest accRates] = gpmtfTestGenesSample2(modelTest, TFs, [], PropDist, mcmcoptions.train);
    toc;
    %
    testGene{n} = samplesTest;
    testaccRates{n} = accRates;
end



                                                                                                                                                                                                                                                                                                                                          demDrosophilaTest3.m                                                                                0000600 0003466 0000024 00000005141 11320651171 014451  0                                                                                                    ustar   mtitsias                        games                                                                                                                                                                                                                  
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


%
for n=1:size(Genes,1)

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
    testGene{n} = samplesTest;
    testaccRates{n} = accRates;
    % if you want to plot genes
    %gpmtfTestPlot(modelTest, Genes(n,:,:), GenesVar(n,:,:),  'n', testGene{n}, TFs, 'Dros', 0);

end
                                                                                                                                                                                                                                                                                                                                                                                                                               demDrosophilaTest3Single.m                                                                          0000600 0003466 0000024 00000005214 11321200452 015605  0                                                                                                    ustar   mtitsias                        games                                                                                                                                                                                                                  
dataName = 'drosophila_dataTest';
expNo = 1;
storeRes = 0;
printPlot = 0;

%%%%%%%%%%%%%%  Load the test genes data  %%%%%%%%%%%%%%%% 
load datasets/drosophila_data;
load datasets/trainScaleDros;
%load datasets/testset;
load drosTrainTotal;

testset.indices = [12197];
%numGenes = 500;
numGenes = 1;
Genes = drosexp.fitmean(testset.indices(1:numGenes), :);
GenesVar = drosexp.fitvar(testset.indices(1:numGenes), :);

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


%
for n=1:size(Genes,1)

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
    testGene{n} = samplesTest;
    testaccRates{n} = accRates;
    % if you want to plot genes
    %gpmtfTestPlot(modelTest, Genes(n,:,:), GenesVar(n,:,:),  'n', testGene{n}, TFs, 'Dros', 0);

end
                                                                                                                                                                                                                                                                                                                                                                                    demDrosophilaTest4.m                                                                                0000600 0003466 0000024 00000005244 11321165516 014462  0                                                                                                    ustar   mtitsias                        games                                                                                                                                                                                                                  
% run only in the top-10 ranked genes for of mef2 and twist from the
% single-Tf ranking

dataName = 'drosophila_dataTest';
expNo = 1;
storeRes = 0;
printPlot = 0;

%%%%%%%%%%%%%%  Load the test genes data  %%%%%%%%%%%%%%%% 
load datasets/drosophila_data;
load datasets/trainScaleDros;
load datasets/testset;
load topranked10GenesMef2Twi;
load drosTrainTotal;

numGenes = 20;
Genes = drosexp.fitmean(G, :);
GenesVar = drosexp.fitvar(G, :);

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


%
for n=1:size(Genes,1)

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
    testGene{n} = samplesTest;
    testaccRates{n} = accRates;
    % if you want to plot genes
    %gpmtfTestPlot(modelTest, Genes(n,:,:), GenesVar(n,:,:),  'n', testGene{n}, TFs, 'Dros', 0);

end
                                                                                                                                                                                                                                                                                                                                                            demMakeTestPlot1.m                                                                                  0000600 0003466 0000024 00000003056 11321171453 014063  0                                                                                                    ustar   mtitsias                        games                                                                                                                                                                                                                  
load datasets/drosophila_data.mat;
load datasets/trainScaleDros;

printPlots = 0; 

G = [];

disp('mef2 Top-ranked genes');
for n=1:size(mef2_r,2)
   % 
   %
   I = find(strcmp(mef2_r{n}.genes, drosexp.fbgns));
   prScore = mean(drosexp.fitmean(I, :) ./ sqrt(drosexp.fitvar(I, :)), 2);   
   % select one probe  
   [val j] = max(prScore);
   I = I(j); 
   
   
   Gene = drosexp.fitmean(I, :)/sc;
   Gene = reshape(Gene,1,12,3);
   GeneVar = drosexp.fitvar(I, :)/(sc.^2); 
   GeneVar = reshape(GeneVar,1,12,3);
   
   % make the plot 
   gpmtfTestPlot(modelTest, Gene, GeneVar,  mef2_r{n}.genes, mef2_r{n}.testGene, TFs, 'dros', printPlots);
   
   fprintf(1,'Gene #%2d Average log likelihood %5f\n',n, mean(mef2_r{n}.testGene.LogL));     
   disp('Press any key to continue');
   
   pause;
   
   G = [G, I]; 
   close all; 
   
   %
   %
end


for n=1:size(twi_r, 2)
   %
   %
   
   I = find(strcmp(twi_r {n}.genes, drosexp.fbgns));
   
   prScore = mean(drosexp.fitmean(I, :) ./ sqrt(drosexp.fitvar(I, :)), 2);   
   % select one probe  
   [val j] = max(prScore);
   I = I(j); 
   
   Gene = drosexp.fitmean(I, :)/sc;
   Gene = reshape(Gene,1,12,3);
   GeneVar = drosexp.fitvar(I, :)/(sc.^2); 
   GeneVar = reshape(GeneVar,1,12,3);
      
   % make the plot 
   gpmtfTestPlot(modelTest, Gene, GeneVar,  twi_r {n}.genes, twi_r{n}.testGene, TFs, 'dros', printPlots);
   
   
   fprintf(1,'Gene #%2d Average log likelihood %5f\n',n, mean(twi_r{n}.testGene.LogL));     
   disp('Press any key to continue');
   
   pause;
   
   G = [G, I]; 
    
   close all; 
   %
   %
end


                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  demRegress1.m                                                                                       0000700 0003466 0000024 00000001632 11275534060 013125  0                                                                                                    ustar   mtitsias                        games                                                                                                                                                                                                                  
dataName = 'regressionToy';
expNo = 1;
storeRes = 0;

% load outputs and inputs
% create some data or load data from somewhere
%[Y, X] = loadDataset(dataName)
X = (0:0.01:1)';
sigma2 = 0.3^2; 
Y = sin(2*pi*X) + cos((5*pi*X)) + sqrt(sigma2)*randn(101,1);
plot(X, Y, '+k', 'lineWidth', 3); 
[n D] = size(X); 
T = size(Y,2);


% model options
options = gpsampOptions('regression'); 

% create the model
model = gpsampCreate(Y, X, options);

mcmcoptions = mcmcOptions('controlPnts');
mcmcoptions.adapt.incrNumContrBy = 1;
% adaption phase
[model PropDist samples accRates] = gpsampControlAdapt(model, mcmcoptions.adapt);
% training/sampling phase
[model PropDist samples accRates] = gpsampControlTrain(model, PropDist, mcmcoptions.train);

if storeRes == 1
    d = date; 
    save(['dem' dataName num2str(experimentNo) d '.mat'], 'model', 'PropDist','samples','accRates');
end

% plot the samples
gpsampPlot(model, samples);
                                                                                                      demRegressArtifIncrDimension.m                                                                      0000600 0003466 0000024 00000004741 11275534103 016515  0                                                                                                    ustar   mtitsias                        games                                                                                                                                                                                                                  % create artificial regression datasets for several input dimensions and apply 
% the control-points algorithm. Keep the hyperparameters fixed and compute the KL
% divergence between exact the GP posterior Gaussian and the one obtained 
% from sampling 

% dimensions 
Dim = 10; 
Ntrain = 200;
% simple exponential kernel 
sigmaf = 1; % kernel variance   
sigma2 = 0.3^2; % nosie variance
ell2 = 0.01; % lengthscale
iter = 1; % how many times to perform the experiment  

% model options
options = gpsampOptions('regression');
mcmcoptions = mcmcOptions('controlPnts');
mcmcoptions.train.StoreEvery = 10;
mcmcoptions.train.T = 30000; 

for it = 1:iter 
%    
  for D=1:Dim 
    % 
    %set the the kernel and noise parameters in the log scale
    kernLogtheta = zeros(1,D+1);
    kernLogtheta(1:D) = log(ell2); 
    kernLogtheta(D+1) = log(sigmaf);
    likLogtheta(1) = log(sigma2);
 
    % randomly generate input data in the unit hypercube 
    Xtrain = rand(Ntrain,D); 
    % if one-dimensional, then sort the input data 
    Xtrain = sort(Xtrain); 
  
    % compute the covariance matrix 
    gptmp.logtheta = kernLogtheta;
    Knn = kernCompute(gptmp, Xtrain, Xtrain);
   
    % average correlation coefficient
    Correl(it,D) = (sum(Knn(:)) - size(Knn,1))/(Ntrain*(Ntrain-1)); 
   
    % generate a sample from the Gaussian process
    mu = zeros(1,Ntrain);
    F = gaussianSample(1, mu, Knn);
    % produce the output data by adding nosie to noise to F 
    Ytrain = F(:) + sqrt(sigma2)*randn(Ntrain,1);
    
    if D == 1
       plot(Xtrain,Ytrain,'+'); 
       hold on; 
       plot(Xtrain,F,'r');
       pause(1);
    end
  
    % set up the model to run MCMC 
    % create the model
    model = gpsampCreate(Ytrain, Xtrain, options);
    % fix the hyperparameter to ground-truth 
    model.GP.logtheta = kernLogtheta;
    model.Likelihood.logtheta = likLogtheta; 
    model.constraints.kernHyper = 'fixed';
    model.constraints.likHyper = 'fixed';
    
    mcmcoptions.adapt.incrNumContrBy = floor(D/3) + 1;
    [model PropDist samples accRates] = gpsampControlAdapt(model, mcmcoptions.adapt);
    [model PropDist samples accRates] = gpsampControlTrain(model, PropDist, mcmcoptions.train);
    samplesControl{D} = samples; 
    
    % compute KL divergence between full GP posteriot and MCMC Gaussian
    KLControl(it,D) = computeKLs(model, samples);
   
    data{D}.Xtrain = Xtrain;
    data{D}.Ytrain = Ytrain;
    numControl(it,D) = size(model.Xu,1);
    latentF{it}(D,:) = F(:)';
    mm{D} = model; 
  %
end
end

                               demRunSparseRegression.m                                                                            0000600 0003466 0000024 00000001075 11312735634 015417  0                                                                                                    ustar   mtitsias                        games                                                                                                                                                                                                                  tfs = {drosTF.names, {'mef2', 'twi', 'bin'}, {'mef2', 'twi', 'tin'}};
priorinds = {1:5, [5, 3, 2], [5, 3, 1]};

for l=1:length(tfs),
  fbgns = {};
  for k=1:length(tfs{l}),
    fbgns{k} = drosTF.fbgns.(tfs{l}{k});
  end
  inputs = drosFindGeneinds(drosexp, fbgns);
  
  regression_res{l} = zeros(length(testset.indices), 2^length(inputs));
  regression_ass{l} = zeros(length(testset.indices), 5);

  for k=1:length(testset.indices),
    [regression_res{l}(k, :), regression_ass{l}(k, :)] = drosSparseRegression(drosexp, inputs, testset.indices(k), priorinds{l});
  end
end
                                                                                                                                                                                                                                                                                                                                                                                                                                                                   demTestModel1.m                                                                                     0000600 0003466 0000024 00000006455 11312753652 013424  0                                                                                                    ustar   mtitsias                        games                                                                                                                                                                                                                  


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

                                                                                                                                                                                                                   demTestModel2.m                                                                                     0000600 0003466 0000024 00000010122 11320465413 013401  0                                                                                                    ustar   mtitsias                        games                                                                                                                                                                                                                  


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
%options.constraints.W0 = 0;
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
    
    W(n,:) = abs(randn(1,numTFs)).*Net(n,:);  
    W0(n) = options.constraints.W0.*randn;
    
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
    TestGenesVar = [];
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
end                                                                                                                                                                                                                                                                                                                                                                                                                                              drosMakeValidationMatrix.m                                                                          0000600 0003466 0000024 00000000705 11313700570 015702  0                                                                                                    ustar   mtitsias                        games                                                                                                                                                                                                                  function M = drosMakeValidationMatrix(chipdistances, genes, threshold),

tforder = {'tin', 'bin', 'twi', 'bap', 'mef2'};

if nargin < 3,
  threshold = 2000;
end

I = zeros(size(tforder));
for k=1:length(tforder),
  I(k) = find(~cellfun('isempty', strfind(chipdistances.labels, tforder{k})));
end

J = drosFindGeneinds(chipdistances, genes, 1);

M = zeros(length(J), length(I));
M(J==0, :) = NaN;
M(J~=0, :) = chipdistances.data(J(J~=0), I) < threshold;
                                                           drosPlotValidation.m                                                                                0000600 0003466 0000024 00000004745 11321625060 014565  0                                                                                                    ustar   mtitsias                        games                                                                                                                                                                                                                  %zthresholds = 0:.1:9;
%for t=1:length(zthresholds),
%  zacc(t) = nanmean(M(results.zscores' > zthresholds(t)));
%end

pplrthresholds = 0:.1:1;
for t=1:length(pplrthresholds),
  val = M(results.pplrs' > pplrthresholds(t));
  pplracc(t) = nanmean(val);
  ppval(t) = 1 - hygecdf(nansum(val)-1, sum(~isnan(M(:))), nansum(M(:)), sum(~isnan(val)));
end

scorethresholds = 10.^[-5:.5:0];
for t=1:length(scorethresholds),
  val = M(scores > scorethresholds(t));
  scoreacc(t) = nanmean(val);
  spval(t) = 1 - hygecdf(nansum(val)-1, sum(~isnan(M(:))), nansum(M(:)), sum(~isnan(val)));
end

subplot(4, 2, 1);

loglog(scorethresholds, spval)
xlabel('score threshold')
title('p-values')
%plot(pplrthresholds, pplracc)
%xlabel('pplr score threshold')
%ylabel('frequency of positive ChIP validation')

subplot(4, 2, 2);
semilogx(scorethresholds, scoreacc)
xlabel('score threshold')
title('accuracies')

% subplot(4, 3, 3);
% regthresholds = 0.9:.001:1;
% for t=1:length(regthresholds),
%   val = M(regression_probs > regthresholds(t));
%   regacc(t) = nanmean(val);
%   rpval(t) = 1 - hygecdf(nansum(val)-1, sum(~isnan(M(:))), nansum(M(:)), sum(~isnan(val)));
% end

% plot(regthresholds, regacc)
% xlabel('regression probability threshold')
% axis tight

TFs = {'bin', 'twi', 'mef2'};

for k=1:3,
  %for t=1:length(zthresholds),
  %  zaccs{k}(t) = nanmean(M(results.zscores(k, :) > zthresholds(t), k));
  %end

  for t=1:length(scorethresholds),
    val = M(scores(:, k) > scorethresholds(t), k);
    scoreaccs{k}(t) = nanmean(val);
    pvals{k}(t) = 1 - hygecdf(nansum(val)-1, sum(~isnan(M(:, k))), nansum(M(:, k)), sum(~isnan(val)));
  end
%  for t=1:length(pplrthresholds),
%    val = M(results.pplrs(k, :) > pplrthresholds(t), k);
%    pplraccs{k}(t) = nanmean(val);
%    pvals{k}(t) = 1 - hygecdf(nansum(val)-1, sum(~isnan(M(:, k))), nansum(M(:, k)), sum(~isnan(val)));
%  end

  subplot(4, 2, 1 + 2*k);
  
  loglog(scorethresholds, pvals{k})
  xlabel('score threshold')
  ylabel(TFs{k})
  %ylabel('frequency of positive ChIP validation')

  subplot(4, 2, 2 + 2*k);

  semilogx(scorethresholds, scoreaccs{k})
  xlabel('score threshold')
  axis tight
  %ylabel('frequency of positive ChIP validation')

%   subplot(4, 3, 1 + 3*k);
  
%   semilogy(pplrthresholds, pvals{k})
%   xlabel('pplr score threshold')
%   ylabel(TFs{k})
%   %ylabel('frequency of positive ChIP validation')

%   subplot(4, 3, 2 + 3*k);

%   plot(pplrthresholds, pplraccs{k})
%   xlabel('pplr score threshold')
%   %ylabel('frequency of positive ChIP validation')

end
                           drosSparseRegression.m                                                                              0000600 0003466 0000024 00000002207 11312735634 015132  0                                                                                                    ustar   mtitsias                        games                                                                                                                                                                                                                  function [res, mlass] = drosSparseRegression(drosexp, inputs, output, priorindices),

spikepriors = 1 - [0.0751 0.1163 0.1729  0.0378  0.2387];

spikepriors = spikepriors(priorindices);

F = drosexp.fitmean(inputs, :)';
y = drosexp.fitmean(output, :)';
yvar = drosexp.fitvar(output, :)';

F = [F ./ repmat(sqrt(var(F)), [size(F, 1), 1]), ones(length(y), 1)];
yscale = sqrt(var(y));
y = y ./ yscale;
yvar = yvar ./ (yscale .^ 2);

A = F' * diag(1./yvar) * F;
Ainv = inv(A);
b = sum(F .* repmat(y ./ yvar, [1, length(inputs) + 1]))';
c = .5 * sum(y.^2 ./ yvar);

lls = zeros(1, 2^length(inputs));

for k=0:(2^length(inputs)-1),
  bitmask = bitget(uint32(k), length(inputs):-1:1);
  prior = double(bitmask) .* (1-spikepriors) + (double(1-bitmask) .* spikepriors);
  sigma_w = [double(bitmask).*1 + double((1-bitmask)).* 0.0001, 1];
  sigmainv = diag(1 ./ sigma_w) + A;
  lls(k+1) = -.5 * log(det(sigmainv)) - .5 * sum(log(sigma_w)) ...
      +.5 * b'*(sigmainv\b) + sum(log(prior));
end

res = exp(lls - max(lls));
res = res ./ sum(res);

[foo, k] = max(res);
mlass0 = bitget(uint32(k-1), length(inputs):-1:1);
mlass = zeros(1, 5);
mlass(priorindices) = mlass0;
                                                                                                                                                                                                                                                                                                                                                                                         gpmtfAdapt.m                                                                                        0000700 0003466 0000024 00000046062 11312474512 013036  0                                                                                                    ustar   mtitsias                        games                                                                                                                                                                                                                  function [model PropDist samples accRates] = gpmtfAdapt(model, AdaptOps)
%[model PropDist samples accRates] = gpmtfAdapt(model, AdaptOps)
%
% Description: Sample the parameters of the Bayesian differential equation 
%              model so that to tune/adapt the proposal distribution (number 
%              of control points, variances of Gaussian proposals etc). 
%                ** This function should always be called before 
%                   gpTFSample. When the gpTFSample will be called the 
%                   proposal distribution (previously adapted)
%                   is kept fixed. **
%
% Inputs: 
%    -- model: the structure that contains the likelihood and GP
%              parameters as well as the priors for all these
%              quantities; write 'help model' for full details. 
%    -- PropDist: a stucture that defines the functional form of the 
%                 proposal distribution. It contains the following fields
%              * qF: contains precomputed variables needed for the 
%                    definition of the poposal distribution (that uses 
%                    control variables) over the GP latent functions 
%              * qKinVars: NumOfGenes x 4 matrix that contains the
%                    variances of the Gaussian proposal distributions
%                    used to sample the logarithm of the kinetic parameters 
%                    (represented in the log space; see modelCreate)
%                    for each gene. These Gaussian proposals have diagonal 
%                    covariance matrices so as each row corresponds to each 
%                    gene, 1 column to the basal rates (B), 2 column to 
%                    the decay rates (D), 3 column to the sensitivities(S)
%                    and 4 column to the initial condition parameter (A). 
%              * qWeigVars: NumOfGenes x M matrix that contains the
%                    variances of the Gaussian proposal distribution
%                    of all the parameters that exist in the activation 
%                    function of the latent GP function (that represent the 
%                    log of the TFs). Each row of the matrix corresponds to 
%                    each gene. The number M of parameters depends on the 
%                    number of TFs used and the functional form of the 
%                    activation function. When the sigmoid activation 
%                    function is used then  M=NumOfTFs+1. When the 
%                    Michaelis-Menten that is valid only for a single 
%                    TF (NumOfTFs=1), then M = 1; the single parameter 
%                    corresponds to the gamma Michaelis-Menten constant
%              * qLengScVars: NumOfTFs x 1 vector that contains all the
%                    variances of the Gaussian proposal distributions used to 
%                    sample the logarithm of the lengthscales of the NumOfTFs
%                    different rbf GP priors (see modelCreate) 
%    -- Genes : NumOfGenes x NumOfTimes x Replicas that stores the
%               gene expressions for all genes, all times and 
%               replicas
%    -- TimesG: The time points where gene expression are evaluated 
%    -- TimesF: The times where the GP function are evaluated
%               TimesF>>TimesG
%    -- trainOps: User defined options about the burn-in and sampling
%               iterations, visualization features and others 
%               **see demos**
%
% Outputs: 
%    -- model: see above or write 'help model' for full details. 
%              The outputed model is updated to contain the parameters 
%              values of the final MCMC iteration
%    -- PropDist: as above. the precomputations in the PropDist.qF field
%              can be different (compared to the ones given at the input) 
%              due the update of the kernel lengthscale that determine 
%              the proposal over the latent GP functions
%    -- samples: the structure that contains the samples. In contains
%              
%              
%
%    -- accRateF: acceptance rates during sampling time (after burn-in) 
%                 for the TFs
%    -- accRateKin: >>  >> for the kinetic parameters 
%    -- accRateW: >>  >> for the interaction weigths and any other 
%              parameters that exist in the activation function of 
%              the TF (e.g. the gamma in the Michaels-Menten activation)  
%    -- accRateLengSc:  >> >> for the lengthscales of the Gaussian
%              kernel 

BurnInIters = AdaptOps.Burnin; 
Iters = AdaptOps.T; 

%%% SIZES 
[NumOfGenes SizG NumReplicas] = size(model.Likelihood.Genes);
% number of genes 
NumOfGenes = model.Likelihood.numGenes;
% number of times the gene expression is evaluated
SizG = model.Likelihood.numTimes;
% number of replicas per gene 
NumReplicas = model.Likelihood.numReplicas;
% number of transcription factors
NumOfTFs = model.Likelihood.numTFs;
SizKin = size(model.Likelihood.kinetics,2);
SizKinTF = 2;
SizF = size(model.Likelihood.TimesF,2);

% initial number of control variables 
M = AdaptOps.initialNumContrPnts;

% initialize the transcription factors
if strcmp(model.constraints.replicas,'free')
    %
    F = zeros(NumOfTFs, SizF, NumReplicas);
    % initialize the control variables 
    Fu = zeros(NumOfTFs, M, NumReplicas);
    for r=1:NumReplicas
       inG = mean(model.Likelihood.Genes(:,:,r),1);       
       F(:,:,r) = mean(inG(:))*ones(NumOfTFs,SizF);
       Fu(:,:,r) = repmat(mean(model.Likelihood.Genes(:,1:M,r)),[NumOfTFs 1]); 
    end
    %
else
    %
    F = zeros(NumOfTFs, SizF);
    % initialize the control variables 
    Fu = zeros(NumOfTFs, M);
    inG = mean(model.Likelihood.Genes(:,:,1),1);       
    F = mean(inG(:))*ones(NumOfTFs,SizF);
    Fu = repmat(mean(model.Likelihood.Genes(:,1:M,1)),[NumOfTFs 1]); 
    %
end


% if the initial values of the TFs is kept fixed, then set the first 
% control point to be zero (-5 even GP function is the log of the TF)
FixedFirst = zeros(1,NumOfTFs);
for j=1:NumOfTFs
if model.constraints.Ft0(j)==0  
   FixedFirst(j) = 1;
   if model.constraints.Ft0_value == 0
       if strcmp(model.Likelihood.singleAct,'lin') & strcmp(model.Likelihood.jointAct,'sigmoid') == 0
           if strcmp(model.constraints.replicas,'free')
           Fu(:,1,:) = zeros(NumOfTFs,NumReplicas);
           else
           Fu(:,1) = zeros(NumOfTFs,1);  
           end    
       else
           if strcmp(model.constraints.replicas,'free')
           Fu(:,1,:) = -5*ones(NumOfTFs,NumReplicas); 
           else
           Fu(:,1) = -5*ones(NumOfTFs,1);     
           end 
       end
   end        
end
end

% if the initial condition of the differential equations is zero, then set accordingly the
% corresponding kinetic parameter (initial condition)
for j=1:NumOfGenes
if model.constraints.InitialConds(j) == 0
   %if model.Constraints.A_value == 0
       model.Likelihood.kinetics(:,4) = model.Likelihood.kinetics(:,1)./model.Likelihood.kinetics(:,2);
   %end        
end
end

% Initial input locations of the control points
% (placed in a regular grid)
TimesF = model.Likelihood.TimesF;
step = (max(TimesF) - min(TimesF))/(M-1);
xu = TimesF(1):step:TimesF(end);
xu = xu(1:M);

% optimize these locations 
%[xu f] = minimize(xu(:), 'trace_CondCov', 50, TimesF(:), model.GP{1}.logtheta, FixedFirst(1)); 

% each TF has its own control inputs locations  
n = SizF;
for j=1:NumOfTFs 
    Xu{j} = xu(:)';
end
M = M*ones(1,NumOfTFs); 

%
for j=1:NumOfTFs 
%    
   U = n+1:n+M(j);
   PropDist.qF{j}.m = zeros(n+M(j),1);
   PropDist.qF{j}.K = kernCompute(model.GP{j}, [TimesF(:); Xu{j}(:)]); 
   %PropDist.qF{1}.K = PropDist.qF{1}.K + 0.1*eye(size(PropDist.qF{1}.K)); 
   L=jitterChol(PropDist.qF{j}.K)';
   PropDist.qF{j}.invL = L\eye(n+M(j)); 
   PropDist.qF{j}.LogDetK = 2*sum(log(diag(L)));
      
   % compute the conditional GP prior given the control variables
   [cmuMinus, cSigma, KInvKu] = gaussianFastConditional(PropDist.qF{j}.m', PropDist.qF{j}.K, 1:n, U);
   [L,er]=jitterChol(cSigma);
   if er>0, L = real(sqrtm(cSigma)); end
   %
   if strcmp(model.constraints.replicas,'free')
       %
       for r=1:NumReplicas
       cmu = cmuMinus + Fu(j,:,r)*KInvKu;
       F(j,:,r) = gaussianFastSample(1, cmu, L);
       end
       %
   else
       %
       cmu = cmuMinus + Fu(j,:)*KInvKu;
       F(j,:) = gaussianFastSample(1, cmu, L);
       %
   end
   %
   PropDist.qF{j}.cmuMinus = cmuMinus;
   PropDist.qF{j}.cSigma = cSigma;
   PropDist.qF{j}.KInvKu = KInvKu;
   PropDist.qF{j}.L = L;
   % compute all the conditional variances for the control Points
   for i=1:M(j)
   % 
       G = [1:i-1, i+1:M(j)];  
       [alpha(i), ku(i), KInvK(i,:)] = gaussianFastConditional(PropDist.qF{j}.m(U)', PropDist.qF{j}.K(U,U), i, G);
   %
   end 
   %
   PropDist.qF{j}.alpha = alpha;
   PropDist.qF{j}.ku = ku;
   PropDist.qF{j}.KInvK = KInvK;
   %
   % precompution
   X = [TimesF(:); Xu{j}(:)];
   model.GP{j}.X2 = -2*X*X' + repmat(sum(X.*X,2)',n+M(j),1) + repmat(sum(X.*X,2),1,n+M(j));
   if strcmp(model.constraints.replicas,'free')
       model.Fu{j} = zeros(M(j),NumReplicas);
       for r=1:NumReplicas
          model.Fu{j}(:,r) = Fu(j,:,r)';
       end
   else
       model.Fu{j} = zeros(M(j),1);
       model.Fu{j}(:,1) = Fu(j,:)';
   end
   %
end % end NumOfTFs loop
%
model.F = F; 
model.Xu = Xu;
model.M = M; 

% Initial proposal Gaussian distribution (with diagonal covariance matrices) 
% for the kinetic parameters interection weights and the lengthscale of the GPs 
PropDist.kin = 0.5*ones(NumOfGenes,SizKin);
% interaction weigths and bias 
PropDist.W = 0.5*ones(NumOfGenes,NumOfTFs+1);
PropDist.LengSc = 0.1*(1/model.prior.GPkernel.lenghtScale.sigma2)*ones(1,NumOfTFs);

% additional proposal distribution for the TF kinetic parameters
if isfield(model.Likelihood,'GenesTF')
  PropDist.TFkin = 0.5*ones(NumOfTFs,2);
end

% useful ranges needed in the adaption of the 
% variances of theese proposal distribution 
qKinBelow = 0.000001; qKinAbove = 2;
qWbelow = 0.000001;   qWabove = 2;
qLengScBelow = 0.0001*PropDist.LengSc(1); 
qLengScAbove = 2*PropDist.LengSc(1);
epsilon = 0.1;

cnt = 0;
%
% do the adaption 
while 1
%
%  
   %tic;
   [model PropDist samples accRates] = gpmtfSample(model, PropDist, AdaptOps);
   %toc;
   [model.GP{1}.logtheta(1)]
   accRateF = accRates.F;
   accRateKin = accRates.Kin;
   accRateW = accRates.W;
   accRateLengSc = accRates.LengSc;
   %
   if isfield(model.Likelihood,'GenesTF')
       accRateTFKin = accRates.TFKin;
   end
   
   if AdaptOps.disp == 1
   fprintf(1,'------ ADAPTION STEP #%2d ------ \n',cnt+1); 
   fprintf(1,'Acceptance Rates for GP functions\n');
   for jj=1:NumOfTFs 
       if strcmp(model.constraints.replicas,'free')
          fprintf(1,'TF function #%2d (rows: #%2d replicas, columns: #%2d control points) \n',jj,NumReplicas,M(jj));       
          disp(accRateF{jj}');
       else
          fprintf(1,'TF function #%2d (coupled replicas) : #%2d control points \n',jj,M(jj));       
          disp(accRateF{jj}'); 
       end
   end    
   fprintf(1,'Acceptance Rates for kinetic parameters (per gene))\n');
   disp(accRateKin);
   
   if isfield(model.Likelihood,'GenesTF')
      fprintf(1,'Acceptance Rates for kinetic parameters (per TF-gene))\n');
      disp(accRateTFKin);
   end

   fprintf(1,'Acceptance Rates for Interaction weights (per gene)\n');
   disp(accRateW);
   fprintf(1,'Acceptance Rates for lengthscales (per GP function)\n');
   disp(accRateLengSc);
   fprintf(1,'Average likelihood value %15.8f',mean(samples.LogL));
   if isfield(model.Likelihood,'GenesTF')
       fprintf(1,' TFGenes LogL %15.8f\n',mean(samples.LogLTF));
   else
       fprintf(1,'\n');
   end
   fprintf(1,'------------------------------- \n',cnt+1);
   end   
       
   if (min(accRateKin(:))>15) & (min(accRateW(:))>15) & (min(accRateLengSc(:))>15)
   %
       allTFs = 0;
       for jj=1:NumOfTFs
           if min(accRateF{jj}(:)) > ((0.2/M(jj))*100)
               allTFs =  allTFs + 1;
           end
       end
       %
       if allTFs == NumOfTFs
          disp('END OF ADAPTION: acceptance rates OK');
          %pause
          break;
       end
   end
   
    
   cnt = cnt + 1;
   % do not allow more than 80 iterations when you adapt the proposal distribution
   if cnt == 100
       warning('END OF ADAPTION: acceptance rates were not all OK');
       break;
   end
   
   %%%%%%%%%%%%%%%%%%%%%%% START ADAPTING CONTROL VARIABLES %%%%%%%%%%%%%%%%
   % adapt the proposal distribution over control variables
   % by adding one control variable 
   for j=1:NumOfTFs 
      accRFj = squeeze(accRateF{j}); 
      if (min(accRFj(:)) < ((0.2/M(j))*100)) & (mod(cnt,2)==0) 
          M(j) = M(j)+1; 
          model.M(j) = M(j); 
          step = (max(TimesF) - min(TimesF))/(M(j)-1);
          xu = TimesF(1):step:TimesF(end);
          xu = xu(1:M(j));
          % optimize control input locations
          %[xu f] = minimize(xu(:), 'trace_CondCov', 50, TimesF(:), model.GP{j}.logtheta, FixedFirst(j)); 
          %
          model.Xu{j} = xu(:)'; 
          % initialize the control variables given the current F
          if strcmp(model.constraints.replicas,'free')
              model.Fu{j} = zeros(M(j),NumReplicas);
          else
              model.Fu{j} = zeros(M(j),1);
          end
          %
          U = n+1:n+M(j);
          %
          PropDist.qF{j}.m = zeros(n+M(j),1);
          PropDist.qF{j}.K = kernCompute(model.GP{j}, [TimesF(:); model.Xu{j}(:)]);     
          L=jitterChol(PropDist.qF{j}.K)';
          PropDist.qF{j}.invL = L\eye(n+M(j)); 
          PropDist.qF{j}.LogDetK = 2*sum(log(diag(L)));
      
          [cmuMinus, cSigma, KInvKu] = gaussianFastConditional(PropDist.qF{j}.m', PropDist.qF{j}.K, U, 1:n);
          [L,er]=jitterChol(cSigma);
          if er>0, L = real(sqrtm(cSigma)); end
          %
          if strcmp(model.constraints.replicas,'free')
              for r=1:NumReplicas
                  cmu = cmuMinus + model.F(j,:,r)*KInvKu;
                  sFu = gaussianFastSample(1, cmu, L);
                  model.Fu{j}(:,r) = sFu';
              end
          else
              cmu = cmuMinus + model.F(j,:)*KInvKu;
              sFu = gaussianFastSample(1, cmu, L);
              model.Fu{j}(:,1) = sFu';
          end
          %
          % compute the conditional GP prior given the control variables
          [cmuMinus, cSigma, KInvKu] = gaussianFastConditional(PropDist.qF{j}.m', PropDist.qF{j}.K, 1:n, U);
          [L,er]=jitterChol(cSigma);
          if er>0, L = real(sqrtm(cSigma)); end
          PropDist.qF{j}.cmuMinus = cmuMinus; 
          PropDist.qF{j}.cSigma = cSigma;
          PropDist.qF{j}.KInvKu = KInvKu;
          PropDist.qF{j}.L = L;
          clear alpha ku KInvK;
          for i=1:M(j)
          %  
             G = [1:i-1, i+1:M(j)];  
             [alpha(i), ku(i), KInvK(i,:)] = gaussianFastConditional(PropDist.qF{j}.m(U)', PropDist.qF{j}.K(U,U), i, G);
          %
          end
          PropDist.qF{j}.alpha = alpha;
          PropDist.qF{j}.ku = ku;
          PropDist.qF{j}.KInvK = KInvK; 
          X = [TimesF(:); model.Xu{j}(:)];
          model.GP{j}.X2 = -2*X*X' + repmat(sum(X.*X,2)',n+M(j),1) + repmat(sum(X.*X,2),1,n+M(j));       
        %          
     end
   %
   end
   %%%%%%%%%%%%%%%%%%%%%%% END of ADAPT CONTROL VARIABLES %%%%%%%%%%%%%%%%
   
   %%%%%%%%%%%%%%%%%%%%%%% START of ADAPT KINETICS PROPOSAL %%%%%%%%%%%%%%%%
   % adapt the proposal over the kinetic parameters (desired acceptance rate: 15-35%)
   for j=1:NumOfGenes
      if accRateKin(j) > 35
         % incease the covariance to reduce the acceptance rate
         PropDist.kin(j,:) = PropDist.kin(j,:) + epsilon*PropDist.kin(j,:);
         if PropDist.kin(j,1) > qKinAbove 
             PropDist.kin(j,:) = qKinAbove*ones(1,SizKin);
         end
      end
      if accRateKin(j) < 15
         % decrease the covariance to incease the acceptance rate
         PropDist.kin(j,:) = PropDist.kin(j,:) - epsilon*PropDist.kin(j,:);    
         if PropDist.kin(j,1) < qKinBelow 
             PropDist.kin(j,:) = qKinBelow*ones(1,SizKin);
         end
         %
      end
       %
   end
   %%%%%%%%%%%%%%%%%%%%%%% END of ADAPT KINETICS PROPOSAL %%%%%%%%%%%%%%%%
   
   
   %%%%%%%%%%%%%%%%%%%%%%% START of ADAPT TF-GENES KINETICS PROPOSAL %%%%%%%%%%%%%%%%
   % adapt the proposal over the kinetic parameters (desired acceptance rate: 15-35%) 
   if isfield(model.Likelihood,'GenesTF')
      for j=1:NumOfTFs
      if accRateTFKin(j) > 35
         % incease the covariance to reduce the acceptance rate
         PropDist.TFkin(j,:) = PropDist.TFkin(j,:) + epsilon*PropDist.TFkin(j,:);
         if PropDist.TFkin(j,1) > qKinAbove 
             PropDist.TFkin(j,:) = qKinAbove*ones(1,SizKinTF);
         end
      end
      if accRateTFKin(j) < 15
         % decrease the covariance to incease the acceptance rate
         PropDist.TFkin(j,:) = PropDist.TFkin(j,:) - epsilon*PropDist.TFkin(j,:);    
         if PropDist.TFkin(j,1) < qKinBelow 
             PropDist.TFkin(j,:) = qKinBelow*ones(1,SizKinTF);
         end
         %
      end
       %
      end 
      % 
   end
   %%%%%%%%%%%%%%%%%%%%%%% END of ADAPT TF-GENES KINETICS PROPOSAL %%%%%%%%%%%%%%%%%%
   
   
   %%%%%%%%%%%%%%%%%%%%%%% START of ADAPT WEIGHTS PROPOSAL %%%%%%%%%%%%%%%%
   % adapt the proposal over the interaction weights (desired acceptance rate: 15-35%)
   for j=1:NumOfGenes
      if accRateW(j) > 35
         % incease the covariance to reduce the acceptance rate
         PropDist.W(j,:) = PropDist.W(j,:) + epsilon*PropDist.W(j,:);
         if PropDist.W(j,1) > qWabove 
             PropDist.W(j,:) = qWabove*ones(1,NumOfTFs+1);
         end
      end
      if accRateW(j) < 15
         % decrease the covariance to incease the acceptance rate
         PropDist.W(j,:) = PropDist.W(j,:) - epsilon*PropDist.W(j,:);    
         if PropDist.W(j,1) < qWbelow 
             PropDist.W(j,:) = qWbelow*ones(1,NumOfTFs+1);
         end
         %
      end
       %
   end
   %%%%%%%%%%%%%%%%%%%%%%% END of ADAPT WEIGHTS PROPOSAL %%%%%%%%%%%%%%%%
   
   
   %%%%%%%%%%%%%%%%%%%%%%% START of ADAPT LENGTHSCALES PROPOSAL %%%%%%%%%%%%%%%%
   % adapt the proposal over the interaction weights (desired acceptance rate: 15-35%)
   for j=1:NumOfTFs
      if accRateLengSc(j) > 35
         % incease the covariance to reduce the acceptance rate
         PropDist.LengSc(j) = PropDist.LengSc(j) + epsilon*PropDist.LengSc(j);
         if PropDist.LengSc(j) > qLengScAbove
             PropDist.LengSc(j) = qLengScAbove;
         end
      end
      if accRateLengSc(j) < 15
         % decrease the covariance to incease the acceptance rate
         PropDist.LengSc(j) = PropDist.LengSc(j) - epsilon*PropDist.LengSc(j);    
         if PropDist.LengSc(j) < qLengScBelow 
             PropDist.LengSc(j) = qLengScBelow;
         end
         %
      end
       %
   end
   %%%%%%%%%%%%%%%%%%%%%%% END of ADAPT LENGTHSCALES PROPOSAL %%%%%%%%%%%%%%%%
%
%
end
                                                                                                                                                                                                                                                                                                                                                                                                                                                                              gpmtfComputeGeneODE.m                                                                               0000700 0003466 0000024 00000004452 11312731117 014542  0                                                                                                    ustar   mtitsias                        games                                                                                                                                                                                                                  function PredGenes = gpmtfComputeGeneODE(LikParams, F, R, Gindex)
% Description: Computes the log likelihood of the linear differential
%              equation model with possibly multiple TFs 
%
% Inputs:  
%        -- LikParams: structure with the likelihood parameters             
%        -- F: GP functions for the replica R 
%        -- R: the replica for whihc you compute the log likelihood  
%        -- Gindex: Indicets the gens for which you compute the log
%           likelihood      
%
% Outputs: 
%        -- loglikval: A vector with log likelihood values corresponding 
%           to individual gene contributions 
%        -- PredGenes: Predicted gene expressions in dense grid of times
%           given by TimesF 
%
% Notes: Numerical integration was done using the Trapezoid rule 


Ngenes = size(Gindex,2);

% This useful to speed up computation when TF are fixed (therefore precomputed once)
if ~isempty(LikParams.TF) 
   F = LikParams.TF(:,:,R);
else
   F = gpmtfComputeTF(LikParams, F, 1:LikParams.numTFs);
end
%

% compute the joint activation function of the TFs 
% i.e. g(f_1(u-tau_j),...,f_M(u-tau_j);w_j) 
fx = jointactFunc(LikParams,F,Gindex);
%fx = jointactFunc2(LikParams,F,Gindex);

uu = LikParams.TimesF(LikParams.startTime:end);
Delta = uu(2)-uu(1); 

PredGenes = zeros(Ngenes,size(uu,2));

for m=1:Ngenes
    j = Gindex(m);
    %
    B = LikParams.kinetics(j,1);
    D = LikParams.kinetics(j,2);
    S = LikParams.kinetics(j,3);   
    A = LikParams.kinetics(j,4);
  
    % Trapezoid rule of numerical integration
    %IntVals = exp(D*uu).*fx(m,:);
    ffx = exp(D*uu).*fx(m,LikParams.Tausindex(j):LikParams.Tausindex(j)+LikParams.sizTime-1);
    IntVals = zeros(size(ffx));
    IntVals(2:end) = .5 * Delta*cumsum(ffx(1:end-1) + ffx(2:end));
    %IntVals = Delta*cumtrapz(IntVals);
    %IntVals = IntVals(comInds);
    
    % Simpson rule of integration 
    %ffx = exp(D*uu).*fx(m,:); 
    %IntVals = ffx;
    %IntVals(2:2:end-1) = 4*IntVals(2:2:end-1);
    %IntVals(3:2:end-2) = 2*IntVals(3:2:end-2);
    %IntVals = cumsum(IntVals);
    %IntVals = (Delta/3)*IntVals;
    %IntVals(1:end-1) = IntVals(1:end-1)-(Delta/3)*ffx(1:end-1);
    
    expD = exp(-uu*D);
    PredGenes(m,:) = B/D  + (A - B/D)*expD + S*(expD.*IntVals);
   %
end    
                                                                                                                                                                                                                      gpmtfComputeTF.m                                                                                    0000600 0003466 0000024 00000001520 11267362461 013650  0                                                                                                    ustar   mtitsias                        games                                                                                                                                                                                                                  function PredTFs = gpmtfComputeTF(LikParams, F, TFindex)
% Description: Computes the log likelihood of the linear differential
%              equation model with possibly multiple TFs 
%
% Inputs:  
%        -- LikParams: structure with the likelihood parameters             
%        -- F: GP functions for the replica R 
%        -- R: the replica for whihc you compute the log likelihood  
%        -- TFindex: Indices of the TF-genes for which you compute the log
%           likelihood      
%
% Outputs: 
%        -- PredTFs: Predicted TF functions at time in TimesF 
%



if ~isfield(LikParams,'GenesTF')
    %
    % apply single activation for the GP functions  
    PredTFs = singleactFunc(LikParams.singleAct,F);
    %
else
    %
    % Translational ODE model for the TFs  
    PredTFs = gpmtfComputeTFODE(LikParams, F, TFindex);
   %
end                                                                                                                                                                                gpmtfComputeTFODE.m                                                                                 0000700 0003466 0000024 00000003026 11305234710 014170  0                                                                                                    ustar   mtitsias                        games                                                                                                                                                                                                                  function PredTFs = gpmtfComputeTFODE(LikParams, F, TFindex)
% Description: Computes the log likelihood of the linear differential
%              equation model with possibly multiple TFs 
%
% Inputs:  
%        -- LikParams: structure with the likelihood parameters             
%        -- F: GP functions for the replica R 
%        -- TFindex: Indices of the TF-genes for which you compute the log
%           likelihood      
%
% Outputs: 
%        -- PredTFs: Predicted TF functions at time in TimesF 
%
% Notes: Numerical integration was done using the Trapezoid rule


Ntfs = size(TFindex,2);

% apply single activation for the GP functions  
fx = singleactFunc(LikParams.singleAct,F);

%
uu = LikParams.TimesF;
Delta = uu(2)-uu(1); 

PredTFs = zeros(size(fx));

for m=1:Ntfs
    j = TFindex(m);
    % 
    D = LikParams.kineticsTF(j,1);
    S = LikParams.kineticsTF(j,2);    
    
    % Trapezoid rule of numerical integration
    ffx = exp(D*uu).*fx(m,:);
    IntVals = zeros(size(ffx));
    IntVals(2:end) = .5 * Delta*cumsum(ffx(1:end-1) + ffx(2:end));
    %IntVals = Delta*cumtrapz(IntVals);
    
    
    % Simpson rule of integration 
    %ffx = exp(D*uu).*fx(m,:); 
    %IntVals = ffx;
    %IntVals(2:2:end-1) = 4*IntVals(2:2:end-1);
    %IntVals(3:2:end-2) = 2*IntVals(3:2:end-2);
    %IntVals = cumsum(IntVals);
    %IntVals = (Delta/3)*IntVals;
    %IntVals(1:end-1) = IntVals(1:end-1)-(Delta/3)*ffx(1:end-1);
    
    expD = exp(-uu*D);
    PredTFs(m,:) = S*(expD.*IntVals);
    %
end
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           gpmtfCreate.m                                                                                       0000700 0003466 0000024 00000034435 11320453720 013206  0                                                                                                    ustar   mtitsias                        games                                                                                                                                                                                                                  function model = gpmtfCreate(Genes, GenesVar, GenesTF, GenesTFVar, TimesG, TimesF, options)
%model = gpmtfCreate(Genes, GenesVar, GenesTF, GenesTFVar, TimesG, TimesF, options)
%
% Description: Creates the structure variable that contains all data, the
%              parameters, the priors over those parameters, the form 
%              of the likelihood etc
%
% Inputs: 
%     -- Genes:  NumOfGenes x NumOfTimes x Replicas matrix that stores the 
%            gene expressions for all genes in all times and replicas
%     -- GenesVar:  NumOfGenes x NumOfTimes x Replicas matrix
%            that specifies the variances of the gene expressions 
%            in the Gaussian likelihood (a separate variance is given 
%            for each element in the Genes). This is useful when 
%            GeneVars have been estimated previously from a package like 
%            PUMA (http://www.bioinf.manchester.ac.uk/resources/puma/).  
%            When this input argument is not given, then GeneVars is 
%            considered a random variable that needs to be inferred by MCMC. 
%            In such case a single variance is estimated for each gene, 
%            i.e. GeneVars(i,:,:) = sigma2_i * ones(NumOfTimes,Replicas)
%            and we sample sigma2_i
%     -- GenesTF:  NumOfTFs x NumOfTimes x Replicas matrix that stores the 
%            gene expressions for the TF genes (if given, otherwize is empty)
%     -- GenesTFVar:  NumOfTFs x NumOfTimes x Replicas matrix
%            that specifies the variances of the TF gene expressions similarly 
%            to GenesVar (if given, otherwize is empty).
%     -- TimesG: The time points where gene expressions are evaluated 
%     -- TimesF: The times where the GP functions (TFs) are evaluated.
%            Note that size(TimesF)>>size(TimesG).
%     -- options: Determines the functional form of the activation 
%            of the TFs in the differential equation. It contains 3 fields   
%            * singleAct = 'lin' or 'exp' or 'logOnePlusExp'. This specifies the
%                  individual activation of the GP functions   
%            * jointAct = 'genHill' or 'sigmoid' or 'michMenten' 
%                          (see the paper for definition)
%            * jointActBin = 0 or 1. Binarize the outputs of the joint
%                  activation function. If 0 the outputs are not binarized.
%                  If 1 the outputs are binarized.
%                  !IMPORTANT NOTE! singleAct = 'exp' should be used only 
%                  with the jointAct='michMentenAct' or 'michMentenRepres'
%                  but not for TFjointAct = 'sigmoid' (SEE DEMOS)  
%            * Constraints: Determines several constraints of the parameter during
%                   sampling. It contains 3 fields
%               - Ft0: 1 x NumOfTFs  (number of TFs; see below) binary vector that 
%                  detetmines is the initial value of the TF is 0 or not: 
%                  if Ft0(j)=0, then the j TF has concentration 0 at time t=0
%                  if Ft0(j)=1, then the j TF is free to take any value
%               - InitialConds: 1 x NumOfGenes binary vector that decides if the
%                  initial conditions of each differential equation for each gene  
%                  are constrained to be zero at the initial time (say t=0)  
%                  if InitialConds(j)=0, then the ODE initial cond for the jth 
%                  gene is 0. Otheriwise if is free to take any value at time t=0
%                - W:  NumOfGenes x NumOfTFs binary matrix that determines constraints 
%                  coming from side information about which TFs do not regulate 
%                  certain genes. This means that certain values in the interaction 
%                  matrix W are constrained to be zero. if w(i,j)=1 (for i gene and j TF) 
%                  then the interaction is allowed and the weight w(i,j) is learned
%                  If w(i,j)=0, then no interaction is allowed and the weights w(i,j) is 
%                  constrained to be zero. 
% 
%                      
% Outputs: 
%     -- model: The model stucture and parameters. In consists of 3 fields 
%            * Likelihood: The likelihood structure (functional form of the 
%                  activation of the TFs; see options input argument) 
%                  and parameters, e.g. kinetics and observation noise 
%            * GP: The GP prior type and kernel parameters. Only the 
%                  Gaussian or 'rbf' kernel function with unit variance 
%                  and unknown lengthscale 'ell^2' is currently supported:
%
%                    cov = exp[-0.5 * (x_i - x_j)^2/ell^2]
%
%                  The lengthscale ell^2 is encoded as 
%                  logth = 0.5*log(ell^2) = log(ell) and a Gaussian prior 
%                  is placed on the logth variable   
%
%            * prior: the prior distribution over all parameters 
%            * constraints: Specifies the status of several variables
%                  that can be constrained ('free') or 'fixed' during 
%                  sampling 
%      

% create the Likelihood model structrure
% differential equation model
model.Likelihood.type = 'diffeq';
% number of genes
[NumOfGenes NumOfTimes NumOfReplicas] = size(Genes);

% number of genes 
model.Likelihood.numGenes = NumOfGenes;
% number of times the gene expression is evaluated
model.Likelihood.numTimes = NumOfTimes;
% number of replicas per gene 
model.Likelihood.numReplicas = NumOfReplicas;
% number of transcription factors
model.Likelihood.numTFs = options.numTFs;

model.Likelihood.Genes = Genes;
% gene variances
if ~isempty(GenesVar)
   model.Likelihood.sigmas = GenesVar;
   model.constraints.sigmas = 'fixed';  % 'free' or 'fixed'
else
   model.Likelihood.sigmas = 0.05*ones(NumOfGenes, NumOfTimes, NumOfReplicas);
   model.constraints.sigmas = 'free';  % 'free' or 'fixed'
end


if ~isempty(GenesTF)
   % 
   model.Likelihood.GenesTF = GenesTF;
   if ~isempty(GenesTFVar)
      model.Likelihood.sigmasTF = GenesTFVar; 
      model.constraints.sigmasTF = 'fixed';  % 'free' or 'fixed'
   else
      model.Likelihood.sigmasTF = 0.05*ones(options.numTFs, NumOfTimes, NumOfReplicas);
      model.constraints.sigmasTF = 'free';  % 'free' or 'fixed'
   end   
   %
end


model.Likelihood.TimesG = TimesG;
model.Likelihood.TimesF = TimesF;
model.Likelihood.step = TimesF(2)-TimesF(1);

model.Likelihood.kineticsReadme = '1st column: basals, 2st decays, 3st sensitivities, 4st initial conds';
% decay rates
model.Likelihood.kinetics(:,1) = rand(NumOfGenes,1) + 0.5; %B: basal rates 
% sensitivity parameters
model.Likelihood.kinetics(:,2) = rand(NumOfGenes,1) + 0.5; %D: decay rates
% basal rates
model.Likelihood.kinetics(:,3) = rand(NumOfGenes,1) + 0.5; %S: sensitivities 
% initial conditions
model.Likelihood.kinetics(:,4) = model.Likelihood.kinetics(:,1)./model.Likelihood.kinetics(:,2);  % A: initial conditions
% Delays 
model.Likelihood.Tausreadme = 'taus are the delays in the ODEs'; 
model.Likelihood.Taus = zeros(1,NumOfGenes); 
model.Likelihood.Tausindex = (size(TimesF,2) - options.sizTime + 1)*ones(1,NumOfGenes); 
model.Likelihood.sizTime = options.sizTime;
model.Likelihood.startTime = (size(TimesF,2) - options.sizTime + 1);
model.Likelihood.tauMax = options.tauMax;

uu = model.Likelihood.TimesF(model.Likelihood.startTime:end);
[commonSlots, comInds] = intersect(uu,model.Likelihood.TimesG);
model.Likelihood.comInds = comInds; 

uu = model.Likelihood.TimesF;
[commonSlots, comIndsTF] = intersect(uu,model.Likelihood.TimesG);
model.Likelihood.comIndsTF = comIndsTF;

%model.Likelihood.numKins_perGene = 4;

%  additional kinetic parameters for the TF Genes
if ~isempty(GenesTF)
   %
   model.Likelihood.kineticsTFReadme = '1st column: decays, 2st sensitivities';
   model.Likelihood.kineticsTF(:,1) = 0.5*ones(options.numTFs,1);%rand(options.numTFs,1) + 0.5; %d: decay rates
   model.Likelihood.kineticsTF(:,2) = ones(options.numTFs,1);%rand(options.numTFs,1) + 0.5; %s: sensitivities
   %
end

% Gene - TFs interaction weights 
% (relevant only for multiple transcription factors)
model.Likelihood.W = 0.05*randn(NumOfGenes,options.numTFs);
% gene-specific bias terms for the Gene - TFs interaction 
% functions
model.Likelihood.W0 = 0.05*randn(NumOfGenes,1);


% single activation
model.Likelihood.singleAct = options.singleAct;
% valid values of jointAct are : 'sigmoid', 'michMenten'
model.Likelihood.jointAct = options.jointAct;
% binarize or not the joint activation function 
model.Likelihood.jointActBin = options.jointActBin;


% *** This feature is only used when 'michMenten' joint activation is
%     considered. When  the 'sigmoid' joint activation is used 
%     this feature remains inactive        
% ***
% -1: repression, 0 no regulation, 1 for activation 
% (default values are to assume that all TFs are activators) 
model.Likelihood.Net_Learn = 'no';
if strcmp(model.Likelihood.jointAct,'michMenten');
model.Likelihood.Net_X = options.Net_X;  
%model.Likelihood.Net_Learn = options.Net_Learn;
end

% this may be used to stored precomputed TF profiles 
% (for saving computations purposes)
model.Likelihood.TF = [];

% create the GP prior model 
timescale = 1.5*(max(TimesG(:))-min(TimesG(:)))/(size(TimesG,2)-1);
lengthscale = (max(TimesG(:))-min(TimesG(:)))/10;%(size(TimesG,2)-1);
%lengthscale = 1.5*(max(TimesG(:))-min(TimesG(:)))/(size(TimesG,2)-1);
for j=1:options.numTFs
   model.GP{j}.type = {'rbf','white'};
   model.GP{j}.TF = j;
   model.GP{j}.logtheta = [log(lengthscale) 0 0.5*log(1e-06)];
end
X = TimesF(:);
[n D] = size(X);
model.GP{1}.X2 = -2*X*X' + repmat(sum(X.*X,2)',n,1) + repmat(sum(X.*X,2),1,n);


% define the prior of the TFs regulation types 
% (repression, non-regulation, activation)
                       % repression, non-regulation, activation
model.prior.Net.type = 'discrete';
model.prior.Net.contraint = 'probability';
model.prior.Net.priorSpace = 'lin'; % it means NoTransform;
model.prior.Net.prob = [0.025 0.95 0.025];
model.prior.Net.readme = 'The 3 prior probabilities correspond to: repression, non-regulation, activation';


% prior for the kinetics parameters
model.prior.kinetics.type = 'normal';
%model.prior.kinetics.readme = 'kinetics are: decays, sensitivities, basals, initial conds'
model.prior.kinetics.contraint = 'positive';
model.prior.kinetics.priorSpace = 'log';
model.prior.kinetics.mu = -0.5; % mean 
model.prior.kinetics.sigma2 = 2; % variance

model.prior.delays.type = 'discrete'; 
a = 1; 
AllTaus = model.Likelihood.tauMax:model.Likelihood.step:-eps;
AllTaus = [AllTaus, 0]; 
model.prior.delays.prob = exp(a*AllTaus)/(sum(exp(a*AllTaus)));

% prior for the interaction bias
model.prior.weight0.type = 'normal'; %'normal' or Laplace
if strcmp(model.Likelihood.jointAct,'michMenten')
model.prior.weight0.constraint = 'positive';
model.prior.weight0.priorSpace = 'log';
model.Likelihood.W0 = rand(NumOfGenes,1)+0.1;
else
model.prior.weight0.constraint = 'real';
model.prior.weight0.priorSpace = 'lin'; % it means NoTransform;   
end
model.prior.weight0.mu = 0;
model.prior.weight0.sigma2 = 1.5;
    
% prior for the interactino weigths (e.g. inside the sigmoid)
if strcmp(options.spikePriorW,'no')
    %
    model.prior.weights.type = 'normal'; %'normal' or Laplace
    if strcmp(model.Likelihood.jointAct,'michMenten') | strcmp(options.constraints.spaceW,'positive'); 
    model.prior.weights.constraint = 'positive';
    model.prior.weights.type = 'truncNormal';
    model.prior.weights.priorSpace = 'lin';
    model.Likelihood.W = rand(NumOfGenes,options.numTFs)+0.1;
    else
    model.prior.weights.constraint = 'real';
    model.prior.weights.priorSpace = 'lin'; % it means NoTransform;   
    end
    model.prior.weights.mu = 0;
    model.prior.weights.sigma2 = 1.5;
    %
else
    %
    % use the two mixture model with the spike 
    model.prior.weights.type = 'spikeNormal';
    if strcmp(model.Likelihood.jointAct,'michMenten') | strcmp(options.constraints.spaceW,'positive'); 
    model.prior.weights.constraint = 'positive';
    model.prior.weights.type = 'spikeTruncNormal';
    model.prior.weights.priorSpace = 'lin';
    model.Likelihood.W = rand(NumOfGenes,options.numTFs)+0.1;
    else    
    model.prior.weights.constraint = 'real';
    model.prior.weights.priorSpace = 'lin';
    end
    model.prior.weights.mu = 0;
    model.prior.weights.sigma2 = 1.5;
    model.prior.weights.spikeMu = 0; 
    model.prior.weights.spikeSigma2 = 0.0001;
    model.prior.weights.pis = options.spikepriors; 
    % binary variables (randomly chosen for each TF and gene)
    model.prior.weights.S = round(rand(NumOfGenes,options.numTFs));
    %
end

% prior for the gene specifc noise variances (the inverce of them) 
model.prior.invsigma2.type = 'gamma';
model.prior.invsigma2.constraint = 'positive';
model.prior.invsigma2.priorSpace = 'lin'; % it means NoTransform;
model.prior.invsigma2.a = 0.1;
model.prior.invsigma2.b = 0.01;

% prior for the lengthscales
model.prior.GPkernel.lenghtScale.type = 'normal';
model.prior.GPkernel.lenghtScale.constraint = 'positive';
model.prior.GPkernel.lenghtScale.priorSpace = 'log';
ok = 1.5*(max(TimesG(:))-min(TimesG(:)))/8;%(size(TimesG,2)-1);
model.prior.GPkernel.lenghtScale.mu = 2*log(ok); % mean 
model.prior.GPkernel.lenghtScale.sigma2 = 2;         % variance
%model.prior.GPkernel.lenghtScale.a = 1;
%model.prior.GPkernel.lenghtScale.b = 1/(ok^2);
 

% initial value of the TFs
% if 0, then the TF has concentration 0 at time t=0
% if 1, then the TF is free to take any value at time t=0
model.constraints.Ft0 = options.constraints.Ft0;
% initial value of the TF (this feature 
% this feature is activated only if model.Constraints.Ft0 = 'fixed')
model.constraints.Ft0_value = 0; 

% initial condition of the differential equation 
% if 0, then the ODE initial cond for the jth gene is 0 at time t=0
% if 1, then the ODE initial cond is free to take any value at time t=0
model.constraints.InitialConds = options.constraints.initialConds; 
% initial value of the TF (this feature 
% this feature is activated only if model.Constraints.InitialConds = 0)
model.constraints.InitialConds_value = 0; 

%
model.constraints.geneTFsensitivity = options.constraints.geneTFsensitivity; 
model.constraints.geneTFsensitivity_value = 1; 

% constraints on the interaction weigths between TF and genes 
model.constraints.W = options.constraints.X;


% constraints on the interaction weigths between TF and genes 
model.constraints.W0 = options.constraints.W0;
% constraint on TF functions across different replicas  
model.constraints.replicas = options.constraints.replicas;
                                                                                                                                                                                                                                   gpmtfDiscretize.m                                                                                   0000700 0003466 0000024 00000002175 11267362461 014117  0                                                                                                    ustar   mtitsias                        games                                                                                                                                                                                                                  function [options, TimesF] = gpmtfDiscretize(TimesG, options)
% DISCRETIZE the GP function
% (around discr times more the number the discrete time points
% we have gene expressions)


% TimesF discretizes the [-tau_max, T] where T = max(TimesG)
% - First discretize in [0,T] (where we have observed genes) 

% the GPs functions is going to be descretize in  discr*(size(TimesG,2)-1)+1
% points. This number of points must be a odd number 
discr=10;

if (discr*(size(TimesG,2)-1))+1 > 200 
   discr = floor(200/(size(TimesG,2)-1)) + 1; 
end 
if mod(discr*(size(TimesG,2)-1)+1,2) == 0
   discr = discr-1;
end 
step = ((max(TimesG) - min(TimesG))/(size(TimesG(:),1)-1))/discr;
TimesF =[]; TimesF(1) = TimesG(1);
for j=1:size(TimesG(:),1)-1
   TimesF = [TimesF, ((TimesG(j)+step):step:TimesG(j+1))];
   if TimesF(end) ~= TimesG(j+1)
      TimesF = [TimesF, TimesG(j+1)];
   end
end

% - Now discretize in [-tau_max,0) (the "delay" part of the TF) 
options.sizTime = size(TimesF,2);
if options.tauMax > 0
DelayP = -step:-step:-options.tauMax; 
options.tauMax = DelayP(end); 
TimesF = [DelayP(end:-1:1) TimesF];    
end
                                                                                                                                                                                                                                                                                                                                                                                                   gpmtfGenerateToydata.m                                                                              0000700 0003466 0000024 00000010335 11306031463 015054  0                                                                                                    ustar   mtitsias                        games                                                                                                                                                                                                                  function [LikParams lengthscale] = gpmtfGenerateToydata(NumGs, TimesG, NumTFs, TimesF)
%
%


% generate NumTFs different GP functions using the exponential kernel

n = size(TimesF(:),1);
SizG = size(TimesG(:),1)
sigma2 = 0.01; 
sigmaf = 1;
lengthscale = [2 1.5 0.5];
contr = [n+1 n+2 n+3 n+4 n+5 n+6 n+7 n+8]; 
FF = [-5 -2 0.5 -1 -2 -3 -3 -4;
      -5 -4 -4 -3 -3 -1 0 -3;
      -5 -4 -4 -1 0 -4 -4 -5];
step = (TimesF(end) - TimesF(1))/(size(contr,2)-1);  
Xu = TimesF(1):step:TimesF(end);
step
Xu
pause
%Xu = TimesG(1);
%FF = (-5)*ones(NumTFs,1);
%contr = size(TimesF(:),1)+1;
FF = FF(1:NumTFs,:);
%lengthscale = 4*rand(NumTFs,1) + 0.5;

un = size(Xu,2);

m = zeros(n+un,1);
for i =1:NumTFs
    covfunParams.logtheta = [0.5*log(lengthscale(i)) 0.5*log(sigmaf) 0.5*log(sigma2)];
    K = kernCompute(covfunParams, [TimesF(:); Xu(:)], [TimesF(:); Xu(:)]);
   
    [cmuMinus, cSigma, KInvKu] = gaussianFastConditional(m', K, 1:n,contr);
    [L,er]=jitterChol(cSigma);
    if er>0, L = real(sqrtm(cSigma)); end
    cmu = cmuMinus + FF(i,:)*KInvKu;
    F(i,:) = gaussianFastSample(1, cmu, L);
    if i ==1
    plot(TimesF,exp(F(i,:)));
    hold on;
    elseif i==2
    plot(TimesF,exp(F(i,:)),'r');
    elseif i==3
    plot(TimesF,exp(F(i,:)),'g');
    else
    plot(TimesF,exp(F(i,:)),'c');    
    end
end

max(F(:))
pause

% produce Random interaction Weights (the last column are the bias terms) 
W = randn(NumGs,NumTFs);
ok = zeros(NumGs,NumTFs);
for j=1:NumGs  
  ch = randperm(NumTFs);
  ok(j,ch(1)) = 1;
  if (NumTFs > 1) & (rand > 0.75)
      ok(j,ch(2)) = 1;
  end
end
W = W.*ok;

Net_X = (-1)*(W<0) + 1*(W>0); 
Net_X = ok;    
  
W0 = randn(NumGs,1);

% produce random kinetic parameters for each gene 
for j=1:NumGs
    %B D S A 
    %ok = randperm(4);% randn;
    %Kinetics(j,:) = rand(1,4).*ok; 
    Kinetics(j,:) = 2*rand(1,4)-1;
    Kinetics(j,:) = exp(Kinetics(j,:));%exp(rand(1,4));
    % decays fixed to one
    %Kinetics(j,2) = 1; 
    Kinetics(j,4) = Kinetics(j,1)/Kinetics(j,2);  
end


% produce random kinetic parameters for each TF-gene 
for j=1:NumTFs
    %D S A 
    KineticsTF(j,:) = 2*rand(1,2)-1;
    KineticsTF(j,:) = exp(KineticsTF(j,:));
    %
end
KineticsTF(:,2) = ones(NumTFs,1);

% solve numerically the differential equation
LikParams.numTFs = NumTFs;
LikParams.F = F;
LikParams.Genes = zeros(NumGs,SizG);
LikParams.GenesTF = zeros(NumTFs,SizG);
LikParams.TimesG = TimesG; 
LikParams.TimesF = TimesF;
LikParams.kinetics = Kinetics;
LikParams.kineticsTF = KineticsTF;
LikParams.Net_X = Net_X;
LikParams.W = W;
LikParams.W0 = W0;
LikParams.sigmas = sigma2*ones(NumGs,SizG);
LikParams.sigmasTF = sigma2*ones(NumTFs,SizG);
LikParams.singleAct = 'logOnePlusExp';
LikParams.jointAct = 'genHill';
LikParams.jointActBin = 0;
LikParams.startTime = 1;

step = TimesF(2)-TimesF(1);
NofTaus = floor(abs(TimesF(1))/step);
NofTaus
ok = 0:step:TimesG(end);
if ok(end) < TimesG(end)
    ok = [ok, TimesG(end)];
end

LikParams.sizTime = size(ok,2);
sizTime = LikParams.sizTime; 
LikParams.startTime = size(TimesF,2) - sizTime + 1;
LikParams.startTime
TimesF(LikParams.startTime)
LikParams.Tausreadme = 'taus are the delays in the ODEs'; 

for j=1:NumGs
  LikParams.Taus(j) = 0; 
  ok = randperm(NofTaus)-1;
  ok = ok(1); 
  size(TimesF,2) - sizTime + 1 - ok;
  LikParams.Tausindex(j) = size(TimesF,2) - sizTime + 1 - ok;
  LikParams.Taus(j) = ok*step; 
end 
LikParams.Taus
LikParams.Tausindex

uu = LikParams.TimesF(LikParams.startTime:end);
[commonSlots, comInds] = intersect(uu,LikParams.TimesG);
LikParams.comInds = comInds; 

[loglikval, Genes] = gpmtfLogLikelihoodgene(LikParams, F, 1, 1:NumGs);
[commonSlots, comInds] = intersect(TimesF,TimesG);
uu = TimesF(LikParams.startTime:end);
[commonSlots, comInds] = intersect(uu,TimesG);
Genes = Genes(:,comInds);
% add random noise 
Genes = Genes + sqrt(sigma2)*randn(size(Genes));

uu = LikParams.TimesF;
[commonSlots, comIndsTF] = intersect(uu,LikParams.TimesG);
LikParams.comIndsTF = comIndsTF;
[loglikval, GenesTF] = gpmtfLogLikelihoodgeneTF(LikParams, F, 1, 1:NumTFs);

% add random noise 
GenesTF = GenesTF + sqrt(sigma2)*randn(size(GenesTF));
uu = TimesF;
[commonSlots, comInds] = intersect(uu,TimesG);
GenesTF = GenesTF(:,comInds);

LikParams.TF = gpmtfComputeTF(LikParams, F,  1:NumTFs);
LikParams.Genes = Genes; 
LikParams.GenesTF = GenesTF;                                                                                                                                                                                                                                                                                                   gpmtfLogLikelihoodGene.m                                                                            0000700 0003466 0000024 00000004015 11312733016 015316  0                                                                                                    ustar   mtitsias                        games                                                                                                                                                                                                                  function [loglikval, PredGenes] = gpmtfLogLikelihoodGene(LikParams, F, R, Gindex)
%function [loglikval, PredGenes] = gpmtfLogLikelihoodgene(LikParams, F, R, Gindex)
%
% Description: Computes the log likelihood of the linear differential
%              equation model with multiple TFs 
%
% Inputs:  
%        -- LikParams: structure with the likelihood parameters             
%        -- F: GP functions for the replica R (from whihc TFs are constructed) 
%        -- R: the replica for which you compute the log likelihood  
%        -- Gindex: Indices of the genes for which you compute the log
%           likelihood      
%
% Outputs: 
%        -- loglikval: A vector with log likelihood values corresponding 
%           to individual gene contributions 
%        -- PredictedGenes: Predicted gene expressions in dense grid of times
%           given by TimesF 


%Genes = LikParams.Genes(:,:,R);
PredGenes = gpmtfComputeGeneODE(LikParams, F, R, Gindex);

%uu = LikParams.TimesF(LikParams.startTime:end);
%[commonSlots, comInds] = intersect(uu,LikParams.TimesG);


%
%for m=1:size(Gindex,2)
%    %
%    j = Gindex(m);
%    loglikval(m) = - 0.5*sum(log(2*pi*LikParams.sigmas(j,:,R)),2)....
%                   - 0.5*sum(((Genes(j,:) - PredGenes(m,LikParams.comInds)).^2)./LikParams.sigmas(j,:,R),2);
%               
%    %totalmse(m) = sum(sum((Genes(j,:) -  PredGenes(m,LikParams.comInds)).^2)); 
%    
%    %[Genes(j,:);PredGenes(m,comInds)]
%    %plot(LikParams.TimesG, LikParams.Genes(j,:),'r');
%    %hold on; 
%    %plot(LikParams.TimesF(31:end), PredGenes(m,:));
%    %pause
%    %hold off
%    %
%end


%loglikval = - 0.5*sum(log(2*pi*LikParams.sigmas(Gindex,:,R)),2)....
%                   - 0.5*sum(((Genes(Gindex,:) - PredGenes(:,LikParams.comInds)).^2)./LikParams.sigmas(Gindex,:,R),2);
loglikval = - 0.5*sum(log(2*pi*LikParams.sigmas(Gindex,:,R)),2)....
                   - 0.5*sum(((LikParams.Genes(Gindex,:,R) - PredGenes(:,LikParams.comInds)).^2)./LikParams.sigmas(Gindex,:,R),2);
             
loglikval = loglikval';
                 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   gpmtfLogLikelihoodGeneTF.m                                                                          0000700 0003466 0000024 00000004073 11306027426 015560  0                                                                                                    ustar   mtitsias                        games                                                                                                                                                                                                                  function [loglikvalTF, PredGenesTF] = gpmtfLogLikelihoodGeneTF(LikParams, F, R, TFindex)
% function [loglikval, PredictedGenes] = logLTFdiffEquation(LikParams, F, R, Gindex, Genes, TimesG, TimesF)
%
% Description: Computes the log likelihood of the linear differential
%              equation model with possibly multiple TFs 
%
% Inputs:  
%        -- LikParams: structure with the likelihood parameters             
%        -- F: GP functions for the replica R 
%        -- R: the replica for whihc you compute the log likelihood  
%        -- Gindex: Indicets the gens for which you compute the log
%           likelihood      
%        -- Genes: Expression of all gene for the replica R
%        -- TimesG: Times you have gene expression measurements 
%        -- TimesF: Times you you have discretize the TFs 
%
% Outputs: 
%        -- loglikval: A vector with log likelihood values corresponding 
%           to individual gene contributions 
%        -- PredictedGenes: Predicted gene expressions in dense grid of times
%           given by TimesF 
%


GenesTF = LikParams.GenesTF(:,:,R);
PredGenesTF = singleactFunc(LikParams.singleAct, F(TFindex,:));

%uu = LikParams.TimesF;
%[commonSlots, comInds] = intersect(uu,LikParams.TimesG);

%for m=1:size(TFindex,2)
%    %
%    j = TFindex(m);
%    loglikvalTF(m) = - 0.5*sum(log(2*pi*LikParams.sigmasTF(j,:,R)),2)....
%                   - 0.5*sum(((GenesTF(j,:) - PredGenesTF(m,LikParams.comInds)).^2)./LikParams.sigmasTF(j,:,R),2);
%    
%                   
%    %totalmse(m) = sum(sum((GenesTF(j,:) -  PredGenesTF(m,LikParams.comInds)).^2)); 
%    
%    %[GenesTF(j,:);PredGenesTF(m,comInds)]
%    %plot(LikParams.TimesG, LikParams.GenesTF(j,:),'r');
%    %hold on; 
%    %plot(LikParams.TimesF, PredGenesTF(m,:));
%    %disp('in TF')
%    %pause
%    %hold off
%    %
%end


loglikvalTF = - 0.5*sum(log(2*pi*LikParams.sigmasTF(TFindex,:,R)),2)....
               - 0.5*sum(((GenesTF(TFindex,:) - PredGenesTF(:,LikParams.comIndsTF)).^2)./LikParams.sigmasTF(TFindex,:,R),2);

loglikvalTF = loglikvalTF';
%sum(abs(loglikvalTF - loglikvalTF1))
%pause
                                                                                                                                                                                                                                                                                                                                                                                                                                                                     gpmtfOptions.m                                                                                      0000700 0003466 0000024 00000011232 11313747237 013437  0                                                                                                    ustar   mtitsias                        games                                                                                                                                                                                                                  function options = gpmtfOptions(Genes,numTFs) 
%
%
%


NumOfGenes = size(Genes,1);

if numTFs == 1
        %
        % NUMBER of TFs
        options.numTFs = 1;
        %
        % INDIVIDUAL transform of the GP function
        % --- either  'logOnePlusExp' or 'loglogOnePlusExp' or 'lin' or 'exp'
        % --- 'exp' or 'logOnePlusExp' must combine with TFjointAct= 'lin' or 'michMenten'
        % --- 'loglogOnePlusExp' or 'lin' with TFjointAct='sigmoid'
        options.singleAct = 'logOnePlusExp';
        %
        % - JOINT transform by passing through a sigmoid type of function.
        options.jointAct = 'michMenten';
        
        options.Net_X = ones(NumOfGenes,1);
        %        
        % Make the output BINARY at the end or not 
        options.jointActBin = 0; % binarize the final output of the joint activation function 
                                 % to 0-1 (where the threshold 0.5)
                                 % if 0 the outputs are not binarilzed
                                 % if 1 the outputs are binarized
        %                         
        % maximum positve value for the gene DELAYS
        % -- 0 means that no delays are allowed 
        options.tauMax = 0;
        %
        % CONSTRAINTS of the initial value of the TF at time t=0. 
        % if 0, then the TF has concentration 0 at time t=0
        % if 1, then the TF is free to take any value at time t=0
        options.constraints.Ft0 = 0; 
        %
        % CONSTRAINTS of the initial conditions of each differential equation 
        % for each gene  
        % if 0, then the ODE initial cond for the jth gene is 0 at time t=0
        % if 1, then the ODE initial cond is free to take any value at time t=0
        options.constraints.initialConds = ones(1,NumOfGenes); 
        options.constraints.geneTFsensitivity = 0;
        %
        %
        options.constraints.X = ones(NumOfGenes,1);
        options.constraints.replicas = 'free'; % 'free' or 'coupled'
        %
else
        %
        % NUMBER of TFs
        options.numTFs = numTFs;
        %
        % Define the activation function that transforms the GP functions. 
        % For multiple TFs the GP functions are linearly mixed and then are passed 
        % through a sigmoid function
        % INDIVIDUAL transform of the GP function
        % --- either  'loglogOnePlusExp' or 'lin' or 'exp'
        % --- 'exp' must combined with TFjointAct= 'lin'
        % --- the rest with TFjointAct='sigmoid'
        options.singleAct = 'logOnePlusExp';
        %
        % - JOINT transform by passing through a sigmoid type of function.
        options.jointAct = 'genHill';
        
        % Make the output BINARY at the end or not 
        options.jointActBin = 0; % binarize the final output of the joint activation function 
                                 % to 0-1 (where the threshold 0.5)
                                 % if 0 the outputs are not binarilzed
                                 % if 1 the outputs are binarized
        %                         
        % maximum positve value for the gene DELAYS
        % -- 0 means that no delays are allowed 
        options.tauMax = 0;
        
        % use or not the two mixture (one spike and one borad) for the interaction weights 
        options.spikePriorW = 'no'; 
        %
        % CONSTRAINTS of the initial value of the TF at time t=0. 
        % if 0, then the TF has concentration 0 at time t=0
        % if 1, then the TF is free to take any value at time t=0
        options.constraints.Ft0 = zeros(1,options.numTFs); 
        %
        % CONSTRAINTS of the initial conditions of each differential equation 
        % for each gene  
        % if 0, then the ODE initial cond for the jth gene is 0 at time t=0
        % if 1, then the ODE initial cond is free to take any value at time t=0
        options.constraints.initialConds = ones(1,NumOfGenes); 
        %
        options.constraints.geneTFsensitivity = zeros(1,numTFs);
        
        % constrain interaction bias to be zero
        options.constraints.W0 = ones(1,NumOfGenes); 
        
        %
        options.constraints.spaceW = 'real'; 
        %
        % CONSTRAINTS coming from side information about which TFs do not regulate 
        % certain genes. This means that certain values in the interaction 
        % matrix W are constrained to be zero. 
        % if w(i,j)=1 (for i gene and j TF) then the interaction is 
        % allowed and the weight w(i,j) is learned
        % if w(i,j)=0, then no interaction is allowed and the weights w(i,j) is 
        % constrained to be zero. 
        options.constraints.X = ones(NumOfGenes,numTFs);
        options.constraints.replicas = 'free'; % 'free' or 'coupled'
        %
end
                                                                                                                                                                                                                                                                                                                                                                      gpmtfPlot.m                                                                                         0000700 0003466 0000024 00000035657 11305744544 012742  0                                                                                                    ustar   mtitsias                        games                                                                                                                                                                                                                  function gpmtfPlot(model, samples, demdata, printResults)
%function gpmtfPlot(model, samples, demdata, printResults)
%
% Description: Creates plots to display the outcome of MCMC 
%
% Inputs: 
%      -- model: Contains the structure of the GP model 
%      -- samples: A structure that conatsint he samples
%      -- Genes: Expression of the genes 
%      -- TimesG: Times when expression measuremtns are available 
%      -- TimesF: Time discretization of the TF (TimesF >> TimesG) 
%      -- demdata: a string that characterizes the experiments, e.g. 'p53' 
%      -- printResults: if 0 the ptols will not be preinted to files 
%                       if 1 the ptols will be printed in the directory resutls/ 
%      -- GeneVars: if you know the Gene varariances (from PUMA)  give them here       
%
% Notes:
%       1. The confidence intervals are 95% (usually obtained with percentiles) 
%       2. In the current version the error bars of the TF profiles do not
%       contain likelihoood noise variances. If GeneVars are given then those
%       are plotted together with the observed gene expressions.    
%

%%%%%%  user defined parameters
%demdata = 'demEcoli';
%order = [1 5 3 4 2];
%order = 1:NumOfGenes;
%order = [1 5 3 4 2]; % for the Barenco data 
%order = 1:1:size(Y,1);  % for ecoli or other data
%DataID = 1; % 1 for Barenco data, 2 for Ecoli data
%bbar = 1;  % draw or not erroor when ploting basal , decay and sensitivities
%dataIDvasknownPuma = 1;  % plot the errors bard of the observed gene expressions if they have been
                         % computed by some package e.g. PUMA
%%%%%%  End of user defined parameters

%if strcmp(model.Likelihood.TFjointAct,'sigmoid') 
%   model.Likelihood.TFsingleAct = 'exp';
%end
dirr = '/usr/local/michalis/mlprojects/gpsamp/tex/diagrams/';

Genes = model.Likelihood.Genes;

if strcmp(model.constraints.sigmas,'fixed')
    GeneVars = model.Likelihood.sigmas;
end
 
TimesG = model.Likelihood.TimesG; 
TimesF = model.Likelihood.TimesF;
if isfield(model.Likelihood,'GenesTF')
    GenesTF = model.Likelihood.GenesTF;
    if strcmp(model.constraints.sigmasTF,'fixed')
        GeneTFVars = model.Likelihood.sigmasTF;
    end
end

TimesFF = TimesF(model.Likelihood.startTime:end);
ok = date;
fileName = [demdata 'MCMC' ok model.Likelihood.singleAct model.Likelihood.jointAct]; 

NumOfTFs = model.Likelihood.numTFs;

if strcmp(model.constraints.sigmas,'free') 
  % plots sigma2s and the lengghscales
  for j=1:model.Likelihood.numGenes
  sigma2j = squeeze(samples.sigmas(j,1,:));
  figure;
  hist(sigma2j,100);   
  if printResults
       print('-depsc', [dirr fileName 'Sigma2']);
  end
  titlestring = 'Observation variance: ';
  titlestring = [titlestring, num2str(j)]; 
  titlestring = [titlestring, ' gene'];
  title(titlestring,'fontsize', 20);
  end
end    

if isfield(model.Likelihood,'sigmasTF')
if strcmp(model.constraints.sigmasTF,'free') 
  % plots sigma2s and the lengghscales
  for j=1:model.Likelihood.numTFs
  sigma2j = squeeze(samples.sigmasTF(j,1,:));
  figure;
  hist(sigma2j,100);   
  if printResults
       print('-depsc', [dirr fileName 'SigmasTF']);
  end
  titlestring = 'Observation variance-TF Genes: ';
  titlestring = [titlestring, num2str(j)]; 
  titlestring = [titlestring, ' TF'];
  title(titlestring,'fontsize', 20);
  end
end  
end

% plot the lengthscales
for j=1:NumOfTFs
    figure;
    hist(squeeze(exp(2*samples.logthetas(j,1,:))),100); 
    %title('Lengthscale','fontsize', 20);
    if printResults
      print('-depsc', [dirr fileName 'LengthSc' 'TF' num2str(j)]);
    end
    titlestring = 'Lengthscale: ';
    titlestring = [titlestring, num2str(j)]; 
    titlestring = [titlestring, ' TF'];
    title(titlestring,'fontsize', 20);
end

NumOfGenes = model.Likelihood.numGenes;
order = 1:NumOfGenes;
NumOfReplicas = model.Likelihood.numReplicas;
NumOfSamples = size(samples.F,2);
SizF = size(samples.F{1},2);
TimesF = TimesF(:); 
for r=1:NumOfReplicas
  % 
  for j=1:NumOfTFs
     % 
     PredTF = zeros(NumOfSamples,SizF);    
     for t=1:NumOfSamples
         %FF(t,:) = samples.F{t}(j,:,r);
         lik = model.Likelihood;
         lik.kinetics = samples.kinetics(:,:,t); 
         if isfield(model.Likelihood,'GenesTF')         
            lik.kineticsTF = samples.kineticsTF(:,:,t);
         end
         PredTF(t,:) = gpmtfComputeTF(lik, samples.F{t}(j,:,r), j);
     end
    
     %mu = mean(PredTF)';
     %stds1 = sqrt(var(PredTF))';
     %stds2 = sqrt(var(PredTF))';
     %if strcmp(model.Likelihood.jointAct,'sigmoid')==1
     %FF = exp(FF);  
     mu = median(PredTF)';
     stds1 = (prctile(PredTF,95,1)'-mu)/2;
     stds2 = (mu-prctile(PredTF,5,1)')/2;
     %end
     
     figure
     plot(TimesF,mu,'b','lineWidth',3);
     hold on;
     fillColor = [0.7 0.7 0.7];
     %fillColor = [0.8 0.8 0.8];  % for the paper
     fill([TimesF; TimesF(end:-1:1)], [mu; mu(end:-1:1)]...
            + 2*[stds1; -stds2(end:-1:1)], fillColor,'EdgeColor',fillColor);
     plot(TimesF,mu,'b','lineWidth',3);
     
     axis([TimesF(1) TimesF(end)+0.1 0 max(mu+2*stds1)+0.1]);
     

    
     % plot the ground truth if exist
     if isfield(model,'groundtr') == 1
     FFgt = model.groundtr.TF(j,:,r);
     %FFgt = feval(model.Likelihood.TFsingleAct,model.GroundTruth.F(j,:,r));
     plot(TimesF,FFgt,'r','lineWidth',3);
     end
     
     if printResults
      print('-depsc', [dirr fileName 'Replica' num2str(r) 'TF' num2str(j)]);
     end 
     titlestring = 'Profile: ';
     titlestring = [titlestring, num2str(r)]; 
     titlestring = [titlestring, ' replica, '];
     titlestring = [titlestring, num2str(j)];
     titlestring = [titlestring, ' TF'];
     title(titlestring,'fontsize', 20);
     %
  end
  %
end


% plot predicted gene expressions 
for r=1:NumOfReplicas
  %  
  for j=1:NumOfGenes
     % 
     GG = zeros(NumOfSamples,model.Likelihood.sizTime);    
     for t=1:NumOfSamples
         GG(t,:) = samples.predGenes{t}(j,:,r);
     end
     
     mu = mean(GG)';
     %mu = mu(21:end);
     %mu = mu(1:2:end);
     stds = sqrt(var(GG))';
     %stds = stds(21:end);
     %stds = stds(1:2:end);
    
     TF = TimesFF'; % TimesFF(1:2:end)';
     %stds(stds>5)=5;
     %stds
     %mu(mu>9)=9;
     %mu(mu<0)=0;
     %pause
     figure
     plot(TF,mu,'b','lineWidth',r);
     hold on;
     fillColor = [0.7 0.7 0.7];
     %fillColor = [0.8 0.8 0.8];  % for the paper
     fill([TF; TF(end:-1:1)], [mu; mu(end:-1:1)]...
            + 2*[stds; -stds(end:-1:1)], fillColor,'EdgeColor',fillColor);
     plot(TF,mu,'b','lineWidth',3);
   
     plot(TimesG,Genes(j,:,r),'rx','markersize', 14','lineWidth', 2);
     if strcmp(model.constraints.sigmas,'fixed')
     errorbar(TimesG,  Genes(j,:,r), 2*sqrt(GeneVars(j,:,r)), 'rx','lineWidth', 1.5);
     end
     axis([min(TimesG(:))-0.1 max(TimesG(:))+0.1 0.95*min(min(Genes(j,:,r))) 1.05*max(max(Genes(j,:,r)))]);
     
     if printResults
      print('-depsc', [dirr fileName 'Replica' num2str(r) 'GeneExp' num2str(j)]);
     end
     titlestring = 'Expressions: ';
     titlestring = [titlestring, num2str(r)]; 
     titlestring = [titlestring, ' replica, '];
     titlestring = [titlestring, num2str(j)];
     titlestring = [titlestring, ' gene'];
     title(titlestring,'fontsize', 20);
     %
  end
  %
end

% predicted TF-Gene expressions 
if isfield(model.Likelihood,'GenesTF')
for r=1:NumOfReplicas
  %  
  for j=1:NumOfTFs
     % 
     GG = zeros(NumOfSamples,size(model.F,2));    
     for t=1:NumOfSamples
         GG(t,:) = samples.predGenesTF{t}(j,:,r);
     end
     
     mu = mean(GG)';
     %mu = mu(1:2:end);
     stds = sqrt(var(GG))';
     %stds = stds(1:2:end);
    
     TF = TimesF; % TimesFF(1:2:end)';
     figure
     plot(TF,mu,'b','lineWidth',r);
     hold on;
     fillColor = [0.7 0.7 0.7];
     %fillColor = [0.8 0.8 0.8];  % for the paper
     fill([TF; TF(end:-1:1)], [mu; mu(end:-1:1)]...
            + 2*[stds; -stds(end:-1:1)], fillColor,'EdgeColor',fillColor);
     plot(TF,mu,'b','lineWidth',3);
   
     plot(TimesG,GenesTF(j,:,r),'rx','markersize', 14','lineWidth', 2);
     if strcmp(model.constraints.sigmasTF,'fixed')
     errorbar(TimesG,  GenesTF(j,:,r), 2*sqrt(GeneTFVars(j,:,r)), 'rx','lineWidth', 1.5);
     end
     axis([min(TimesG(:))-0.1 max(TimesG(:))+0.1 0.95*min(min(GenesTF(j,:,r))) 1.05*max(max(GenesTF(j,:,r)))]);
     
     if printResults
      print('-depsc', [dirr fileName 'Replica' num2str(r) 'GeneTFExp' num2str(j)]);
     end
     titlestring = 'Expressions: ';
     titlestring = [titlestring, num2str(r)]; 
     titlestring = [titlestring, ' replica, '];
     titlestring = [titlestring, num2str(j)];
     titlestring = [titlestring, ' TF-gene'];
     title(titlestring,'fontsize', 20);
     %
  end
  %
end    
end

ok = mean(samples.kinetics,3);  
BB = squeeze(samples.kinetics(:,1,:))';
DD = squeeze(samples.kinetics(:,2,:))';
SS = squeeze(samples.kinetics(:,3,:))';
AA = squeeze(samples.kinetics(:,4,:))';
modelB = median(BB,1);
modelS = median(SS,1);
modelD = median(DD,1);
modelA = median(AA,1);

stdBB1 = prctile(BB,5);
stdDD1 = prctile(DD,5);
stdSS1 = prctile(SS,5);
stdBB2 = prctile(BB,95);
stdDD2 = prctile(DD,95);
stdSS2 = prctile(SS,95);
stdAA1 = prctile(AA,5);
stdAA2 = prctile(AA,95);



% Plot first basal transcription rates.
figure;
if isfield(model,'groundtr') == 1
bar([modelB(order); model.groundtr.kinetics(:,1)']', 0.7); colormap([0.9 0.9 0.9; 0 0 0]);
else
bar(modelB(order)', 0.7); colormap([0.9 0.9 0.9]);
end
hold on;
errorbar([1:NumOfGenes]-0.14, modelB(order), modelB(order)-stdBB1(order), stdBB2(order)-modelB(order),'.'); 
%errorbar([1:NumOfGenes], modelB(order), modelB(order)-stdBB1(order), stdBB2(order)-modelB(order),'.'); 
title('Basal rates','fontsize', 20);

if printResults
      print('-depsc', [dirr fileName 'Basal']);
end
     

% Plot the sensitivities.
figure;
if isfield(model,'groundtr') == 1
bar([modelS(order); model.groundtr.kinetics(:,3)']', 0.7); colormap([0.9 0.9 0.9; 0 0 0]);
else
bar(modelS(order)', 0.7); colormap([0.9 0.9 0.9]);
end
hold on;
errorbar([1:NumOfGenes]-0.14, modelS(order), modelS(order)-stdSS1(order), stdSS2(order)-modelS(order),'.');
%errorbar([1:NumOfGenes], modelS(order), modelS(order)-stdSS1(order), stdSS2(order)-modelS(order),'.');
title('Sensitivities','fontsize', 20);

if printResults
      print('-depsc', [dirr fileName 'Sensitivity']);
end

figure;
% plot degradation rates
if isfield(model,'groundtr') == 1
bar([modelD(order); model.groundtr.kinetics(:,2)']', 0.7); colormap([0.9 0.9 0.9; 0 0 0]);
else
bar(modelD(order)', 0.7); colormap([0.9 0.9 0.9]);
end
hold on;
errorbar([1:NumOfGenes]-0.14, modelD(order), modelD(order)-stdDD1(order), stdDD2(order)-modelD(order),'.');
%errorbar([1:NumOfGenes], modelD(order), modelD(order)-stdDD1(order), stdDD2(order)-modelD(order),'.');
title('Decays','fontsize', 20);

if printResults
      print('-depsc', [dirr fileName 'Decay']);
end

figure;
% plot initial conditions 
if isfield(model,'groundtr') == 1
bar([modelA(order); model.groundtr.kinetics(:,4)']', 0.7); colormap([0.9 0.9 0.9; 0 0 0]);
else
bar(modelA(order)', 0.7); colormap([0.9 0.9 0.9]);
end
hold on;
errorbar([1:NumOfGenes]-0.14, modelA(order), modelA(order)-stdAA1(order), stdAA2(order)-modelA(order),'.');
%errorbar([1:NumOfGenes], modelA(order), modelA(order)-stdAA1(order), stdAA2(order)-modelA(order),'.');
title('Initial conditions','fontsize', 20);

if printResults
      print('-depsc', [dirr fileName 'InitCond']);
end


if isfield(model.Likelihood,'GenesTF')
%
orderTF = 1:NumOfTFs;
ok = mean(samples.kineticsTF,3);  
DD = squeeze(samples.kineticsTF(:,1,:))';
SS = squeeze(samples.kineticsTF(:,2,:))';
modelS = median(SS,1);
modelD = median(DD,1);

stdDD1 = prctile(DD,5);
stdSS1 = prctile(SS,5);
stdDD2 = prctile(DD,95);
stdSS2 = prctile(SS,95);

figure;
% plot degradation rates
if isfield(model,'groundtr') == 1
bar([modelD(orderTF); model.groundtr.kineticsTF(:,1)']', 0.7); colormap([0.9 0.9 0.9; 0 0 0]);
else
bar(modelD(orderTF)', 0.7); colormap([0.9 0.9 0.9]);
end
hold on;
errorbar([1:NumOfTFs]-0.14, modelD(orderTF), modelD(orderTF)-stdDD1(orderTF), stdDD2(orderTF)-modelD(orderTF),'.');
%errorbar([1:NumOfGenes], modelD(order), modelD(order)-stdDD1(order), stdDD2(order)-modelD(order),'.');
title('TF-Genes Decays','fontsize', 20);

if printResults
      print('-depsc', [dirr fileName 'TFGeneDecay']);
end

% Plot the sensitivities.
figure;
if isfield(model,'groundtr') == 1
bar([modelS(orderTF); model.groundtr.kineticsTF(:,2)']', 0.7); colormap([0.9 0.9 0.9; 0 0 0]);
else
bar(modelS(orderTF)', 0.7); colormap([0.9 0.9 0.9]);
end
hold on;
errorbar([1:NumOfTFs]-0.14, modelS(orderTF), modelS(orderTF)-stdSS1(orderTF), stdSS2(orderTF)-modelS(orderTF),'.');
%errorbar([1:NumOfGenes], modelS(order), modelS(order)-stdSS1(order), stdSS2(order)-modelS(order),'.');
title('TF-Genes Sensitivities','fontsize', 20);

if printResults
      print('-depsc', [dirr fileName 'TFGeneSensitivity']);
end
%
end



for j=1:NumOfTFs
W1 = squeeze(samples.Weights(:,j,:))';
modelW1 = mean(W1,1);
stdW1_1 = sqrt(var(W1));
stdW1_2 = sqrt(var(W1));
% Plot first basal transcription rates.
figure;
if isfield(model,'groundtr') == 1
bar([modelW1(order); model.groundtr.W(:,j)']', 0.7); colormap([0.9 0.9 0.9; 0 0 0]);
else 
bar(modelW1(order)', 0.7); colormap([0.9 0.9 0.9]);
end
hold on;
errorbar([1:NumOfGenes]-0.14, modelW1(order), modelW1(order)-stdW1_1(order), stdW1_2(order)-modelW1(order),'.'); 
%errorbar([1:NumOfGenes], modelB(order), modelB(order)-stdBB1(order), stdBB2(order)-modelB(order),'.'); 
%title(j,'fontsize', 20);
if printResults
      print('-depsc', [dirr fileName 'IntWeights' 'TF' num2str(j)]);
end
titlestring = 'Interaction weights: '; 
titlestring = [titlestring, num2str(j)];
titlestring = [titlestring, ' TF'];
title(titlestring,'fontsize', 20);
end

W0 = samples.Weights0';
modelW0 = mean(W0,1);
stdW0_1 = sqrt(var(W0));
stdW0_2 = sqrt(var(W0));
figure;
% plot initial conditions
if isfield(model,'groundtr') == 1
bar([modelW0(order); model.groundtr.W0']', 0.7); colormap([0.9 0.9 0.9; 0 0 0]);
else
bar(modelW0(order)', 0.7); colormap([0.9 0.9 0.9]);
end
hold on;
errorbar([1:NumOfGenes]-0.14, modelW0(order), modelW0(order)-stdW0_1(order), stdW0_2(order)-modelW0(order),'.');
%errorbar([1:NumOfGenes], modelA(order), modelA(order)-stdAA1(order), stdAA2(order)-modelA(order),'.');
%title('W0','fontsize', 20);
if printResults
      print('-depsc', [dirr fileName 'IntBias']);
end
title('Interaction biases','fontsize', 20);


% plot delayes if they were inferred
if abs(model.Likelihood.tauMax) > 0
Taus = samples.Taus';
modelTaus = mean(Taus,1);
stdTaus_1 = prctile(Taus,5);%sqrt(var(Taus));
stdTaus_2 = prctile(Taus,95);%sqrt(var(Taus));
figure;
% plot initial conditions
if isfield(model,'groundtr') == 1
bar([modelTaus(order); model.groundtr.Taus]', 0.7); colormap([0.9 0.9 0.9; 0 0 0]);
else
bar(modelTaus(order)', 0.7); colormap([0.9 0.9 0.9]);
end
hold on;
errorbar([1:NumOfGenes]-0.14, modelTaus(order), modelTaus(order)-stdTaus_1(order), stdTaus_2(order)-modelTaus(order),'.');
%errorbar([1:NumOfGenes], modelA(order), modelA(order)-stdAA1(order), stdAA2(order)-modelA(order),'.');
%title('W0','fontsize', 20);
if printResults
      print('-depsc', [dirr fileName 'Delays']);
end
title('Delays','fontsize', 20);
end                                                                                 gpmtfSample.m                                                                                       0000700 0003466 0000024 00000101675 11312474302 013225  0                                                                                                    ustar   mtitsias                        games                                                                                                                                                                                                                  function [model PropDist samples accRates] = gpmtfSample(model, PropDist, trainOps)
% Description: Draw a set of samples from the Bayesian differential
%              equation model
%
% Inputs: 
%         -- model: the structure that contains the likelihood and GP
%                    parameters as well as the priors for all these
%                    quantities
%         -- PropDist: a stucture that defines the functional form of the proposal distribution
%         -- trainOps: user defined options about the burn-in and sampling iterations
%                      and others (see demos)
%
% Outputs: model: 
%         -- model: as above. The outputed model is updated to contain the
%                   parameters values of the final MCMC iteration
%                   parameters as well as the priors
%         -- PropDist: as above. PropDist can be updated (compared to the input one) 
%                     due to the update of the kernel parameters that
%                     influence the proposal 
%         -- samples: the structure that contrains the samples 
%         -- accRates: acceptance rates 
%

%ppi = 0.5;
BurnInIters = trainOps.Burnin; 
Iters = trainOps.T; 
StoreEvery = trainOps.StoreEvery;

Genes = model.Likelihood.Genes;
TimesG = model.Likelihood.TimesG; 
TimesF = model.Likelihood.TimesF; 
SizF = size(TimesF,2);
[NumOfGenes SizG NumOfReplicas] = size(Genes);


NumOfTFs = model.Likelihood.numTFs;
istart = ones(NumOfTFs,1);
for j=1:NumOfTFs
if model.constraints.Ft0(j)==0
   % do not sample the first control point so that the function 
   % will be fixed at the time t=0
   istart(j)=2;
end
end

%  check if the initial condition is fixed
fixInitCond = 0;
if strcmp(model.constraints.InitialConds_value,'fixed')==1 
    fixInitCond = 1;
end

% check if the observation noise is known/fixed
fixsigma2 = 0;
if strcmp(model.constraints.sigmas,'fixed')
    fixsigma2 = 1;
end

% check if the observation noise is known/fixed for the TF genes
if isfield(model.Likelihood,'GenesTF')
fixsigma2TF = 0;  
SizTFKin = size(model.Likelihood.kineticsTF,2);
if strcmp(model.constraints.sigmasTF,'fixed')
    fixsigma2TF = 1;
end
end


% check if the interaction weigths in the connectivity network are
% constrained to be positive
posw = 0; 
if strcmp(model.prior.weights.constraint,'positive')
    posw = 1;
end

% Construct the Comp binary matrix when you learn the 
% structure with Gibbs sampling --> not included in this version
netLearn = 0; 

% take the initial likelihood-kinetics parameters (defined out of this function)
LikParams = model.Likelihood;
SizKin = size(LikParams.kinetics,2);

% the latent function values are also special parameters that appear in both 
% the likelihood and the GP prior
F = model.F; 
%F = model.groundtr.F;

% store the control variables
Fu = model.Fu; % function values  
Xu = model.Xu; % locations  
M = model.M;
n = SizF;

% compute initial values for the log likelihood 
oldLogLik = zeros(NumOfReplicas, NumOfGenes);
% perform an evaluation of the log likelihood log p(Genes | F) 
if strcmp(model.constraints.replicas,'free')
   for r=1:NumOfReplicas
   %
   % evaluate the likelihood 
   [oldLogLik(r,:) predgen] = gpmtfLogLikelihoodGene(model.Likelihood, F(:,:,r), r, 1:NumOfGenes);
   PredictedGenes(:,:,r) = predgen;
  
   % additional likelihood when you have observations for the TF genes  
   if isfield(model.Likelihood,'GenesTF')      
      [oldLogLikTF(r,:) predgen] = gpmtfLogLikelihoodGeneTF(model.Likelihood, F(:,:,r), r, 1:NumOfTFs);
      PredictedGenesTF(:,:,r) = predgen;
   end
   %     
   end
else
   %
   % evaluate the likelihood for the first replica
   [oldLogLik(1,:) predgen] = gpmtfLogLikelihoodGene(model.Likelihood, F, 1, 1:NumOfGenes);
   PredictedGenes = predgen;
   % compute fast the additional likelihood when you have observations for the TF genes  
   if isfield(model.Likelihood,'GenesTF')      
      [oldLogLikTF(1,:) predgen] = gpmtfLogLikelihoodGeneTF(model.Likelihood, F, 1, 1:NumOfTFs);
      PredictedGenesTF = predgen;
   end
   %
   % the predicted genes are the same for the remaining coupled replicas
   for r=2:NumOfReplicas
      % compute fast the likelihood terms for the remaining replicas  
      oldLogLik(r,:) = remainRepsLikelihood(LikParams,  PredictedGenes, r, 1:NumOfGenes);
      % compute fast the additional likelihood when you have observations for the TF genes  
      if isfield(model.Likelihood,'GenesTF')
         oldLogLikTF(r,:) = remainRepsLikelihoodTF(LikParams, PredictedGenesTF, r, 1:NumOfTFs);
      end
   end
   %              
   %
end
%save ok PredictedGenesTF PredictedGenes;
%oldLogLikTF
%sum(oldLogLik(:))
%model.Likelihood.kineticsTF 
%model.groundtr.kineticsTF
%model.Likelihood.kinetics 
%model.groundtr.kinetics
%[model.Likelihood.Taus; model.groundtr.Taus]
%pause
% evaluation of the log prior for the kinetic parameters
lnpriorKin = ['ln',model.prior.kinetics.type,'pdf'];
TrspaceKin = model.prior.kinetics.priorSpace; 
Likkin = feval(TrspaceKin, LikParams.kinetics+eps);
oldLogPriorKin = feval(lnpriorKin, Likkin, model.prior.kinetics);
if isfield(model.Likelihood,'GenesTF')
  LikkinTF = feval(TrspaceKin, LikParams.kineticsTF+eps); 
  oldLogPriorKinTF = feval(lnpriorKin, LikkinTF, model.prior.kinetics);
end      
% evaluation of the prior for the interaction weights
lnpriorW = ['ln',model.prior.weights.type,'pdf'];
TrspaceW = model.prior.weights.priorSpace; 
LikW = feval(TrspaceW, [LikParams.W, LikParams.W0]);
oldLogPriorW = feval(lnpriorW, LikW, model.prior.weights);
% log prior of the lengthscale lengthscale 
lnpriorLengSc = ['ln',model.prior.GPkernel.lenghtScale.type,'pdf'];
%oldLogPriorLengSc = feval(lnpriorLengSc, exp(2*model.GP.logtheta(:,1)), model.prior.GPkernel.lenghtScale.a, model.prior.GPkernel.lenghtScale.b);
for j=1:NumOfTFs
oldLogPriorLengSc(j) = feval(lnpriorLengSc, 2*model.GP{j}.logtheta(1), model.prior.GPkernel.lenghtScale);
end

cnt = 0;
for j=1:NumOfTFs
    %
    if strcmp(model.constraints.replicas,'free')
      acceptF{j} = zeros(M(j),NumOfReplicas);
    else
      acceptF{j} = zeros(M(j),1);   
    end
    %
end
acceptKin = zeros(1,NumOfGenes);
acceptTFKin = zeros(1,NumOfTFs);
acceptW = zeros(1,NumOfGenes); 
acceptLengSc = zeros(1,NumOfTFs); 
%
for it = 1:(BurnInIters + Iters) 
    %
    %
    %F = model.groundtr.F;
    % Sample the TFs -----------------------------------------
    if 1 
    for j=1:NumOfTFs
    if strcmp(model.constraints.replicas,'free')   
         for r=1:NumOfReplicas
         %
         %
         Fold = F(j,:,r);
         Fuold = Fu{j}(:,r);
         % iterate between control points one-at-a-time 
         for i=istart(j):M(j)
         %
             % sample new control point i  
             %if abs(ppi) >= 1
                Fui = randn.*sqrt(PropDist.qF{j}.ku(i)) + PropDist.qF{j}.KInvK(i,:)*Fu{j}([1:i-1, i+1:end],r);    
             %else
             %   % use the underlaxed schme to sample the control point
             %   mu = PropDist.qF{j}.KInvK(i,:)*Fu{j}([1:i-1, i+1:end],r);
             %   Fui = mu  + ppi*(Fu{j}(i,r) -  mu)  + sqrt(1-ppi^2)*(randn.*sqrt(PropDist.qF{j}.ku(i)));
             %end    
             
             Funew = Fu{j}(:,r)';
             Funew(i) = Fui;
    
             % sample the remaining points 
             cmu = PropDist.qF{j}.cmuMinus + Funew*PropDist.qF{j}.KInvKu;
             Fnew = gaussianFastSample(1, cmu, PropDist.qF{j}.L);
   
             FFnew = F(:,:,r);
             FFnew(j,:) = Fnew;
   
             if ~isfield(model.Likelihood,'GenesTF')
             % perform an evaluation of the likelihood p(Genes | F) 
             [newLogLik predgen] = gpmtfLogLikelihoodGene(LikParams, FFnew, r, 1:NumOfGenes);
       
             % Metropolis-Hastings to accept-reject the proposal
             [accept, uprob] = metropolisHastings(sum(newLogLik(:)),sum(oldLogLik(r,:),2), 0, 0);
             else
             % perform an evaluation of the likelihood p(Genes | F) 
             [newLogLik predgen] = gpmtfLogLikelihoodGene(LikParams, FFnew, r, 1:NumOfGenes);        
             [newLogLikTF predgenTF] = gpmtfLogLikelihoodGeneTF(LikParams, FFnew, r, j); 
           
             % Metropolis-Hastings to accept-reject the proposal
             newL = sum(newLogLik(:)) + newLogLikTF;
             oldL = sum(oldLogLik(r,:),2) + oldLogLikTF(r,j);
             [accept, uprob] = metropolisHastings(newL, oldL, 0, 0);
             end
             %
             if (it<=BurnInIters) & trainOps.disp & (mod(it,50) == 0) 
             % 
                 visualize(model, F, Fu, FFnew, Funew, i, j, r);
             %  
             end
       
             %
             if (it > BurnInIters)
                 acceptF{j}(i,r) = acceptF{j}(i,r) + accept; 
             end
    
             % update protein F
            if accept == 1
                 F(j,:,r) = Fnew;
                 Fu{j}(:,r) = Funew';
                 PredictedGenes(:,:,r) = predgen;
                 oldLogLik(r,:) = newLogLik;
         
                 if isfield(model.Likelihood,'GenesTF')      
                    oldLogLikTF(r,j) = newLogLikTF;
                    PredictedGenesTF(j,:,r) = predgenTF;
                 end
            %   
            end
         end % num of control points loop
         %
         end % num Replicas loop
    else
         % Replicas are coupled
         % ---------------------------------------------
         % iterate between control points one-at-a-time 
         Fold = F(j,:);
         Fuold = Fu{j}(:,1);
         for i=istart(j):M(j)
         %
             % sample new control point i  
             %if abs(ppi) >= 1
                Fui = randn.*sqrt(PropDist.qF{j}.ku(i)) + PropDist.qF{j}.KInvK(i,:)*Fu{j}([1:i-1, i+1:end],1);        
             %else
             %   % use the underlaxed schme to sample the control point
             %   mu = PropDist.qF{j}.KInvK(i,:)*Fu{j}([1:i-1, i+1:end],1);
             %   Fui = mu  + ppi*(Fu{j}(i,1) -  mu)  + sqrt(1-ppi^2)*(randn.*sqrt(PropDist.qF{j}.ku(i)));
             %end    
       
             Funew = Fu{j}(:,1)';
             Funew(i) = Fui;
             
             % sample the remaining points 
             cmu = PropDist.qF{j}.cmuMinus + Funew*PropDist.qF{j}.KInvKu;
             Fnew = gaussianFastSample(1, cmu, PropDist.qF{j}.L);
             FFnew = F;
             FFnew(j,:) = Fnew;          
             
             newLogLik = zeros(NumOfReplicas,NumOfGenes);
             if ~isfield(model.Likelihood,'GenesTF')
                 % perform an evaluation of the likelihood p(Genes | F)      
                 [newLogLik(1,:) predgen] = gpmtfLogLikelihoodGene(LikParams, FFnew, 1, 1:NumOfGenes);
              
                 % computed faster the remaining likelihood terms  
                 for r=2:NumOfReplicas
                      newLogLik(r,:) = remainRepsLikelihood(LikParams,  predgen, r, 1:NumOfGenes);
                 end                 
                 % Metropolis-Hastings to accept-reject the proposal
                 [accept, uprob] = metropolisHastings(sum(newLogLik(:)),sum(oldLogLik(:)), 0, 0);
             else
                 % perform an evaluation of the likelihood p(Genes | F)      
                 [newLogLik(1,:) predgen] = gpmtfLogLikelihoodGene(LikParams, FFnew, 1, 1:NumOfGenes);   
                 [newLogLikTF(1) predgenTF] = gpmtfLogLikelihoodGeneTF(LikParams, FFnew, 1, j); 
           
                 % computed faster the remaining likelihood terms  
                 for r=2:NumOfReplicas
                     newLogLik(r,:) = remainRepsLikelihood(LikParams,  predgen, r, 1:NumOfGenes);
                     newLogLikTF(r) = remainRepsLikelihoodTF(LikParams, predgenTF, r, j);
                 end    
             
                 % Metropolis-Hastings to accept-reject the proposal
                 newL = sum(newLogLik(:)) + sum(newLogLikTF(:));
                 oldL = sum(oldLogLik(:)) + sum(oldLogLikTF(:,j),1);
                 [accept, uprob] = metropolisHastings(newL, oldL, 0, 0);
                 %
             end
             
             %
             if (it<=BurnInIters) & trainOps.disp & (mod(it,50) == 0) 
             % 
                 visualize(model, F, Fu, FFnew, Funew, i, j, 1);
             %  
             end
       
             %
             if (it > BurnInIters)
                 acceptF{j}(i) = acceptF{j}(i) + accept; 
             end
    
             % update protein F
            if accept == 1
                 F(j,:) = Fnew;
                 Fu{j} = Funew';
                 PredictedGenes = predgen;
                 oldLogLik = newLogLik;
         
                 if isfield(model.Likelihood,'GenesTF')      
                    oldLogLikTF(:,j) = newLogLikTF(:);
                    PredictedGenesTF(j,:) = predgenTF;
                 end
            %   
            end
         end % num of control points loop
         %
        end % num Replicas loop
        %
    end % num TFs loop
    end % zero-one if
   
    
    if 1
    % sample new kinetic parameters for each TF-gene  
    if isfield(model.Likelihood,'GenesTF')
    for j=1:NumOfTFs
        TFKineticsNew = randn(1,SizTFKin).*sqrt(PropDist.TFkin(j,:)) + log(LikParams.kineticsTF(j,:));
        TFKineticsNew(TFKineticsNew<-10) =-10; 
        TFKineticsNew(TFKineticsNew>10) = 10;
        TFKineticsNew = exp(TFKineticsNew); 
        %
        if model.constraints.geneTFsensitivity(j) == 0
        TFKineticsNew(2) = 1; 
        end
        
        LikParams1 = LikParams;
        LikParams1.kineticsTF(j,:)=TFKineticsNew; 
        newLogLik = []; 
        if strcmp(model.constraints.replicas,'free')   
           for r=1:NumOfReplicas
               % perform an evaluation of the likelihood p(GENES | TFs) 
               newLogLik(r,:) = gpmtfLogLikelihoodGene(LikParams1, F(:,:,r), r, 1:NumOfGenes);
               % 
           end
        else
           %
               % perform an evaluation of the likelihood p(GENES | TFs)
               % only for the first replica
               [newLogLik(1,:), predgen] = gpmtfLogLikelihoodGene(LikParams1, F, 1, 1:NumOfGenes);
               % computed faster the remaining likelihood terms  
               for r=2:NumOfReplicas
                   newLogLik(r,:) = remainRepsLikelihood(LikParams1, predgen, r, 1:NumOfGenes);
               end 
           %
        end
        %        
        Likkin = feval(TrspaceKin, TFKineticsNew);
        LogPriorKinNew = feval(lnpriorKin, Likkin, model.prior.kinetics);
        
        % Metropolis-Hastings to accept-reject the proposal
        newP = sum(newLogLik(:)) + sum(LogPriorKinNew(:));
        oldP = sum(oldLogLik(:)) + sum(oldLogPriorKinTF(j,:),2);
        [accept, uprob] = metropolisHastings(newP, oldP, 0, 0);
        %[accept, uprob] = metropolisHastings(newP, oldP, newLogProp, oldLogProp);
        if accept == 1
           LikParams.kineticsTF(j,:) = TFKineticsNew;
           oldLogLik = newLogLik; 
           oldLogPriorKinTF(j,:) = LogPriorKinNew;
        end
        %
        if (it > BurnInIters) 
           acceptTFKin(j) = acceptTFKin(j) + accept;
        end
        %
    end
    end
    end
    
    
    
    if 1
    % sample new kinetic parameters for each gene separately 
    for j=1:NumOfGenes
        KineticsNew = randn(1,SizKin).*sqrt(PropDist.kin(j,:)) + log(LikParams.kinetics(j,:));
        KineticsNew(KineticsNew<-10) =-10; 
        KineticsNew(KineticsNew>10) = 10;
        KineticsNew = exp(KineticsNew); 
        %
        % set the initial condition to be zero
        % if you know that it should be zero
        if model.constraints.InitialConds(j) == 0
        KineticsNew(4) = KineticsNew(1)/KineticsNew(2); 
        end
    
        LikParams1 = LikParams;
        LikParams1.kinetics(j,:)=KineticsNew; 
        newLogLik = zeros(1,NumOfReplicas);
        if strcmp(model.constraints.replicas,'free')
           %
           for r=1:NumOfReplicas
              % call the function only with j gene expressions  
              newLogLik(r) = gpmtfLogLikelihoodGene(LikParams1, F(:,:,r), r, j);
           end
        else
            % call the function only for the first replica
            [newLogLik(1), predgen] = gpmtfLogLikelihoodGene(LikParams1, F, 1, j);
            % computed faster the remaining likelihood terms  
            for r=2:NumOfReplicas
                newLogLik(r) = remainRepsLikelihood(LikParams1, predgen, r, j);
            end 
           %
        end
        %        
        Likkin = feval(TrspaceKin, KineticsNew);
        LogPriorKinNew = feval(lnpriorKin, Likkin, model.prior.kinetics);
        % Metropolis-Hastings to accept-reject the proposal
        oldP = sum(oldLogLik(:,j),1) + sum(oldLogPriorKin(j,:),2);
        newP = sum(newLogLik(:))+ sum(LogPriorKinNew(:)); 
        [accept, uprob] = metropolisHastings(newP, oldP, 0, 0);
        %[accept, uprob] = metropolisHastings(newP, oldP, newLogProp, oldLogProp);
        if accept == 1
           LikParams.kinetics(j,:) = KineticsNew;
           oldLogLik(:,j) = newLogLik(:); 
           oldLogPriorKin(j,:) = LogPriorKinNew;
        end
        %
        if (it > BurnInIters) 
           acceptKin(j) = acceptKin(j) + accept;
        end
        %
    end
    %
    end % if 0
    
    
    
    if 1 
    % sample the interaction weights 
    for j=1:NumOfGenes
        % 
        %  
        if posw == 1
            Wnew = randn(1,NumOfTFs+1).*sqrt(PropDist.W(j,:)) + log([LikParams.W(j,:), LikParams.W0(j)]+eps);     
            Wnew = exp(Wnew);
        else
            Wnew = randn(1,NumOfTFs+1).*sqrt(PropDist.W(j,:)) + [LikParams.W(j,:), LikParams.W0(j)];     
        end
        Wnew(1:NumOfTFs) = Wnew(1:NumOfTFs).*model.constraints.W(j,:);
       
        LikParams1 = LikParams;
        LikParams1.W(j,:) = Wnew(1:NumOfTFs);
        LikParams1.W0(j)=Wnew(end);
      
        newLogLik = zeros(1,NumOfReplicas);
        if strcmp(model.constraints.replicas,'free')
           for r=1:NumOfReplicas
               % call the function only with j gene expressions  
               newLogLik(r) = gpmtfLogLikelihoodGene(LikParams1, F(:,:,r), r, j);
           end
        else
           % call the function only for the first replica
           [newLogLik(1), predgen] = gpmtfLogLikelihoodGene(LikParams1, F, 1, j);
           % computed faster the remaining likelihood terms  
           for r=2:NumOfReplicas
               newLogLik(r) = remainRepsLikelihood(LikParams1, predgen, r, j);
           end 
        end
        
        LikW = feval(TrspaceW, Wnew+eps);
        LogPriorWnew = feval(lnpriorW, LikW, model.prior.weights);
        % Metropolis-Hastings to accept-reject the proposal
        oldP = sum(oldLogLik(:,j),1) + sum(oldLogPriorW(j,:),2);
        newP = sum(newLogLik(:)) + sum(LogPriorWnew(:)); 
        %
        [accept, uprob] = metropolisHastings(newP, oldP, 0, 0);
        if accept == 1
           LikParams.W(j,:) = Wnew(1:NumOfTFs);
           LikParams.W0(j) = Wnew(end);
           oldLogLik(:,j) = newLogLik(:); 
           oldLogPriorW(j,:) = LogPriorWnew;
        end
        %
        if (it > BurnInIters) 
           acceptW(j) = acceptW(j) + accept;
        end
        %
    end
    end % if 0
    
    
    % sample the noise of the likelihood when is free parameter.
    % The posterior is gamma so this step involves exact simualtion from
    % the gamma distribution
    if fixsigma2 == 0
       sumSquerrors1 = oldLogLik + 0.5*SizG*log(2*pi*repmat(LikParams.sigmas(:,1,1)',NumOfReplicas,1));
       sumSquerrors1 = -2*repmat(LikParams.sigmas(:,1,1)',NumOfReplicas,1).*sumSquerrors1;
       sumSquerrors = sum(sumSquerrors1,1);
                
       anew = model.prior.invsigma2.a + 0.5*NumOfReplicas*SizG;
       bnew = model.prior.invsigma2.b + 0.5*sumSquerrors;
       newinvsigma2 = gamrnd(anew,1./bnew);
       Nnewsigma2 = 1./newinvsigma2;
       LikParams.sigmas = repmat(Nnewsigma2(:),[1 SizG NumOfReplicas]);
       %for r=1:NumOfReplicas
       %oldLogLik(r,:) = genesLogLikelihood(LikParams, F(:,:,r), r, 1:NumOfGenes, Genes(:,:,r), TimesG, TimesF);
       %end  
       okk = repmat(LikParams.sigmas(:,1,1)',NumOfReplicas,1);
       oldLogLik = - 0.5*SizG*log(2*pi*okk) - (0.5./okk).*sumSquerrors1;
       % 
       %
    end
    
    % sample the noise of the TF-Genes likelihood when is free parameter.
    if isfield(model.Likelihood,'GenesTF')
    if fixsigma2TF == 0
       sumSquerrors1 = oldLogLikTF + 0.5*SizG*log(2*pi*repmat(LikParams.sigmasTF(:,1,1)',NumOfReplicas,1));
       sumSquerrors1 = -2*repmat(LikParams.sigmasTF(:,1,1)',NumOfReplicas,1).*sumSquerrors1;
       sumSquerrors = sum(sumSquerrors1,1);
                
       anew = model.prior.invsigma2.a + 0.5*NumOfReplicas*SizG;
       bnew = model.prior.invsigma2.b + 0.5*sumSquerrors;
       newinvsigma2 = gamrnd(anew,1./bnew);
       Nnewsigma2 = 1./newinvsigma2;
       LikParams.sigmasTF = repmat(Nnewsigma2(:),[1 SizG NumOfReplicas]);
       %for r=1:NumOfReplicas
       %oldLogLik(r,:) = genesLogLikelihood(LikParams, F(:,:,r), r, 1:NumOfGenes, Genes(:,:,r), TimesG, TimesF);
       %end  
       okk = repmat(LikParams.sigmasTF(:,1,1)',NumOfReplicas,1);
       oldLogLikTF = - 0.5*SizG*log(2*pi*okk) - (0.5./okk).*sumSquerrors1;
       % 
       %
    end
    end

    if 1 
    if model.Likelihood.tauMax < 0
    % sample the delay parameters in the ODEs 
    for j=1:NumOfGenes
        LikParams1 = LikParams;
        if LikParams1.Tausindex(j) == 1
            % propose to decrease tau
            LikParams1.Tausindex(j) = 2;
            LikParams1.Taus(j) = LikParams1.Taus(j) - LikParams1.step;
            logbias = log(0.5);
            %
        elseif LikParams1.Tausindex(j) == LikParams1.startTime
            % propose to increase tau 
            LikParams1.Tausindex(j) = LikParams1.startTime-1;
            LikParams1.Taus(j) = LikParams1.Taus(j) + LikParams1.step;
            logbias = -log(0.5);
            %
        else
            %
            % propose to decrease or increase with probability 0.5
            ud = round(rand); 
            ud(ud==0)=-1;
            logbias = 0;
            LikParams1.Tausindex(j) = LikParams1.Tausindex(j)-ud;
            LikParams1.Taus(j) = LikParams1.Taus(j) + ud*LikParams1.step;
        end
        
        newLogLik = zeros(1,NumOfReplicas);
        if strcmp(model.constraints.replicas,'free')
           for r=1:NumOfReplicas
               % call the function only with j gene expressions  
               newLogLik(r) = gpmtfLogLikelihoodGene(LikParams1, F(:,:,r), r, j);
           end
        else
           % call the function only for the first replica
           [newLogLik(1), predgen] = gpmtfLogLikelihoodGene(LikParams1, F, 1, j);
           % computed faster the remaining likelihood terms  
           for r=2:NumOfReplicas
               newLogLik(r) = remainRepsLikelihood(LikParams1, predgen, r, j);
           end  
        end
        
        oldP = log(model.prior.delays.prob(LikParams.Tausindex(j))) + sum(oldLogLik(:,j),1) - logbias;
        newP = log(model.prior.delays.prob(LikParams1.Tausindex(j))) + sum(newLogLik(:)); 
        [accept, uprob] = metropolisHastings(newP, oldP, 0, 0);
        if accept == 1
           LikParams.Tausindex(j) = LikParams1.Tausindex(j);
           LikParams.Taus(j) = LikParams1.Taus(j);
           oldLogLik(:,j) = newLogLik(:);
        end
        %
        %
    end
    end
    end
   
   
    %
    %
    if 1
    % sample the lenghtscale of the TFs 
    for j=1:NumOfTFs
        %
        % to samples the lengthscale you need to evaluate 
        % the GP pior compute the new and the current/old GP prior 
        newlogEll = randn.*sqrt(PropDist.LengSc(j)) + 2*model.GP{j}.logtheta(1);
        %newlogEll = model.GP{j}.logtheta(1);
        newEll2 = exp(newlogEll);  
        newK = exp(-(0.5/newEll2)*model.GP{j}.X2) + ...
               exp(2*model.GP{j}.logtheta(end))*eye(size(model.GP{j}.X2,1)); 
        % compute the Cholesky decomposition of the new K
        [newL,er]=jitterChol(newK);
        newL = newL';
        % evaluate the new log GP prior value 
        invnewL = newL\eye(SizF+M(j));
        newLogDetK = 2*sum(log(diag(newL)));
        if strcmp(model.constraints.replicas,'free')
            %
            newlogGP = - 0.5*NumOfReplicas*newLogDetK;
            oldlogGP = - 0.5*NumOfReplicas*PropDist.qF{j}.LogDetK;
            for r=1:NumOfReplicas     
                temp = invnewL*([F(j,:,r), Fu{j}(:,r)']'); 
                newlogGP = newlogGP - 0.5*temp'*temp;
                temp = PropDist.qF{j}.invL*([F(j,:,r), Fu{j}(:,r)']'); 
                oldlogGP = oldlogGP - 0.5*temp'*temp;
            end
            %
        else
            newlogGP = - 0.5*newLogDetK;
            oldlogGP = - 0.5*PropDist.qF{j}.LogDetK;     
            temp = invnewL*([F(j,:), Fu{j}(:,1)']'); 
            newlogGP = newlogGP - 0.5*temp'*temp;
            temp = PropDist.qF{j}.invL*([F(j,:), Fu{j}(:,1)']'); 
            oldlogGP = oldlogGP - 0.5*temp'*temp;
        end
            
        LogPriorLengScnew = feval(lnpriorLengSc, newlogEll, model.prior.GPkernel.lenghtScale);
        % Metropolis-Hastings to accept-reject the proposal
        oldlogGP = oldlogGP + oldLogPriorLengSc(j);
        newlogGP = newlogGP + LogPriorLengScnew; 
        
        
        %YnewlogEll = 2*model.GP{j}.logtheta(1);
        %%newlogEll = model.GP{j}.logtheta(1);
        %YnewEll2 = exp(YnewlogEll);  
        %YnewK = kernCompute(model.GP{j}, [TimesF(:); Xu(:)]);
        %%ok     = exp(-(0.5/YnewEll2)*model.GP{1}.X2) + ...
        %%        exp(2*model.GP{j}.logtheta(end))*eye(size(model.GP{1}.X2,1));
            
           
        %% compute the Cholesky decomposition of the new K
        %[YnewL,er]=jitterChol(YnewK);
        %YnewL = YnewL';
        %% evaluate the new log GP prior value 
        %YinvnewL = YnewL\eye(SizF+M);
        %YnewLogDetK = 2*sum(log(diag(YnewL)));
        %YnewlogGP = - 0.5*NumOfReplicas*YnewLogDetK;
        %for r=1:NumOfReplicas     
        %   temp = YinvnewL*([F(j,:,r), Fu(j,:,r)]'); 
        %   YnewlogGP = YnewlogGP - 0.5*temp'*temp;
        %end
        %
        %YLogPriorLengScnew = feval(lnpriorLengSc, YnewlogEll, model.prior.GPkernel.lenghtScale.a, model.prior.GPkernel.lenghtScale.b); 
        %YnewlogGP = YnewlogGP + oldLogPriorLengSc(j);
        
        %if mod(it,100) == 0
        %[oldlogGP  newlogGP]
        %[oldLogPriorLengSc(j) LogPriorLengScnew]
        %[exp(2*model.GP{j}.logtheta(1)) newEll2]
        %end
        %
        
        [accept, uprob] = metropolisHastings(newlogGP, oldlogGP, 0, 0);
        %%%%%%%%%%%%%%%%  start accept/update proposal for the lengthscale %%%%%%%%%%%% 
        if accept == 1
           U = n+1:n+M(j);
           model.GP{j}.logtheta(1) = newlogEll/2;
           oldLogPriorLengSc(j) = LogPriorLengScnew;
           %newEll2
           %pause
           PropDist.qF{j}.K = newK;
           PropDist.qF{j}.invL = invnewL; 
           PropDist.qF{j}.LogDetK = newLogDetK;
           
           % change the proposal distribution for the TF
           % compute the conditional GP prior given the control variables
           [cmuMinus, cSigma, KInvKu] = gaussianFastConditional(PropDist.qF{j}.m', newK, 1:n, U);
           [L,er]=jitterChol(cSigma);
           if er>0, L = real(sqrtm(cSigma)); end
           PropDist.qF{j}.cmuMinus = cmuMinus; 
           PropDist.qF{j}.cSigma = cSigma;
           PropDist.qF{j}.KInvKu = KInvKu;
           PropDist.qF{j}.L = L;
           for i=1:M(j)
           %  
              G = [1:i-1, i+1:M(j)];  
              [alpha(i), ku(i), KInvK(i,:)] = gaussianFastConditional(PropDist.qF{j}.m(U)', PropDist.qF{j}.K(U,U), i, G);
           %
           end
           PropDist.qF{j}.alpha = alpha;
           PropDist.qF{j}.ku = ku;
           PropDist.qF{j}.KInvK = KInvK;
           clear alpha  ku  KInvK;
        end
        % 
        if (it > BurnInIters) 
           acceptLengSc(j) = acceptLengSc(j) + accept;
        end
        %%%%%%%%%%%%%%%%%%%%%%% end accept/update proposal for the lengthscale %%%%%%%%%%%%%%%%
        %
    end
    end
    
    %
    %
    % keep samples after burn in
    if (it > BurnInIters)  & (mod(it,StoreEvery) == 0)
        %
        cnt = cnt + 1;
        samples.F{cnt} = F;
        samples.Fu{cnt} = Fu;
        %samples.predGenes{cnt} = PredictedGenes;
        samples.kinetics(:,:,cnt) = LikParams.kinetics;
        samples.Weights(:,:,cnt) = LikParams.W;
        samples.Weights0(:,cnt) = LikParams.W0;
        %
        if model.Likelihood.tauMax < 0
           samples.Taus(:,cnt) = LikParams.Taus(:);
           samples.Tausindex(:,cnt) = LikParams.Tausindex(:);
        end
        %
        if isfield(model.Likelihood,'GenesTF')
            samples.kineticsTF(:,:,cnt) = LikParams.kineticsTF;
            samples.LogLTF(cnt) = sum(oldLogLikTF(:));
            %samples.predGenesTF{cnt} = PredictedGenesTF;
            if fixsigma2TF == 0
                samples.sigmasTF(:,:,cnt) = LikParams.sigmasTF(:,:,1);
            end
        end
        if fixsigma2 == 0
           samples.sigmas(:,:,cnt) = LikParams.sigmas(:,:,1);
        end
        %if netLearn == 1
        %    samples.NetX(:,:,cnt) = LikParams.Net_X;
        %end
        for jin=1:NumOfTFs
            samples.logthetas(jin,cnt) = model.GP{jin}.logtheta(1);
        end
        samples.LogL(cnt) = sum(oldLogLik(:));
        %save(trainOps.ResFile,'samples','model');
        %
    end
    %
    %        
end

% Before you return store the final state of the Markov chain to 
% the model structure
model.Likelihood.kinetics = LikParams.kinetics;
model.Likelihood.W = LikParams.W;
model.Likelihood.W0 = LikParams.W0;
model.Likelihood.Tausindex = LikParams.Tausindex;
model.Likelihood.Taus = LikParams.Taus;
if isfield(model.Likelihood,'GenesTF')
    model.Likelihood.kineticsTF = LikParams.kineticsTF;
    model.Likelihood.sigmasTF = LikParams.sigmasTF;
end

model.Likelihood.sigmas = LikParams.sigmas;
if netLearn == 1
model.Likelihood.Net_x = LikParams.Net_X;
end
model.F = F;
model.Fu = Fu;
%
for j=1:NumOfTFs
    accRates.F{j} = (acceptF{j}/Iters)*100; 
    if istart == 2
       if strcmp(model.constraints.replicas,'free')  
          accRates.F{j}(1,:) = 100*ones(1,NumOfReplicas);
       else
          accRates.F{j}(1) = 100; 
       end
    end
end
%

accRates.Kin = (acceptKin/Iters)*100;
accRates.W = (acceptW/Iters)*100;
accRates.LengSc = (acceptLengSc/Iters)*100;
if isfield(model.Likelihood,'GenesTF')
   accRates.TFKin = (acceptTFKin/Iters)*100;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% visualization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
function visualize(model,F,Fu,FFnew,Funew,i,j,r)
%
trform = 'lin'; % model.Likelihood.singleAct; % 'exp' or 'linear'
if strcmp(model.constraints.replicas,'free')
subplot(1,model.Likelihood.numReplicas,r);    
end
Futmp = feval(model.Likelihood.singleAct,Fu{j}(:,r)');
Ftmp = feval(model.Likelihood.singleAct,F(j,:,r));
Funewtmp = feval(model.Likelihood.singleAct,Funew);
FFnewtmp = feval(model.Likelihood.singleAct,FFnew);
%
plot(model.Xu{j}, feval(trform,Futmp),'or','MarkerSize', 14,'lineWidth', 3);
hold on;
plot(model.Likelihood.TimesF, feval(trform,Ftmp),'g','lineWidth',4);
if isfield(model,'groundtr') == 1
    GrFtmp = feval(model.Likelihood.singleAct,model.groundtr.F(j,:));
    plot(model.Likelihood.TimesF, feval(trform,GrFtmp),'k','lineWidth',4);
end
title(j);
pause(0.3);
plot(model.Xu{j}(i), feval(trform,Futmp(i)),'oy','MarkerSize', 14,'lineWidth', 3);
plot(model.Xu{j}(i), feval(trform,Funewtmp(i)), 'md','MarkerSize', 14, 'lineWidth',3);
plot(model.Likelihood.TimesF, feval(trform,FFnewtmp(j,:)), '--b', 'lineWidth', 4);
pause(0.5);
hold off;
%%%%%%%%%%%%%%%%%%%%%%%%%%% end of visualization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function loglikval = remainRepsLikelihood(LikParams,  PredictedGenes, r, Gindex)
%
    Gns = LikParams.Genes(:,:,r);
    loglikval = - 0.5*sum(log(2*pi*LikParams.sigmas(Gindex,:,r)),2)....
                - 0.5*sum(((Gns(Gindex,:) - PredictedGenes(:,LikParams.comInds)).^2)./LikParams.sigmas(Gindex,:,r),2);
    loglikval = loglikval';
% 
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function loglikvalTF = remainRepsLikelihoodTF(LikParams, PredictedGenesTF, r, TFindex)
%    
  
  GnsTF = LikParams.GenesTF(:,:,r);
  loglikvalTF = - 0.5*sum(log(2*pi*LikParams.sigmasTF(TFindex,:,r)),2)....
                       - 0.5*sum(((GnsTF(TFindex,:) - PredictedGenesTF(:,LikParams.comIndsTF)).^2)./LikParams.sigmasTF(TFindex,:,r),2);
  loglikvalTF = loglikvalTF';
%                                                                     gpmtfSummariseResults.m                                                                             0000600 0003466 0000024 00000002520 11321414236 015317  0                                                                                                    ustar   mtitsias                        games                                                                                                                                                                                                                  function results = gpmtfSummariseResults(testGenes, method),

if nargin < 2,
  method='zscore';
end

if strcmp(method, 'loglike'),
  results = zeros(1, length(testGenes));
  for k=1:length(testGenes),
    results(k) = mean(testGenes{k}.LogL);
  end
  return;
elseif strcmp(method, 'harmmeanlike'),
  results = zeros(1, length(testGenes));
  for k=1:length(testGenes),
    results(k) = 1./mean(1 ./ testGenes{k}.LogL);
  end
  return;
elseif strcmp(method, 'meansigma'),
  results = zeros(1, length(testGenes));
  if isfield(testGenes{1}, 'sigmas'),
    for k=1:length(testGenes),
      results(k) = mean(testGenes{k}.sigmas);
    end
  else
    for k=1:length(testGenes),
      results(k) = mean(testGenes{k}.kinetics(2, :));
    end
  end
  return;
end

results = zeros(size(testGenes{1}.Weights, 1), length(testGenes));

for k=1:length(testGenes),
  W = testGenes{k}.Weights;
  d = size(W, 1);
  switch method,
   case 'meansigmaweight',
    if isfield(testGenes{k}, 'sigmas'),
      results(:, k) = mean(W .* repmat(testGenes{K}.sigmas, [d, 1]), 2);
    else
      results(:, k) = mean(W .* repmat(testGenes{K}.kinetics(2, :), [d, 1]), 2);
    end
   case 'meanweight',
    results(:, k) = mean(W, 2);
   case 'zscore',
    results(:, k) = mean(W, 2) ./ std(W, 0, 2);
   case 'pplr2',
    results(:, k) = max(mean(W < -.02, 2), mean(W > .02, 2));
  end
end
                                                                                                                                                                                gpmtfTestGenesAdapt.m                                                                               0000600 0003466 0000024 00000020116 11311755521 014650  0                                                                                                    ustar   mtitsias                        games                                                                                                                                                                                                                  function [model PropDist samples accRates] = gpmtfTestGenesAdapt(model, AdaptOps)
%[model PropDist samples accRates] = gpmtfAdapt(model, AdaptOps)
%
% Description: Sample the parameters of the Bayesian differential equation 
%              model so that to tune/adapt the proposal distribution (number 
%              of control points, variances of Gaussian proposals etc). 
%                ** This function should always be called before 
%                   gpTFSample. When the gpTFSample will be called the 
%                   proposal distribution (previously adapted)
%                   is kept fixed. **
%
% Inputs: 
%    -- model: the structure that contains the likelihood and GP
%              parameters as well as the priors for all these
%              quantities; write 'help model' for full details. 
%    -- PropDist: a stucture that defines the functional form of the 
%                 proposal distribution. It contains the following fields
%              * qF: contains precomputed variables needed for the 
%                    definition of the poposal distribution (that uses 
%                    control variables) over the GP latent functions 
%              * qKinVars: NumOfGenes x 4 matrix that contains the
%                    variances of the Gaussian proposal distributions
%                    used to sample the logarithm of the kinetic parameters 
%                    (represented in the log space; see modelCreate)
%                    for each gene. These Gaussian proposals have diagonal 
%                    covariance matrices so as each row corresponds to each 
%                    gene, 1 column to the basal rates (B), 2 column to 
%                    the decay rates (D), 3 column to the sensitivities(S)
%                    and 4 column to the initial condition parameter (A). 
%              * qWeigVars: NumOfGenes x M matrix that contains the
%                    variances of the Gaussian proposal distribution
%                    of all the parameters that exist in the activation 
%                    function of the latent GP function (that represent the 
%                    log of the TFs). Each row of the matrix corresponds to 
%                    each gene. The number M of parameters depends on the 
%                    number of TFs used and the functional form of the 
%                    activation function. When the sigmoid activation 
%                    function is used then  M=NumOfTFs+1. When the 
%                    Michaelis-Menten that is valid only for a single 
%                    TF (NumOfTFs=1), then M = 1; the single parameter 
%                    corresponds to the gamma Michaelis-Menten constant
%              * qLengScVars: NumOfTFs x 1 vector that contains all the
%                    variances of the Gaussian proposal distributions used to 
%                    sample the logarithm of the lengthscales of the NumOfTFs
%                    different rbf GP priors (see modelCreate) 
%    -- Genes : NumOfGenes x NumOfTimes x Replicas that stores the
%               gene expressions for all genes, all times and 
%               replicas
%    -- TimesG: The time points where gene expression are evaluated 
%    -- TimesF: The times where the GP function are evaluated
%               TimesF>>TimesG
%    -- trainOps: User defined options about the burn-in and sampling
%               iterations, visualization features and others 
%               **see demos**
%
% Outputs: 
%    -- model: see above or write 'help model' for full details. 
%              The outputed model is updated to contain the parameters 
%              values of the final MCMC iteration
%    -- PropDist: as above. the precomputations in the PropDist.qF field
%              can be different (compared to the ones given at the input) 
%              due the update of the kernel lengthscale that determine 
%              the proposal over the latent GP functions
%    -- samples: the structure that contains the samples. In contains
%              
%              
%
%    -- accRateF: acceptance rates during sampling time (after burn-in) 
%                 for the TFs
%    -- accRateKin: >>  >> for the kinetic parameters 
%    -- accRateW: >>  >> for the interaction weigths and any other 
%              parameters that exist in the activation function of 
%              the TF (e.g. the gamma in the Michaels-Menten activation)  
%    -- accRateLengSc:  >> >> for the lengthscales of the Gaussian
%              kernel 

BurnInIters = AdaptOps.Burnin; 
Iters = AdaptOps.T; 

%%% SIZES 
[NumOfGenes SizG NumReplicas] = size(model.Likelihood.Genes);
% number of genes 
NumOfGenes = model.Likelihood.numGenes;
% number of times the gene expression is evaluated
SizG = model.Likelihood.numTimes;
% number of replicas per gene 
NumReplicas = model.Likelihood.numReplicas;
% number of transcription factors
NumOfTFs = model.Likelihood.numTFs;
SizKin = size(model.Likelihood.kinetics,2);
SizF = size(model.Likelihood.TimesF,2);


% Initial proposal Gaussian distribution (with diagonal covariance matrices) 
% for the kinetic parameters interection weights and the lengthscale of the GPs 
PropDist.kin = 0.05*ones(NumOfGenes,SizKin);
% interaction weigths and bias 
PropDist.W = 0.05*ones(NumOfGenes,NumOfTFs+1);

% useful ranges needed in the adaption of the 
% variances of theese proposal distribution 
qKinBelow = 0.000001; qKinAbove = 2;
qWbelow = 0.000001;   qWabove = 2;
epsilon = 0.1;

cnt = 0;
%
% do the adaption 
while 1
%
%  
   [model PropDist samples accRates] = gpmtfTestGenesSample(model, PropDist, AdaptOps);
 
   accRateKin = accRates.Kin;
   accRateW = accRates.W;
   
   if AdaptOps.disp == 1 
   fprintf(1,'------ ADAPTION STEP #%2d ------ \n',cnt+1); 
   fprintf(1,'Acceptance Rates for kinetic parameters (per gene))\n');
   disp(accRateKin);
 
   fprintf(1,'Acceptance Rates for Interaction weights (per gene)\n');
   disp(accRateW);
   fprintf(1,'Average likelihood value %15.8f',mean(samples.LogL));
   fprintf(1,'\n');
   fprintf(1,'------------------------------- \n',cnt+1);
   end 
   
   % if you got a good acceptance rate, then stop
   if  (min(accRateKin(:))>15) & (min(accRateW(:))>15) 
        disp('END OF ADAPTION: acceptance rates OK');
        %pause
        break;
   end
    
   cnt = cnt + 1;
   % do not allow more than 50 iterations when you adapt the proposal distribution
   if cnt == 50
       warning('END OF ADAPTION: acceptance rates were not all OK');
       break;
   end
   
   %%%%%%%%%%%%%%%%%%%%%%% START of ADAPT KINETICS PROPOSAL %%%%%%%%%%%%%%%%
   % adapt the proposal over the kinetic parameters (desired acceptance rate: 15-35%)
   for j=1:NumOfGenes
      if accRateKin(j) > 35
         % incease the covariance to reduce the acceptance rate
         PropDist.kin(j,:) = PropDist.kin(j,:) + epsilon*PropDist.kin(j,:);
         if PropDist.kin(j,1) > qKinAbove 
             PropDist.kin(j,:) = qKinAbove*ones(1,SizKin);
         end
      end
      if accRateKin(j) < 15
         % decrease the covariance to incease the acceptance rate
         PropDist.kin(j,:) = PropDist.kin(j,:) - epsilon*PropDist.kin(j,:);    
         if PropDist.kin(j,1) < qKinBelow 
             PropDist.kin(j,:) = qKinBelow*ones(1,SizKin);
         end
         %
      end
       %
   end
   %%%%%%%%%%%%%%%%%%%%%%% END of ADAPT KINETICS PROPOSAL %%%%%%%%%%%%%%%%
         
   %%%%%%%%%%%%%%%%%%%%%%% START of ADAPT WEIGHTS PROPOSAL %%%%%%%%%%%%%%%%
   % adapt the proposal over the interaction weights (desired acceptance rate: 15-35%)
   for j=1:NumOfGenes
      if accRateW(j) > 35
         % incease the covariance to reduce the acceptance rate
         PropDist.W(j,:) = PropDist.W(j,:) + epsilon*PropDist.W(j,:);
         if PropDist.W(j,1) > qWabove 
             PropDist.W(j,:) = qWabove*ones(1,NumOfTFs+1);
         end
      end
      if accRateW(j) < 15
         % decrease the covariance to incease the acceptance rate
         PropDist.W(j,:) = PropDist.W(j,:) - epsilon*PropDist.W(j,:);    
         if PropDist.W(j,1) < qWbelow 
             PropDist.W(j,:) = qWbelow*ones(1,NumOfTFs+1);
         end
         %
      end
       %
   end
   %%%%%%%%%%%%%%%%%%%%%%% END of ADAPT WEIGHTS PROPOSAL %%%%%%%%%%%%%%%%
   
%
%
end
                                                                                                                                                                                                                                                                                                                                                                                                                                                  gpmtfTestGenesAdapt2.m                                                                              0000600 0003466 0000024 00000021567 11312010607 014733  0                                                                                                    ustar   mtitsias                        games                                                                                                                                                                                                                  function [model PropDist samples accRates] = gpmtfTestGenesAdapt2(model, TFs, simMat,  AdaptOps)
%[model PropDist samples accRates] = gpmtfAdapt(model, AdaptOps)
%
% Description: Sample the parameters of the Bayesian differential equation 
%              model so that to tune/adapt the proposal distribution (number 
%              of control points, variances of Gaussian proposals etc). 
%                ** This function should always be called before 
%                   gpTFSample. When the gpTFSample will be called the 
%                   proposal distribution (previously adapted)
%                   is kept fixed. **
%
% Inputs: 
%    -- model: the structure that contains the likelihood and GP
%              parameters as well as the priors for all these
%              quantities; write 'help model' for full details. 
%    -- PropDist: a stucture that defines the functional form of the 
%                 proposal distribution. It contains the following fields
%              * qF: contains precomputed variables needed for the 
%                    definition of the poposal distribution (that uses 
%                    control variables) over the GP latent functions 
%              * qKinVars: NumOfGenes x 4 matrix that contains the
%                    variances of the Gaussian proposal distributions
%                    used to sample the logarithm of the kinetic parameters 
%                    (represented in the log space; see modelCreate)
%                    for each gene. These Gaussian proposals have diagonal 
%                    covariance matrices so as each row corresponds to each 
%                    gene, 1 column to the basal rates (B), 2 column to 
%                    the decay rates (D), 3 column to the sensitivities(S)
%                    and 4 column to the initial condition parameter (A). 
%              * qWeigVars: NumOfGenes x M matrix that contains the
%                    variances of the Gaussian proposal distribution
%                    of all the parameters that exist in the activation 
%                    function of the latent GP function (that represent the 
%                    log of the TFs). Each row of the matrix corresponds to 
%                    each gene. The number M of parameters depends on the 
%                    number of TFs used and the functional form of the 
%                    activation function. When the sigmoid activation 
%                    function is used then  M=NumOfTFs+1. When the 
%                    Michaelis-Menten that is valid only for a single 
%                    TF (NumOfTFs=1), then M = 1; the single parameter 
%                    corresponds to the gamma Michaelis-Menten constant
%              * qLengScVars: NumOfTFs x 1 vector that contains all the
%                    variances of the Gaussian proposal distributions used to 
%                    sample the logarithm of the lengthscales of the NumOfTFs
%                    different rbf GP priors (see modelCreate) 
%    -- Genes : NumOfGenes x NumOfTimes x Replicas that stores the
%               gene expressions for all genes, all times and 
%               replicas
%    -- TimesG: The time points where gene expression are evaluated 
%    -- TimesF: The times where the GP function are evaluated
%               TimesF>>TimesG
%    -- trainOps: User defined options about the burn-in and sampling
%               iterations, visualization features and others 
%               **see demos**
%
% Outputs: 
%    -- model: see above or write 'help model' for full details. 
%              The outputed model is updated to contain the parameters 
%              values of the final MCMC iteration
%    -- PropDist: as above. the precomputations in the PropDist.qF field
%              can be different (compared to the ones given at the input) 
%              due the update of the kernel lengthscale that determine 
%              the proposal over the latent GP functions
%    -- samples: the structure that contains the samples. In contains
%              
%              
%
%    -- accRateF: acceptance rates during sampling time (after burn-in) 
%                 for the TFs
%    -- accRateKin: >>  >> for the kinetic parameters 
%    -- accRateW: >>  >> for the interaction weigths and any other 
%              parameters that exist in the activation function of 
%              the TF (e.g. the gamma in the Michaels-Menten activation)  
%    -- accRateLengSc:  >> >> for the lengthscales of the Gaussian
%              kernel 

BurnInIters = AdaptOps.Burnin; 
Iters = AdaptOps.T; 

%%% SIZES 
[NumOfGenes SizG NumReplicas] = size(model.Likelihood.Genes);
% number of genes 
NumOfGenes = model.Likelihood.numGenes;
% number of times the gene expression is evaluated
SizG = model.Likelihood.numTimes;
% number of replicas per gene 
NumReplicas = model.Likelihood.numReplicas;
% number of transcription factors
NumOfTFs = model.Likelihood.numTFs;
SizKin = size(model.Likelihood.kinetics,2);
SizF = size(model.Likelihood.TimesF,2);


% initialize the transcription factors 
if strcmp(model.constraints.replicas,'free')
    %
    % this variable will never used (its jsut to cheat the fucntino compute likelihoods
    % which will use the precomputed TFs from the Fs )
    F = zeros(NumOfTFs, SizF, NumReplicas);
    %
else
    %
    F = zeros(NumOfTFs, SizF);
    %
end

gPerm = randperm(size(TFs,2));       
ch = gPerm(1); 
TFindex = ch; 
model.Likelihood.TF = TFs{ch};

% fake variable
model.F = F;

% index of the TF fucntion form the training samples 
model.TFindex = TFindex; 

% Initial proposal Gaussian distribution (with diagonal covariance matrices) 
% for the kinetic parameters interection weights and the lengthscale of the GPs 
PropDist.kin = 0.05*ones(NumOfGenes,SizKin);
% interaction weigths and bias 
PropDist.W = 0.05*ones(NumOfGenes,NumOfTFs+1);

% useful ranges needed in the adaption of the 
% variances of theese proposal distribution 
qKinBelow = 0.000001; qKinAbove = 2;
qWbelow = 0.000001;   qWabove = 2;
epsilon = 0.1;

cnt = 0;
%
% do the adaption 
nextbreak = 0; 
while 1
%
%  
   [model PropDist samples accRates] = gpmtfTestGenesSample2(model, TFs, simMat, PropDist, AdaptOps);
 
   accRateKin = accRates.Kin;
   accRateW = accRates.W;
   accRateF = accRates.F;
   if AdaptOps.disp == 1
   fprintf(1,'------ ADAPTION STEP #%2d ------ \n',cnt+1); 
   fprintf(1,'Acceptance Rates for GP functions\n');
   disp(accRateF);
       
   fprintf(1,'Acceptance Rates for kinetic parameters (per gene))\n');
   disp(accRateKin);
 
   fprintf(1,'Acceptance Rates for Interaction weights (per gene)\n');
   disp(accRateW);
   fprintf(1,'Average likelihood value %15.8f',mean(samples.LogL));
   fprintf(1,'\n');
   fprintf(1,'------------------------------- \n',cnt+1);
   end
   
   
   % if you got a good acceptance rate, then stop
   if  (min(accRateKin(:))>20) & (min(accRateW(:))>20)  & (max(accRateKin(:))<35) & (max(accRateW(:))<35) 
      if nextbreak == 1
          disp('END OF ADAPTION: acceptance rates OK');
          break;
      else
          nextbreak = 1;
      end
   end
    
   cnt = cnt + 1;
   % do not allow more than 50 iterations when you adapt the proposal distribution
   if cnt == 150   
       disp('END OF ADAPTION: acceptance rates OK');
       break;
       %  
   end
   
   %%%%%%%%%%%%%%%%%%%%%%% START of ADAPT KINETICS PROPOSAL %%%%%%%%%%%%%%%%
   % adapt the proposal over the kinetic parameters (desired acceptance rate: 15-35%)
   for j=1:NumOfGenes
      if accRateKin(j) > 35
         % incease the covariance to reduce the acceptance rate
         PropDist.kin(j,:) = PropDist.kin(j,:) + epsilon*PropDist.kin(j,:);
         if PropDist.kin(j,1) > qKinAbove 
             PropDist.kin(j,:) = qKinAbove*ones(1,SizKin);
         end
      end
      if accRateKin(j) < 20
         % decrease the covariance to incease the acceptance rate
         PropDist.kin(j,:) = PropDist.kin(j,:) - epsilon*PropDist.kin(j,:);    
         if PropDist.kin(j,1) < qKinBelow 
             PropDist.kin(j,:) = qKinBelow*ones(1,SizKin);
         end
         %
      end
       %
   end
   %%%%%%%%%%%%%%%%%%%%%%% END of ADAPT KINETICS PROPOSAL %%%%%%%%%%%%%%%%
         
   %%%%%%%%%%%%%%%%%%%%%%% START of ADAPT WEIGHTS PROPOSAL %%%%%%%%%%%%%%%%
   % adapt the proposal over the interaction weights (desired acceptance rate: 15-35%)
   for j=1:NumOfGenes
      if accRateW(j) > 35
         % incease the covariance to reduce the acceptance rate
         PropDist.W(j,:) = PropDist.W(j,:) + epsilon*PropDist.W(j,:);
         if PropDist.W(j,1) > qWabove 
             PropDist.W(j,:) = qWabove*ones(1,NumOfTFs+1);
         end
      end
      if accRateW(j) < 20
         % decrease the covariance to incease the acceptance rate
         PropDist.W(j,:) = PropDist.W(j,:) - epsilon*PropDist.W(j,:);    
         if PropDist.W(j,1) < qWbelow 
             PropDist.W(j,:) = qWbelow*ones(1,NumOfTFs+1);
         end
         %
      end
       %
   end
   %%%%%%%%%%%%%%%%%%%%%%% END of ADAPT WEIGHTS PROPOSAL %%%%%%%%%%%%%%%%
   
%
%
end


                                                                                                                                         gpmtfTestGenesSample.m                                                                              0000600 0003466 0000024 00000015447 11310244603 015044  0                                                                                                    ustar   mtitsias                        games                                                                                                                                                                                                                  function [model PropDist samples accRates] = gpmtfTestGenesSample(model, PropDist, trainOps)
% Description: Draw a set of samples from the Bayesian differential
%              equation model
%
% Inputs: 
%         -- model: the structure that contains the likelihood and GP
%                    parameters as well as the priors for all these
%                    quantities
%         -- PropDist: a stucture that defines the functional form of the proposal distribution
%         -- trainOps: user defined options about the burn-in and sampling iterations
%                      and others (see demos)
%
% Outputs: model: 
%         -- model: as above. The outputed model is updated to contain the
%                   parameters values of the final MCMC iteration
%                   parameters as well as the priors
%         -- PropDist: as above. PropDist can be updated (compared to the input one) 
%                     due to the update of the kernel parameters that
%                     influence the proposal 
%         -- samples: the structure that contrains the samples 
%         -- accRates: acceptance rates 
%


BurnInIters = trainOps.Burnin; 
Iters = trainOps.T; 
StoreEvery = trainOps.StoreEvery;

Genes = model.Likelihood.Genes;
TimesG = model.Likelihood.TimesG; 
TimesF = model.Likelihood.TimesF; 
SizF = size(TimesF,2);
[NumOfGenes SizG NumOfReplicas] = size(Genes);
NumOfTFs = model.Likelihood.numTFs;

%  check if the initial condition is fixed
fixInitCond = 0;
if strcmp(model.constraints.InitialConds_value,'fixed')==1 
    fixInitCond = 1;
end

% check if the interaction weigths in the connectivity network are
% constrained to be positive
posw = 0; 
if strcmp(model.prior.weights.constraint,'positive')
    posw = 1;
end


% take the initial likelihood-kinetics parameters (defined out of this function)
LikParams = model.Likelihood;
SizKin = size(LikParams.kinetics,2);

n = SizF;
F = model.F;

% compute initial values for the log likelihood 
oldLogLik = zeros(NumOfReplicas, NumOfGenes);
% perform an evaluation of the log likelihood log p(Genes | F) 
for r=1:NumOfReplicas
  %
  % evaluate the likelihood 
  [oldLogLik(r,:) predgen] = gpmtfLogLikelihoodGene(model.Likelihood, F(:,:,r), r, 1:NumOfGenes);
  PredictedGenes(:,:,r) = predgen;  
  %     
end

% evaluation of the log prior for the kinetic parameters
lnpriorKin = ['ln',model.prior.kinetics.type,'pdf'];
TrspaceKin = model.prior.kinetics.priorSpace; 
Likkin = feval(TrspaceKin, LikParams.kinetics+eps);
oldLogPriorKin = feval(lnpriorKin, Likkin, model.prior.kinetics.a, model.prior.kinetics.b);
% evaluation of the prior for the interaction weights
lnpriorW = ['ln',model.prior.weights.type,'pdf'];
TrspaceW = model.prior.weights.priorSpace; 
LikW = feval(TrspaceW, [LikParams.W, LikParams.W0]+eps);
oldLogPriorW = feval(lnpriorW, LikW, model.prior.weights.mu, model.prior.weights.sigma2);

cnt = 0;
acceptKin = zeros(1,NumOfGenes);
acceptW = zeros(1,NumOfGenes); 
%
for it = 1:(BurnInIters + Iters) 
    %
    %
    % sample new kinetic parameters for each gene separately 
    for j=1:NumOfGenes
        KineticsNew = randn(1,SizKin).*sqrt(PropDist.kin(j,:)) + log(LikParams.kinetics(j,:));
        KineticsNew(KineticsNew<-10) =-10; 
        KineticsNew(KineticsNew>10) = 10;
        KineticsNew = exp(KineticsNew); 
        %
        % set the initial condition to be zero
        % if you know that it should be zero
        if model.constraints.InitialConds(j) == 0
        KineticsNew(4) = KineticsNew(1)/KineticsNew(2); 
        end
  
        LikParams1 = LikParams;
        LikParams1.kinetics(j,:)=KineticsNew; 
        newLogLik = [];
        for r=1:NumOfReplicas
          % call the function only with j gene expressions  
          newLogLik(r) = gpmtfLogLikelihoodGene(LikParams1, F(:,:,r), r, j);
        end
        %        
        Likkin = feval(TrspaceKin, KineticsNew);
        LogPriorKinNew = feval(lnpriorKin, Likkin, model.prior.kinetics.a, model.prior.kinetics.b);
        % Metropolis-Hastings to accept-reject the proposal
        oldP = sum(oldLogLik(:,j),1) + sum(oldLogPriorKin(j,:),2);
        newP = sum(newLogLik(:))+ sum(LogPriorKinNew(:)); 
        [accept, uprob] = metropolisHastings(newP, oldP, 0, 0);
        if accept == 1
           LikParams.kinetics(j,:) = KineticsNew;
           oldLogLik(:,j) = newLogLik(:); 
           oldLogPriorKin(j,:) = LogPriorKinNew;
        end
        %
        if (it > BurnInIters) 
           acceptKin(j) = acceptKin(j) + accept;
        end
        %
    end
    %
   
    % sample the interaction weights 
    for j=1:NumOfGenes
        % 
        %  
        if posw == 1
            Wnew = randn(1,NumOfTFs+1).*sqrt(PropDist.W(j,:)) + log([LikParams.W(j,:), LikParams.W0(j)]+eps);     
            Wnew = exp(Wnew);
        else
            Wnew = randn(1,NumOfTFs+1).*sqrt(PropDist.W(j,:)) + [LikParams.W(j,:), LikParams.W0(j)];     
        end
        Wnew(1:NumOfTFs) = Wnew(1:NumOfTFs).*model.constraints.W(j,:);
       
        LikParams1 = LikParams;
        LikParams1.W(j,:) = Wnew(1:NumOfTFs);
        LikParams1.W0(j)=Wnew(end);
      
        newLogLik = [];
        for r=1:NumOfReplicas
          % call the function only with j gene expressions  
          newLogLik(r) = gpmtfLogLikelihoodGene(LikParams1, F(:,:,r), r, j);
        end
        
        LikW = feval(TrspaceW, Wnew+eps);
        LogPriorWnew = feval(lnpriorW, LikW, model.prior.weights.mu, model.prior.weights.sigma2);
        % Metropolis-Hastings to accept-reject the proposal
        oldP = sum(oldLogLik(:,j),1) + sum(oldLogPriorW(j,:),2);
        newP = sum(newLogLik(:)) + sum(LogPriorWnew(:)); 
        %
        [accept, uprob] = metropolisHastings(newP, oldP, 0, 0);
        if accept == 1
           LikParams.W(j,:) = Wnew(1:NumOfTFs);
           LikParams.W0(j) = Wnew(end);
           oldLogLik(:,j) = newLogLik(:); 
           oldLogPriorW(j,:) = LogPriorWnew;
        end
        %
        if (it > BurnInIters) 
           acceptW(j) = acceptW(j) + accept;
        end
        %
    end
    %
    % keep samples after burn in
    if (it > BurnInIters)  & (mod(it,StoreEvery) == 0)
        %
        cnt = cnt + 1;
        samples.predGenes{cnt} = PredictedGenes;
        samples.kinetics(:,:,cnt) = LikParams.kinetics;
        samples.Weights(:,:,cnt) = LikParams.W;
        samples.Weights0(:,cnt) = LikParams.W0;
        samples.LogL(cnt) = sum(oldLogLik(:));
        %
    end
    %
    %        
end

% Before you return store the final state of the Markov chain to 
% the model structure
model.Likelihood.kinetics = LikParams.kinetics;
model.Likelihood.W = LikParams.W;
model.Likelihood.W0 = LikParams.W0;

accRates.Kin = (acceptKin/Iters)*100;
accRates.W = (acceptW/Iters)*100;

                                                                                                                                                                                                                         gpmtfTestGenesSample2.m                                                                             0000600 0003466 0000024 00000041420 11320643332 015117  0                                                                                                    ustar   mtitsias                        games                                                                                                                                                                                                                  function [model PropDist samples accRates] = gpmtfTestGenesSample2(model, TFs, simMat, PropDist, trainOps)
% Description: Draw a set of samples from the Bayesian differential
%              equation model
%
% Inputs: 
%         -- model: the structure that contains the likelihood and GP
%                    parameters as well as the priors for all these
%                    quantities
%         -- PropDist: a stucture that defines the functional form of the proposal distribution
%         -- trainOps: user defined options about the burn-in and sampling iterations
%                      and others (see demos)
%
% Outputs: model: 
%         -- model: as above. The outputed model is updated to contain the
%                   parameters values of the final MCMC iteration
%                   parameters as well as the priors
%         -- PropDist: as above. PropDist can be updated (compared to the input one) 
%                     due to the update of the kernel parameters that
%                     influence the proposal 
%         -- samples: the structure that contrains the samples 
%         -- accRates: acceptance rates 
%


BurnInIters = trainOps.Burnin; 
Iters = trainOps.T; 
StoreEvery = trainOps.StoreEvery;

Genes = model.Likelihood.Genes;
TimesG = model.Likelihood.TimesG; 
TimesF = model.Likelihood.TimesF; 
SizF = size(TimesF,2);
[NumOfGenes SizG NumOfReplicas] = size(Genes);
NumOfTFs = model.Likelihood.numTFs;

num_stored = floor(Iters/StoreEvery);
samples.TFindex = zeros(1, num_stored);
%samples.predGenes = zeros(NumOfReplicas, SizF, num_stored);
samples.kinetics = zeros(4, num_stored);
samples.Weights = zeros(NumOfTFs, num_stored);
samples.Weights0 = zeros(1, num_stored);
samples.LogL = zeros(1, num_stored);

%  check if the initial condition is fixed
fixInitCond = 0;
if strcmp(model.constraints.InitialConds_value,'fixed')==1 
    fixInitCond = 1;
end

% check if the observation noise is known/fixed
fixsigma2 = 0;
if strcmp(model.constraints.sigmas,'fixed')
    fixsigma2 = 1;
end

% check if the interaction weigths in the connectivity network are
% constrained to be positive
posw = 0; 
if strcmp(model.prior.weights.constraint,'positive')
    posw = 1;
end


% take the initial likelihood-kinetics parameters (defined out of this function)
LikParams = model.Likelihood;
SizKin = size(LikParams.kinetics,2);

n = SizF;
F = model.F;
TFindex = model.TFindex;


% compute initial values for the log likelihood 
oldLogLik = zeros(NumOfReplicas, NumOfGenes);
if strcmp(model.constraints.replicas,'free')
   for r=1:NumOfReplicas
   %
   % evaluate the likelihood 
   [oldLogLik(r,:) predgen] = gpmtfLogLikelihoodGene(model.Likelihood, F(:,:,r), r, 1:NumOfGenes);
   PredictedGenes(r,:) = predgen;
   %     
   end
else
   %
   % evaluate the likelihood for the first replica
   [oldLogLik(1,:) predgen] = gpmtfLogLikelihoodGene(model.Likelihood, F, 1, 1:NumOfGenes);
   PredictedGenes = predgen;
   % compute fast the additional likelihood when you have observations for the TF genes  
   %
   % the predicted genes are the same for the remaining coupled replicas
   for r=2:NumOfReplicas
      % compute fast the likelihood terms for the remaining replicas  
      oldLogLik(r,:) = remainRepsLikelihood(LikParams,  PredictedGenes, r, 1:NumOfGenes);
      % compute fast the additional likelihood when you have observations
      % for the TF genes  
   end
   %              
   %
end

% evaluation of the log prior for the kinetic parameters
lnpriorKin = ['ln',model.prior.kinetics.type,'pdf'];
TrspaceKin = model.prior.kinetics.priorSpace; 
Likkin = feval(TrspaceKin, LikParams.kinetics+eps);
oldLogPriorKin = feval(lnpriorKin, Likkin, model.prior.kinetics);

% evaluation of the prior for the interaction bias
lnpriorW0 = ['ln',model.prior.weight0.type,'pdf'];
TrspaceW0 = model.prior.weight0.priorSpace; 
LikW0 = feval(TrspaceW0, LikParams.W0);
oldLogPriorW0 = feval(lnpriorW0, LikW0, model.prior.weight0);

% evaluation of the prior for the interaction weights
lnpriorW = ['ln',model.prior.weights.type,'pdf'];
TrspaceW = model.prior.weights.priorSpace; 
LikW = feval(TrspaceW, LikParams.W);
oldLogPriorW = feval(lnpriorW, LikW, model.prior.weights);

cnt = 0;

%if strcmp(model.constraints.replicas,'free')
%   acceptF = zeros(NumOfTFs,NumOfReplicas);
%else
%   acceptF = zeros(NumOfTFs,1);   
%end
acceptF = 0; 

acceptKin = zeros(1,NumOfGenes);
acceptW = zeros(1,NumOfGenes);
numSamples = size(TFs,2);
%
for it = 1:(BurnInIters + Iters) 
    %
    % choose one sample for the TFS from the training set  
    
    % choose a training sample 
    ch = round(rand*numSamples) + 1; 
    ch(ch>numSamples) = numSamples;
    
    LikParams1 = LikParams;
    % store the TF in the LikeParams to save computations 
    LikParams1.TF = TFs{ch};
    
    
    newLogProp = 0;  % forward Hastings Q(s_t+1 | s_t)
    oldLogProp = 0;  % backward Hastings Q(s_t| s_t+1) 
    
    newLogLik = zeros(NumOfReplicas, NumOfGenes);
    predgen = zeros(NumOfReplicas, model.Likelihood.sizTime);   
    if strcmp(model.constraints.replicas,'free') 
       for r=1:NumOfReplicas
       %  
       
       %for j=1:NumOfTFs
       %     
           
           %%%%%%%%%%%%%%%%%%  bit of code to be changed 
           % choose randomly among all samples 
           %ch = round(rand*numSamples) + 1; 
           %ch(ch>numSamples) = numSamples;
           
           %newLogProp = 0;  % forward Hastings Q(s_t+1 | s_t)
           %oldLogProp = 0;  % backward Hastings Q(s_t| s_t+1) 
           
           %% draw from the geometric distribution
           %kk = geornd(model.geo) + 1;
           %% propose the kkth nearest neighbor of the next TF
           %ch = simMat{j,r}(TFindex(j,r),kk);
           %
           % 
           %% forward Hastings Q(s_t+1 | s_t)
           %newLogProp = log(geopdf(kk-1,model.geo));
           %
           %% backward Hastings Q(s_t| s_t+1) 
           %bkk = find(simMat{j,r}(ch,:)==TFindex(j,r));
           %oldLogProp = log(geopdf(bkk-1,model.geo)); 
          
           %%%%%%%%%%%%%%%%%% end of bit of code to be changed 
           
           %LikParams1 = LikParams;
           % store the TF in the LikeParams to save computations 
           %LikParams1.TF(j,:,r) = TFs{ch}(j,:,r);
           
           % perform an evaluation of the likelihood p(Genes | F) 
           [newLogLik(r,:) predgen(r,:)] = gpmtfLogLikelihoodGene(LikParams1, F(:,:,r), r, 1:NumOfGenes);
       %end % num TFs loop
           
       end % num Replicas loop
       %
       % Metropolis-Hastings to accept-reject the proposal
       [accept, uprob] = metropolisHastings(sum(newLogLik(:)),sum(oldLogLik(:)), newLogProp, oldLogProp);
           
       if (it > BurnInIters)
           acceptF = acceptF + accept; 
       end
    
       % update protein F
       if accept == 1
           %LikParams.TF(j,:,r) = TFs{ch}(j,:,r);
           LikParams.TF = TFs{ch};
           %TFindex(1:NumOfTFs,1:NumOfReplicas) = ch;
           TFindex = ch;
           PredictedGenes = predgen;
           oldLogLik = newLogLik;
        %   
       end
    %
    else
         % Replicas are coupled
         for j=1:NumOfTFs
             
           %%%%%%%%%%%%%%%%%%  bit of code to be changed 
           % choose randomly among all samples 
           gPerm = randperm(size(TFs,2));
           ch = gPerm(1); 
           
           newLogProp = 0;
           oldLogProp = 0; 
           %%%%%%%%%%%%%%%%%% end of bit of code to be changed 
           
           LikParams1 = LikParams;
           % store the TF in the LikeParams to save computations 
           LikParams1.TF(j,:) = TFs{ch}(j,:);
             
           newLogLik = zeros(NumOfReplicas,NumOfGenes);
           % perform an evaluation of the likelihood p(Genes | F)      
           [newLogLik(1,:) predgen] = gpmtfLogLikelihoodGene(LikParams1, F, 1, 1:NumOfGenes);
              
           % computed faster the remaining likelihood terms  
           for r=2:NumOfReplicas
               newLogLik(r,:) = remainRepsLikelihood(LikParams1,  predgen, r, 1:NumOfGenes);
           end                 
           % Metropolis-Hastings to accept-reject the proposal
           [accept, uprob] = metropolisHastings(sum(newLogLik(:)),sum(oldLogLik(:)), newLogProp, oldLogProp);
          
           %
           if (it > BurnInIters)
               acceptF(j) = acceptF(j) + accept; 
           end
    
           % update protein F
           if accept == 1
               LikParams.TF(j,:) = TFs{ch}(j,:);
               TFindex(j) = ch;  
               PredictedGenes = predgen;
               oldLogLik = newLogLik;
            %   
           end
         end % num TFs loop
        %
    end % if end
    
    
    %
    % sample new kinetic parameters for each gene separately 
    for j=1:NumOfGenes
        KineticsNew = randn(1,SizKin).*sqrt(PropDist.kin(j,:)) + log(LikParams.kinetics(j,:));
        KineticsNew(KineticsNew<-10) =-10; 
        KineticsNew(KineticsNew>10) = 10;
        KineticsNew = exp(KineticsNew); 
        %
        % set the initial condition to be zero
        % if you know that it should be zero
        if model.constraints.InitialConds(j) == 0
        KineticsNew(4) = KineticsNew(1)/KineticsNew(2); 
        end
  
        LikParams1 = LikParams;
        LikParams1.kinetics(j,:)=KineticsNew; 
        newLogLik = zeros(1,NumOfReplicas);
        if strcmp(model.constraints.replicas,'free')
           %
           for r=1:NumOfReplicas
              % call the function only with j gene expressions  
              newLogLik(r) = gpmtfLogLikelihoodGene(LikParams1, F(:,:,r), r, j);
           end
        else
            % call the function only for the first replica
            [newLogLik(1), predgen] = gpmtfLogLikelihoodGene(LikParams1, F, 1, j);
            % computed faster the remaining likelihood terms  
            for r=2:NumOfReplicas
                newLogLik(r) = remainRepsLikelihood(LikParams1, predgen, r, j);
            end 
           %
        end
       
        %        
        Likkin = feval(TrspaceKin, KineticsNew);
        LogPriorKinNew = feval(lnpriorKin, Likkin, model.prior.kinetics);
        % Metropolis-Hastings to accept-reject the proposal
        oldP = sum(oldLogLik(:,j),1) + sum(oldLogPriorKin(j,:),2);
        newP = sum(newLogLik(:))+ sum(LogPriorKinNew(:)); 
        [accept, uprob] = metropolisHastings(newP, oldP, 0, 0);
        if accept == 1
           LikParams.kinetics(j,:) = KineticsNew;
           oldLogLik(:,j) = newLogLik(:); 
           oldLogPriorKin(j,:) = LogPriorKinNew;
        end
        %
        if (it > BurnInIters) 
           acceptKin(j) = acceptKin(j) + accept;
        end
        %
    end
    %
   
    % sample the interaction weights 
    for j=1:NumOfGenes
        % 
        %   
        newLogProp = 0;  % forward Hastings Q(w_t+1 | w_t)
        oldLogProp = 0;  % backward Hastings Q(w_t| w_t+1)  
        if posw == 1
            %
            % sample from a truncated Gaussian using rejection
            % sampling 
            while 1
                trW = randn(1,NumOfTFs).*sqrt(PropDist.W(j,1:NumOfTFs)) + LikParams.W(j,:);  
                if min(trW) >= 0
                    break; 
                end
            end 
            Wnew(1:NumOfTFs) = trW;
            % sample also the bias which is allowed to be negative
            Wnew0 = randn.*sqrt(PropDist.W(j,NumOfTFs+1)) +  LikParams.W0(j);
            Wnew = [trW, Wnew0]; 
            
            trprior.sigma2 = PropDist.W(j,1:NumOfTFs); 
            trprior.mu = LikParams.W(j,:); 
        
            newLogProp = sum(lntruncNormalpdf(trW, trprior)); 
            
            trprior.mu = trW; 
            oldLogProp  = sum(lntruncNormalpdf(LikParams.W(j,:), trprior)); 
            
            %
        else
            Wnew = randn(1,NumOfTFs+1).*sqrt(PropDist.W(j,:)) + [LikParams.W(j,:), LikParams.W0(j)];     
        end
        Wnew(1:NumOfTFs) = Wnew(1:NumOfTFs).*model.constraints.W(j,:);
        
        if model.constraints.W0(j) == 0
           Wnew(NumOfTFs + 1) = 0;
        end
       
        LikParams1 = LikParams;
        LikParams1.W(j,:) = Wnew(1:NumOfTFs);
        LikParams1.W0(j)=Wnew(end);
      
        newLogLik = zeros(1,NumOfReplicas);
        if strcmp(model.constraints.replicas,'free')
           for r=1:NumOfReplicas
               % call the function only with j gene expressions  
               newLogLik(r) = gpmtfLogLikelihoodGene(LikParams1, F(:,:,r), r, j);
           end
        else
           % call the function only for the first replica
           [newLogLik(1), predgen] = gpmtfLogLikelihoodGene(LikParams1, F, 1, j);
           % computed faster the remaining likelihood terms  
           for r=2:NumOfReplicas
               newLogLik(r) = remainRepsLikelihood(LikParams1, predgen, r, j);
           end 
        end
        
        
        % evaluation of the prior for the interaction bias 
        LikW0 = feval(TrspaceW0, LikParams.W0);
        LogPriorWnew0 = feval(lnpriorW0, Wnew(end), model.prior.weight0);
        % >>>  interaction weights
        LikW = feval(TrspaceW, Wnew(1:NumOfTFs));
        LogPriorWnew = feval(lnpriorW, LikW, model.prior.weights);
        
        % Metropolis-Hastings to accept-reject the proposal
        oldP = sum(oldLogLik(:,j),1) + sum(oldLogPriorW(j,:),2) + oldLogPriorW0(j);
        newP = sum(newLogLik(:)) + sum(LogPriorWnew(:)) + LogPriorWnew0; 
         
        [accept, uprob] = metropolisHastings(newP, oldP, newLogProp, oldLogProp);
        if accept == 1
           LikParams.W(j,:) = Wnew(1:NumOfTFs);
           LikParams.W0(j) = Wnew(end);
           oldLogLik(:,j) = newLogLik(:); 
           oldLogPriorW(j,:) = LogPriorWnew;
           oldLogPriorW0(j) = LogPriorWnew0;
        end
        %
        if (it > BurnInIters) 
           acceptW(j) = acceptW(j) + accept;
        end
        %
    end
    
    % sample the noise of the likelihood when is free parameter.
    % The posterior is gamma so this step involves exact simualtion from
    % the gamma distribution
    if fixsigma2 == 0
       sumSquerrors1 = oldLogLik + 0.5*SizG*log(2*pi*repmat(LikParams.sigmas(:,1,1)',NumOfReplicas,1));
       sumSquerrors1 = -2*repmat(LikParams.sigmas(:,1,1)',NumOfReplicas,1).*sumSquerrors1;
       sumSquerrors = sum(sumSquerrors1,1);
                
       anew = model.prior.invsigma2.a + 0.5*NumOfReplicas*SizG;
       bnew = model.prior.invsigma2.b + 0.5*sumSquerrors;
       newinvsigma2 = gamrnd(anew,1./bnew);
       Nnewsigma2 = 1./newinvsigma2;
       
       %[Nnewsigma2 LikParams.sigmas(1,1,1)]
       
       LikParams.sigmas = repmat(Nnewsigma2(:),[1 SizG NumOfReplicas]);
       % 
       okk = repmat(LikParams.sigmas(:,1,1)',NumOfReplicas,1); 
       
       oldLogLik = - 0.5*SizG*log(2*pi*okk) - (0.5./okk).*sumSquerrors1;
       
       % 
       %
    end
    
    
    %
    % keep samples after burn in
    if (it > BurnInIters)  & (mod(it,StoreEvery) == 0)
        %
        cnt = cnt + 1;
        samples.TFindex(cnt) = TFindex;
        %samples.predGenes(:,:,cnt) = PredictedGenes;
        samples.kinetics(:,cnt) = LikParams.kinetics;
        samples.Weights(:,cnt) = LikParams.W;
        samples.Weights0(:,cnt) = LikParams.W0;
        if fixsigma2 == 0
           samples.sigmas(cnt) = LikParams.sigmas(1,1,1);
        end
        samples.LogL(cnt) = sum(oldLogLik(:));
        %
    end
    %
    %        
end

% Before you return store the final state of the Markov chain to 
% the model structure
model.Likelihood.kinetics = LikParams.kinetics;
model.Likelihood.W = LikParams.W;
model.Likelihood.W0 = LikParams.W0;

model.TFindex = TFindex; 

accRates.F = (acceptF/Iters)*100;
accRates.Kin = (acceptKin/Iters)*100;
accRates.W = (acceptW/Iters)*100;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function loglikval = remainRepsLikelihood(LikParams,  PredictedGenes, r, Gindex)
%
    Gns = LikParams.Genes(:,:,r);
    loglikval = - 0.5*sum(log(2*pi*LikParams.sigmas(Gindex,:,r)),2)....
                - 0.5*sum(((Gns(Gindex,:) - PredictedGenes(:,LikParams.comInds)).^2)./LikParams.sigmas(Gindex,:,r),2);
    loglikval = loglikval';
% 
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function loglikvalTF = remainRepsLikelihoodTF(LikParams, PredictedGenesTF, r, TFindex)
%    
  
  GnsTF = LikParams.GenesTF(:,:,r);
  loglikvalTF = - 0.5*sum(log(2*pi*LikParams.sigmasTF(TFindex,:,r)),2)....
                       - 0.5*sum(((GnsTF(TFindex,:) - PredictedGenesTF(:,LikParams.comIndsTF)).^2)./LikParams.sigmasTF(TFindex,:,r),2);
  loglikvalTF = loglikvalTF';
%                                                                                                                                                                                                                                                  gpmtfTestPlot.m                                                                                     0000600 0003466 0000024 00000011701 11320465160 013550  0                                                                                                    ustar   mtitsias                        games                                                                                                                                                                                                                  function gpmtfTestPlot(model, Genes, GeneVars, fbgn, samples, TFs, demdata, printResults, Grtruth)

dirr = '/usr/local/michalis/mlprojects/gpsamp/tex/diagrams/';
dirrhtml = '/usr/local/michalis/mlprojects/gpsamp/html/';

%Genes = model.Likelihood.Genes;
gg = 1;

%if strcmp(model.constraints.sigmas,'fixed')
%    GeneVars = model.Likelihood.sigmas;
%end
 
TimesG = model.Likelihood.TimesG; 
TimesF = model.Likelihood.TimesF;
if isfield(model.Likelihood,'GenesTF')
    GenesTF = model.Likelihood.GenesTF;
    if strcmp(model.constraints.sigmasTF,'fixed')
        GeneTFVars = model.Likelihood.sigmasTF;
    end
end

TimesFF = TimesF(model.Likelihood.startTime:end);
ok = date;
fileName = [demdata 'Test' 'MCMC' ok model.Likelihood.singleAct model.Likelihood.jointAct char(fbgn)]; 

NumOfTFs = model.Likelihood.numTFs;


NumOfGenes = model.Likelihood.numGenes;
order = 1:NumOfGenes;
NumOfReplicas = model.Likelihood.numReplicas;
NumOfSamples = size(samples.LogL,2);
TimesF = TimesF(:);
SizF = size(TimesF,1);

% plot predicted gene expressions 
for r=1:NumOfReplicas
  % 
  GG = zeros(NumOfSamples,model.Likelihood.sizTime);    
  for t=1:NumOfSamples
      LikParams = model.Likelihood; 
      LikParams.kinetics = samples.kinetics(:,t)';
      LikParams.W = samples.Weights(:,t)';
      LikParams.W0 = samples.Weights0(t);
      LikParams.TF = TFs{samples.TFindex(t)};  
      predgen = gpmtfComputeGeneODE(LikParams, zeros(NumOfTFs, SizF), r, 1);
      GG(t,:) = predgen;
  end
     
  mu = mean(GG)';
  stds = sqrt(var(GG))';
    
  TF = TimesFF'; % TimesFF(1:2:end)';
  figure
  plot(TF,mu,'b','lineWidth',r);
  hold on;
  fillColor = [0.7 0.7 0.7];
  %fillColor = [0.8 0.8 0.8];  % for the paper
  fill([TF; TF(end:-1:1)], [mu; mu(end:-1:1)]...
        + 2*[stds; -stds(end:-1:1)], fillColor,'EdgeColor',fillColor);
  plot(TF,mu,'b','lineWidth',3);
   
  plot(TimesG,Genes(1,:,r),'rx','markersize', 14','lineWidth', 2);
  if strcmp(model.constraints.sigmas,'fixed')
     errorbar(TimesG,  Genes(1,:,r), 2*sqrt(GeneVars(1,:,r)), 'rx','lineWidth', 1.5);
  end
  axis([min(TimesG(:))-0.1 max(TimesG(:))+0.1 0.95*min(min(Genes(1,:,r))) 1.05*max(max(Genes(1,:,r)))]);
     
  if printResults
     print('-depsc', [dirr fileName 'Replica' num2str(r) 'GeneExp' num2str(gg)]);
     print('-dpng', [dirrhtml fileName 'Replica' num2str(r) 'GeneExp' num2str(gg)]);
  end
  titlestring = 'Expressions: ';
  titlestring = [titlestring, num2str(r)]; 
  titlestring = [titlestring, ' replica, '];
  titlestring = [titlestring, num2str(gg)];
  titlestring = [titlestring, ' gene'];
  title(titlestring,'fontsize', 20);
  %
  %
end

figure;
h = hist(samples.kinetics(1,:));
hist(samples.kinetics(1,:));
title('Basal rates','fontsize', 20);
% ground truth 
if nargin == 9
    %
    hold on;
    plot([Grtruth.kinetics(1) Grtruth.kinetics(1)], [0 max(h)], 'LineWidth', 5,'Color', 'r');
    %
end
    

if printResults
      print('-depsc', [dirr fileName 'Basal']);
      print('-dpng', [dirrhtml fileName 'Basal']);
end
figure;   
h = hist(samples.kinetics(3,:));
hist(samples.kinetics(3,:));
title('Sensitivities','fontsize', 20);
% ground truth 
if nargin == 9
    %
    hold on;
    plot([Grtruth.kinetics(3) Grtruth.kinetics(3)], [0 max(h)], 'LineWidth', 5,'Color', 'r');
    %
end



if printResults
      print('-depsc', [dirr fileName 'Sensitivity']);
      print('-dpng', [dirrhtml fileName 'Sensitivity']);
end
figure;
h = hist(samples.kinetics(2,:));
hist(samples.kinetics(2,:));
title('Decays','fontsize', 20);
% ground truth 
if nargin == 9
    %
    hold on;
    plot([Grtruth.kinetics(2) Grtruth.kinetics(2)], [0 max(h)], 'LineWidth', 5,'Color', 'r');
    %
end


if printResults
      print('-depsc', [dirr fileName 'Decay']);
      print('-dpng', [dirrhtml fileName 'Decay']);
end
figure;
h = hist(samples.kinetics(4,:));
hist(samples.kinetics(4,:));
title('Initial conditions','fontsize', 20);
% ground truth 
if nargin == 9
    %
    hold on;
    plot([Grtruth.kinetics(4) Grtruth.kinetics(4)], [0 max(h)], 'LineWidth', 5,'Color', 'r');
    %
end

if printResults
      print('-depsc', [dirr fileName 'InitCond']);
      print('-dpng', [dirrhtml fileName 'InitCond']);
end


for j=1:NumOfTFs
W1 = samples.Weights(j,:);
figure;
h = hist(W1);
hist(W1);
% ground truth 
if nargin == 9
    %
    hold on;
    plot([Grtruth.W(j) Grtruth.W(j)], [0 max(h)], 'LineWidth', 5,'Color', 'r');
    %
end

if printResults
     print('-depsc', [dirr fileName 'IntWeights' 'TF' num2str(j)]);
     print('-dpng', [dirrhtml fileName 'IntWeights' 'TF' num2str(j)]);
end
titlestring = 'Interaction weights: '; 
titlestring = [titlestring, num2str(j)];
titlestring = [titlestring, ' TF'];
title(titlestring,'fontsize', 20);
end

W0 = samples.Weights0';

figure;
h = hist(W0);
hist(W0);
if nargin == 9
    %
    hold on;
    plot([Grtruth.W0 Grtruth.W0], [0 max(h)], 'LineWidth', 5,'Color', 'r');
    %
end

if printResults
      print('-depsc', [dirr fileName 'IntBias']);
      print('-dpng', [dirrhtml fileName 'IntBias']);
end
title('Interaction biases','fontsize', 20);
                                                               gpsampControlAdapt.m                                                                                0000700 0003466 0000024 00000015512 11275534152 014552  0                                                                                                    ustar   mtitsias                        games                                                                                                                                                                                                                  function [model PropDist samples accRates] = gpsampControlAdapt(model, ops)
% adaptive MCMc to specify the number and location of control points 
% as well as proposal distribution for kernel hyperparameters 
%
%

BurnInIters = ops.Burnin; 
Iters = ops.T; 
% initial number of control variables 
M = ops.initialNumContrPnts;
if M < 3, 
    M = 3; 
end
Init = ops.initialNumContrPnts;
IncM = ops.incrNumContrBy;

[n D] = size(model.X);
perm = randperm(n);
Xu =[];
X = model.X;


% initial input location of the control variables
if strcmp(model.constraints.kernHyper, 'fixed') 
 %
 % keep adding control points till the trace falls under the 10% 
 model.K = kernCompute(model.GP, X); 
   while 1 
      if isempty(Xu) 
         Xu = X(perm(1:M),:);
      else
         Xu = [Xu; X(perm(M+1:M+IncM),:)];
      end
      M = size(Xu,1);
      [Xu f] = minimize(Xu(:), 'trace_CondCov', 200, X, 0.5*model.GP.logtheta);
      Xu = reshape(Xu,M,D);
      if f(end) < (0.1*sum(diag(model.K)))
          break;
      end
   end
   %     
else
% 
  Xu = X(perm(1:M),:);
  [Xu f] = minimize(Xu(:), 'trace_CondCov', 200, X, 0.5*model.GP.logtheta); 
  Xu = reshape(Xu,M,D);
%
end

M = size(Xu,1); 
Fu = zeros(1, M);

% proposal distribution for F
U = n+1:n+M;  
PropDist.qF.m = zeros(n+M,1);
PropDist.qF.K = kernCompute(model.GP, [X; Xu]); 
L=jitterChol(PropDist.qF.K)';
PropDist.qF.invL = L\eye(n+M); 
PropDist.qF.LogDetK = 2*sum(log(diag(L)));      
% compute the conditional GP prior given the control variables
[cmuMinus, cSigma, KInvKu] = gaussianFastConditional(PropDist.qF.m', PropDist.qF.K, 1:n, U);
[L,er]=jitterChol(cSigma);
if er>0, L = real(sqrtm(cSigma)); end
cmu = cmuMinus + Fu*KInvKu;
F = gaussianFastSample(1, cmu, L);
model.F = F;
PropDist.qF.cmuMinus = cmuMinus;
PropDist.qF.cSigma = cSigma;
PropDist.qF.KInvKu = KInvKu;
PropDist.qF.L = L;
% compute all the conditional variances for the control Points
for i=1:M
% 
   G = [1:i-1, i+1:M];
   [alpha(i), ku(i), KInvK(i,:)] = gaussianFastConditional(PropDist.qF.m(U)', PropDist.qF.K(U,U), i, G);
%
end 
%
PropDist.qF.alpha = alpha;
PropDist.qF.ku = ku;
PropDist.qF.KInvK = KInvK;

% proposal for the kernel hyperparameters
PropDist.kern = 0.2*(1/model.prior.kernParams.b)*ones(1,model.GP.nParams);
if model.Likelihood.nParams > 0
  PropDist.lik = 0.2*(1/model.prior.likParams.b)*ones(1,model.Likelihood.nParams);
end

nextbreak = 0;
epsilon = 0.1;
cnt = 0;
model.F = F;
model.Fu = Fu;
model.Xu = Xu;

% do the adaption  
while 1
    %
    model.Xu = Xu;
    [model PropDist samples accRates] = gpsampControlTrain(model, PropDist, ops);
    
    accRateF = accRates.F;
    accRateKern = accRates.kern;
    accRateLik = accRates.lik;
 
    fprintf(1,'------ ADAPTION STEP #%2d, Number of Control Points %2d ------ \n',cnt+1,M); 
    if ops.disp == 1
       fprintf(1,'Acceptance Rates for GP function\n');
       fprintf(1,' (Rates per control point) \n');       
       disp(accRateF);    
       if ~strcmp(model.constraints.kernHyper, 'fixed')
           fprintf(1,'Acceptance Rates for kernel hyperparameters\n');
           disp(accRateKern);
       end
       if ~strcmp(model.constraints.likHyper, 'fixed') & (model.Likelihood.nParams > 0)
           fprintf(1,'Acceptance Rates for likelihood hyperparameters\n');
           disp(accRateLik);
       end
       fprintf(1,'Average likelihood value %15.8f\n',mean(samples.LogL));
    end

    % if you got a descent acceptance rate, then stop
    if ((min(accRateF(:)) > ((0.2/M)*100)) &  (accRateKern>15) & (accRateLik>15) )
        if nextbreak == 1
           disp('END OF ADAPTION: acceptance rates OK');
           break; 
        else
           nextbreak = 1;
        end
    end
    
    cnt = cnt + 1;
    % do not allow more than 80 iterations when you adapt the proposal distribution
    if cnt == 80
        warning('END OF ADAPTION: acceptance rates were not all OK');
        break;
    end

    %  increase the number of control points (add controls in every second iteration)
    if (mod(cnt,2)==0) & (min(accRateF(:)) < ((0.2/M)*100)) 
       %
       %
       Xu = [Xu; X(perm(M+1:M+IncM),:)]; 
       M = M + IncM;
       
       % optimize over control input location  by minimizing the trace term
       [Xu f] = minimize(Xu(:), 'trace_CondCov', 200, X, 0.5*model.GP.logtheta);
       Xu = reshape(Xu,M,D);
       if D == 1
           Xu = sort(Xu);
       end
       model.Xu = Xu;
       
       % initialize the new control variables given the current F
       model.Fu = zeros(1,M);
       U = n+1:n+M;
       % update proposal distribution 
       PropDist.qF.m = zeros(n+M,1);
       PropDist.qF.K = kernCompute(model.GP, [X; Xu]);     
       L=jitterChol(PropDist.qF.K)';
       PropDist.qF.invL = L\eye(n+M); 
       PropDist.qF.LogDetK = 2*sum(log(diag(L)));
      
       [cmuMinus, cSigma, KInvKu] = gaussianFastConditional(PropDist.qF.m', PropDist.qF.K, U, 1:n);
       [L,er]=jitterChol(cSigma);
       if er>0, L = real(sqrtm(cSigma)); end
       cmu = cmuMinus + model.F*KInvKu;
       model.Fu = gaussianFastSample(1, cmu, L);
 
       % compute the conditional GP prior given the control variables
       [cmuMinus, cSigma, KInvKu] = gaussianFastConditional(PropDist.qF.m', PropDist.qF.K, 1:n, U);
       [L,er]=jitterChol(cSigma);
       if er>0, L = real(sqrtm(cSigma)); end
       PropDist.qF.cmuMinus = cmuMinus; 
       PropDist.qF.cSigma = cSigma;
       PropDist.qF.KInvKu = KInvKu;
       PropDist.qF.L = L;
       clear alpha ku KInvK;
       for i=1:M
       %  
         G = [1:i-1, i+1:M];  
         [alpha(i), ku(i), KInvK(i,:)] = gaussianFastConditional(PropDist.qF.m(U)', PropDist.qF.K(U,U), i, G);
       %
       end
       PropDist.qF.alpha = alpha;
       PropDist.qF.ku = ku;
       PropDist.qF.KInvK = KInvK; 
    end
   
    
    % adapt likelihood hyperparameters (everything apart from F)
    if ~strcmp(model.constraints.likHyper, 'fixed') & (model.Likelihood.nParams > 0)
       if accRateLik > 35
        % incease the covariance to reduce the acceptance rate
        PropDist.lik = PropDist.lik + epsilon*PropDist.lik;
       end
       if accRateLik < 15
        % decrease the covariance to incease the acceptance rate
        PropDist.lik = PropDist.lik - epsilon*PropDist.lik;    
        %
       end
    end
    
    
    % adapt kernel hyperparameters proposal 
    if ~strcmp(model.constraints.kernHyper, 'fixed') 
       if accRateKern > 35
        % incease the covariance to reduce the acceptance rate
        PropDist.kern = PropDist.kern + epsilon*PropDist.kern;
       end
       if accRateKern < 15
        % decrease the covariance to incease the acceptance rate
        PropDist.kern = PropDist.kern - epsilon*PropDist.kern;    
        %
       end
    end
    %
    %
end
                                                                                                                                                                                      gpsampControlTrain.m                                                                                0000700 0003466 0000024 00000021607 11275533400 014573  0                                                                                                    ustar   mtitsias                        games                                                                                                                                                                                                                  function [model PropDist samples accRates] = gpsampControlTrain(model, PropDist, trainOps)
% Inputs: 
%         -- model: the structure that contains the likelihood and GP
%                    parameters as well as the priors for all these
%                    quantities
%         -- PropDist: a stucture that defines the functional form of the proposal distribution
%         -- trainOps: user defined options about the burn-in and sampling iterations
%                      and others (see demos)
%
% Outputs: model: 
%         -- model: as above. The outputed model is updated to contain the
%                   parameters values of the final MCMC iteration
%                   parameters as well as the priors
%         -- PropDist: as above. PropDist can be updated (compared to the input one) 
%                     due to the update of the kernel parameters that
%                     influence the proposal 
%         -- samples: the structure that contrains the samples 
%         -- accRates: acceptance rates 
%


BurnInIters = trainOps.Burnin; 
Iters = trainOps.T; 
StoreEvery = trainOps.StoreEvery;
[n D] = size(model.X);
num_stored = floor(Iters/StoreEvery);
samples.F = zeros(num_stored, n);
samples.Fu = zeros(num_stored, n);
samples.LogL = zeros(1, num_stored);
X = model.X;
Y = model.y;
F = model.F; 
Fu = model.Fu; % function values  
Xu = model.Xu; % control locations  
M = size(Fu,2);
n = size(X,1);
U = n+1:n+M; 

% compute the initial values of the likelihood p(Y | F)
loglikHandle = str2func(['logL' model.Likelihood.type]);
oldLogLik = loglikHandle(model.Likelihood, Y, F);
oldLogLik = sum(oldLogLik(:)); 
%if noise == 1
%vars = exp(model.Likelihood.logtheta);
%oldLogLik1 = -0.5*n*log(2*pi*vars) - (0.5/vars)*sum((Y-F(:)).^2);
%elseif noise == 2
%   [pp oldLogLik] = cumGauss(Y,F(:));
%   oldLogLik = sum(oldLogLik(:));
%end


% evaluation of the log prior for the kernel hyperparameters
if strcmp(model.constraints.kernHyper, 'free') 
   lnpriorK = ['ln', model.prior.kernParams.type,'pdf'];
   oldLogPriorK = feval(lnpriorK, model.GP.logtheta, model.prior.kernParams.a, model.prior.kernParams.b);
end

% evaluation of the log prior for the likelihood hyperparameters
if strcmp(model.constraints.likHyper, 'free') & (model.Likelihood.nParams > 0)
   lnpriorLik = ['ln', model.prior.likParams.type,'pdf'];
   oldLogPriorLik = feval(lnpriorLik, model.Likelihood.logtheta, model.prior.likParams.a, model.prior.likParams.b);       
end     


cnt = 0;
acceptF = zeros(1,M);
acceptK = 0;
acceptL = 0;

for it = 1:(BurnInIters + Iters) 
    %
     % sample the control points one-at-a-time 
    Fold = F;
    Fuold =Fu;
    for i=1:M
    %
       % sample the i control point given the rest (Gibbs-like)       
       Fui = randn.*sqrt(PropDist.qF.ku(i)) + PropDist.qF.KInvK(i,:)*Fu([1:i-1, i+1:end])';    
    
       Funew = Fu;
       Funew(i) = Fui;
    
       % resample the whole function
       cmu = PropDist.qF.cmuMinus + Funew*PropDist.qF.KInvKu;
       Fnew = gaussianFastSample(1, cmu, PropDist.qF.L);
    
       % perform an evaluation of the likelihood p(Y | F) 
       newLogLik = loglikHandle(model.Likelihood, Y, Fnew(:));
       newLogLik = sum(newLogLik(:)); 
       %
       
       % Metropolis-Hastings to accept-reject the proposal
       newP = newLogLik;
       oldP = oldLogLik;
       [accept, uprob] = metropolisHastings(newP, oldP, 0, 0);
    
       % visualization
       if trainOps.disp & (mod(it,100) == 0) 
       %[newLogLik oldLogLik]
          if (D==1) 
             trform = 'lin';
             if strcmp(model.Likelihood.type,'Poisson') 
             trform = 'exp';
             end
             plot(X, Y, '+k', 'lineWidth', 2);
             hold on;
             if strcmp(model.Likelihood.type,'Sigmoid') | strcmp(model.Likelihood.type,'Probit')
                 plot(X, zeros(size(X,1)), 'k:')
             end
             plot(X, feval(trform, F), 'g', 'lineWidth', 4);
             %plot(X, F,'or','MarkerSize', 14,'lineWidth', 3);
             pause(0.2);
             plot(Xu, feval(trform, Fu),'or','MarkerSize', 14,'lineWidth', 3);
             plot(Xu(i), feval(trform, Fu(i)),'oy','MarkerSize', 14,'lineWidth', 3);
             %legend(hh,'Current state','Control points');
             set(gca,'FontSize',16);
             plot(Xu(i), feval(trform, Funew(i)), 'md','MarkerSize', 14, 'lineWidth',3);
             plot(X, feval(trform, Fnew), '--b', 'lineWidth', 4); 
             pause(0.5)
             %plot(Xu, Funew, 'md','MarkerSize', 14, 'lineWidth', 3);
             %legend(hh,'Current state','Control points','Proposed state');
             set(gca,'FontSize',16);
             clear hh;
             hold off;
          end
       %
       end
    
       if (it > BurnInIters) 
          acceptF(i) = acceptF(i) + accept;
       end
    
       if accept == 1
          F = Fnew;
          Fu = Funew;
          oldLogLik = newLogLik;
       end
       % 
    end % end control point loop 
    
    
    % sample gp likelihood parameters 
    if strcmp(model.constraints.likHyper, 'free') & (model.Likelihood.nParams > 0)
       % 
       newlogLik = randn.*sqrt(PropDist.lik) + model.Likelihood.logtheta;
       
       Lik1 = model.Likelihood;
       Lik1.logtheta = newlogLik;
       % perform an evaluation of the likelihood p(Y | F) 
       newLogLik = loglikHandle(Lik1, Y, F(:));
       newLogLik = sum(newLogLik(:)); 
       
       newLogPriorLik = feval(lnpriorLik, newlogLik, model.prior.likParams.a, model.prior.likParams.b);
       
       % Metropolis-Hastings to accept-reject the proposal
       oldlogP = oldLogLik + sum(oldLogPriorLik(:));
       newlogP = newLogLik + sum(newLogPriorLik(:)); 
       %
       [accept, uprob] = metropolisHastings(newlogP, oldlogP, 0, 0);
       if accept == 1
        %
          model.Likelihood.logtheta = newlogLik;
          oldLogPriorLik = newLogPriorLik;
          oldLogLik = newLogLik;
        %
       end
        %
       if (it > BurnInIters) 
           acceptL = acceptL + accept;
       end
       %
    end
    
    
    % sample the hyperparameters 
    if strcmp(model.constraints.kernHyper, 'free')
       % 
       newlogK = randn.*sqrt(PropDist.kern) + model.GP.logtheta;
       GPtmp = model.GP; 
       GPtmp.logtheta = newlogK;
       newK = kernCompute(GPtmp, [X; Xu]);    
       [newL,er]=jitterChol(newK);
       newL = newL';
       % evaluate the new log GP prior value 
       invnewL = newL\eye(n+M);
       newLogDetK = 2*sum(log(diag(newL)));
      
       newlogGP = - 0.5*newLogDetK;
       oldlogGP = - 0.5*PropDist.qF.LogDetK;
       temp = invnewL*[F, Fu]'; 
       newlogGP = newlogGP - 0.5*temp'*temp;
       temp = PropDist.qF.invL*[F, Fu]'; 
       oldlogGP = oldlogGP - 0.5*temp'*temp;
       newLogPriorK = feval(lnpriorK, newlogK, model.prior.kernParams.a, model.prior.kernParams.b);
       
       % Metropolis-Hastings to accept-reject the proposal
       oldlogGP = oldlogGP + sum(oldLogPriorK(:));
       newlogGP = newlogGP + sum(newLogPriorK(:)); 
       %
       [accept, uprob] = metropolisHastings(newlogGP, oldlogGP, 0, 0);
       if accept == 1
        %
          model.GP.logtheta = newlogK;
          oldLogPriorK = newLogPriorK;
          % update proposal for F
          PropDist.qF.K = newK;
          PropDist.qF.invL = invnewL; 
          PropDist.qF.LogDetK = newLogDetK;
          [cmuMinus, cSigma, KInvKu] = gaussianFastConditional(PropDist.qF.m', newK, 1:n, U);
          [L,er]=jitterChol(cSigma);
          if er>0, L = real(sqrtm(cSigma)); end
          PropDist.qF.cmuMinus = cmuMinus; 
          PropDist.qF.cSigma = cSigma;
          PropDist.qF.KInvKu = KInvKu;
          PropDist.qF.L = L;
          for i=1:M
          %  
             G = [1:i-1, i+1:M];  
             [alpha(i), ku(i), KInvK(i,:)] = gaussianFastConditional(PropDist.qF.m(U)', PropDist.qF.K(U,U), i, G);
          %
          end
          PropDist.qF.alpha = alpha;
          PropDist.qF.ku = ku;
          PropDist.qF.KInvK = KInvK;  
          %
        end
        % 
        %[accept acceptH oldlogGP]
        
        if (it > BurnInIters) 
           acceptK = acceptK + accept;
        end
        %
    end
  
    % keep samples after burn in 
    if (it > BurnInIters)  & (mod(it,StoreEvery) == 0)
    %
        cnt = cnt + 1;
        samples.F(cnt,:) = F;
        samples.Fu(cnt,:) = F;
        if strcmp(model.constraints.kernHyper, 'free')
           samples.kernLogtheta(cnt,:) = model.GP.logtheta;    
        end 
        if strcmp(model.constraints.likHyper, 'free') & (model.Likelihood.nParams >0)
           samples.likLogtheta(cnt,:) = model.Likelihood.logtheta;  
        end
        samples.LogL(cnt) = oldLogLik;
    %
    end
    %        
end

%
model.F = F;
model.Fu = Fu;
accRates.F = (acceptF/Iters)*100;
if strcmp(model.constraints.kernHyper, 'free')
   accRates.kern = (acceptK/Iters)*100;
else
   accRates.kern = 100;
end

if strcmp(model.constraints.likHyper, 'free') & (model.Likelihood.nParams >0)
   accRates.lik = (acceptL/Iters)*100;
else
   accRates.lik = 100;
end
                                                                                                                         gpsampCreate.m                                                                                      0000700 0003466 0000024 00000003523 11275534210 013355  0                                                                                                    ustar   mtitsias                        games                                                                                                                                                                                                                  function model = gpsampCreate(y, X, options)
%
%

model.type = 'gpmodel';
model.Likelihood.type = options.Likelihood;

[numData D] = size(X);

model.y = y; 
model.X = X;
model.numData = numData;
model.D = D;

switch model.Likelihood.type
    case 'Gaussian' % standard regression
         %
         model.Likelihood.nParams = 1; % parameters (excluding gp function F)
         model.Likelihood.logtheta = log(0.05); % log(sigma2) 
         %
    case 'Probit'  % binary classifcation
         %
         model.Likelihood.nParams = 0;
         %
    case 'Sigmoid' % binary classification
         %
         model.Likelihood.nParams = 0;
         %
    case 'Poisson' % for counts data      
         %
         model.Likelihood.nParams = 0;  
         %
    case 'ODE'
         % %%
end     

model.constraints.kernHyper = options.constraints.kernHyper;
model.constraints.likHyper = options.constraints.likHyper;

model.GP.type = 'rbf';
% kernel hyperparameters 
model.GP.logtheta = [2*log((max(X) - min(X))*0.2) 0];
model.GP.nParams = D+1;

% prior over the likelihood parameters  
% (all are assumed to take non-negative values, so they are represented
% in the log space and prior is define there)
model.prior.likParams.type = 'normal';
model.prior.likParams.constraint = 'positive';
model.prior.likParams.priorSpace = 'log';
model.prior.likParams.a = 0; % mean 
model.prior.likParams.b = 2; % variance
 

% prior over GP kernel hyperparameters 
% (all are assumed to take non-negative values, so they are represented
% in the log space and prior is define there)
model.prior.kernParams.type = 'normal';
model.prior.kernParams.constraint = 'positive';
model.prior.kernParams.priorSpace = 'log';
model.prior.kernParams.a = 0; % mean 
model.prior.kernParams.b = 2; % variance
 
% GP latent function values needed to define the likelihood
model.F = zeros(1,model.numData);
                                                                                                                                                                             gpsampOptions.m                                                                                     0000700 0003466 0000024 00000001112 11275534171 013603  0                                                                                                    ustar   mtitsias                        games                                                                                                                                                                                                                  function options = gpsampOptions(gplik) 
%
%

options.kern = 'rbf';
options.constraints.kernHyper = 'free'; % 'free' or 'fixed'

switch gplik 
    case 'regression'
        % Gaussian likelihood 
        options.Likelihood = 'Gaussian';
    case 'classification' % binary classification 
        % Probit likelihood
        options.Likelihood = 'Sigmoid';
    case 'poissonRegr' % binary classification 
        % Probit likelihood
        options.Likelihood = 'Poisson';
    case 'ODE' 
        % not included 
end

options.constraints.likHyper = 'free'; % 'free' or 'fixed'

        
                                                                                                                                                                                                                                                                                                                                                                                                                                                      gpsampPlot.m                                                                                        0000600 0003466 0000024 00000002605 11275534204 013072  0                                                                                                    ustar   mtitsias                        games                                                                                                                                                                                                                  function gpsampPlot(model, samples)
%
% 



% plot likelihood hyperparameters 
if isfield(samples,'likLogtheta')
   %
   for i=1:model.Likelihood.nParams; 
       figure;
       hist(exp(samples.likLogtheta(:,i)),50); 
       titlestring = 'likelihood hyperp: ';
       titlestring = [titlestring, num2str(i)]; 
       title(titlestring,'fontsize', 20);
   end
   %
end


% plot kernel hyperparameters 
if isfield(samples, 'kernLogtheta')
   %
   for i=1:model.GP.nParams; 
       figure;
       hist(exp(samples.kernLogtheta(:,i)),50); 
       titlestring = 'kernel hyperp: ';
       titlestring = [titlestring, num2str(i)]; 
       title(titlestring,'fontsize', 20);
   end
   %
end


% plot GP function 
if model.D == 1
   %
   %[sortX inX] = sort(model.X);
  
   if strcmp(model.Likelihood.type,'Poisson')
      mu = mean(exp(samples.F))';
      stds = sqrt(var(exp(samples.F)))';
   else
      mu = mean(samples.F)';
      stds = sqrt(var(samples.F))';
   end
   
   %stdsWithnoise = sqrt(var(samples.F))';
    
   figure
   hold on;
   fillColor = [0.7 0.7 0.7];
   %fillColor = [0.8 0.8 0.8];  % for the paper
   fill([model.X; model.X(end:-1:1)], [mu; mu(end:-1:1)]...
            + 2*[stds; -stds(end:-1:1)], fillColor,'EdgeColor',fillColor);
   plot(model.X, mu,'b','lineWidth',3);
     
   %if strcmp(model.Likelihood.type,'Gaussian')
   plot(model.X, model.y, '+k', 'lineWidth', 1);
   %end
   %
end

                                                                                                                           stats.m                                                                                             0000600 0003466 0000024 00000000312 11312724454 012073  0                                                                                                    ustar   mtitsias                        games                                                                                                                                                                                                                  


for n=1:size(testGene,2)
    %
    for j=1:5
        ok = testGene{n}.Weights(j,:);
        ok(ok>=-0.0316 & ok<=0.0316) = 0; 
        ok(ok~=0)=1;
        prob(n,j) = sum(ok)/3000;
    end
    %
end                                                                                                                                                                                                                                                                                                                      toolbox/computeKLs.m                                                                                0000600 0003466 0000024 00000001045 11275533153 014517  0                                                                                                    ustar   mtitsias                        games                                                                                                                                                                                                                  function KLdiv = computeKLs(model, samples)
% auxiliary function
%
%
%

X = model.X;
Y = model.y;
  
[N D] = size(X);
  
sigma2n = exp(model.Likelihood.logtheta(1));
sigmaf = exp(model.GP.logtheta(D+1));

% GP posterior 
Knn = kernCompute(model.GP, X, X);

jitter = 0.000001*mean(diag(Knn));

L = chol(Knn + sigma2n*eye(N))';
alpha = L'\(L\Y);
muGP = Knn*alpha;   
v = L\Knn;  
CovarGP = Knn - v'*v + sigma2n*eye(N); 
     
Fs = samples.F;
muMCMC = mean(Fs);
CovarMCMC = cov(Fs,1)+sigma2n*eye(N);
    
KLdiv = kl(muGP, CovarGP, muMCMC, CovarMCMC);

                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           toolbox/gaussianConditional.m                                                                       0000700 0003466 0000024 00000001255 11275533155 016435  0                                                                                                    ustar   mtitsias                        games                                                                                                                                                                                                                  function [cmu, cSigma] = gaussianConditional(X, condInds, mu, Sigma)
% Description: It finds the conditional distribution parameters
%              (mean and covariance) of  P(X_rest| X_condInds) 
%
%
% OUTPUT: The conditional mean: cmu, the conditional covariance: cSigma  

d = size(Sigma,1);
d1 = size(condInds(:)',2);
d2 = round(d-d1);

nonCondInds = setdiff([1:d], condInds);
X1 = X(condInds); 
X2 = X(nonCondInds);
mu1 = mu(condInds); 
mu2 = mu(nonCondInds); 

Sigma = Sigma([nonCondInds, condInds],:);
Sigma = Sigma(:,[nonCondInds, condInds]);

B = Sigma(d2+1:end,d2+1:end);
C = Sigma(1:d2, (d2+1):end);

ok = B\C';
cmu = mu2 + (X1 - mu1)*ok;
cSigma = Sigma(1:d2,1:d2) - C*ok;

                                                                                                                                                                                                                                                                                                                                                   toolbox/gaussianConditionalDisplay.m                                                                0000700 0003466 0000024 00000001115 11275533155 017756  0                                                                                                    ustar   mtitsias                        games                                                                                                                                                                                                                  function e = gaussianConditionalDisplay(Nc, mu, K)
%
% Illustration of sampling from Gaussian Process conditional 
%

e = 1;

d = size(K,1);


Kb = 0*eye(d);
F = sampleGaussian(1,mu,K+Kb);


M = floor(d/Nc);
obs = [1:M:d]; 
grid = setdiff([1:d], obs);
obs

% compute the conditional Gaussian
[cmu, cSigma] = gaussianConditional(F, obs, mu, K+Kb);

plot(F);
hold on;
plot(obs,F(obs),'ko','MarkerSize',10);

FF = F;
while 1
    FF(grid) = sampleGaussian(1,cmu,cSigma);
    plot(FF,'r')
    plot(F,'b','LineWidth',2);
    plot(obs,F(obs),'go','MarkerSize',10, 'LineWidth', 2);
    pause;
end
                                                                                                                                                                                                                                                                                                                                                                                                                                                   toolbox/gaussianFastConditional.m                                                                   0000700 0003466 0000024 00000002366 11275533156 017260  0                                                                                                    ustar   mtitsias                        games                                                                                                                                                                                                                  function [cmu, cSigma, KInvK] =  gaussianFastConditional(mu, K, Express, Given)
% Partial computation of the conditional X(Express) | X(Given)
% Inputs: 
%      -- mu the mean of the Gaussian as a row vector 
%      -- K the covariance matrix of the Gaussian  
%      -- Express and Given are disjoint sets of indices pointing 
%         in the random vector 
%  
% The function writes the Gaussian N(mu,K) in the form 
%
%   [X(Express)]      mu(Express)  [ A  C  ...]
%   [X(Given)  ]  = N(mu(Given),   [ C^T B ...])
%   [X(Rest)   ]      mu(Rest)     [ ...   ...]
% 
% and then computes elements of the conditional:
%        N(mu(Express) + C*B^(-1)*[X(Given) - mu(Given)], A - C*B^(-1)*C') 
%
% Outputs:
%       -- cmu    =  mu(Express)' - mu(Given)'*(B^(-1)*C') 
%       -- cSigma =  A - C*B^(-1)*C'
%       -- KInvK  =  B^(-1)*C'
%

jitter = 10e-8;

SizF = size(mu,2);
SizE = size(Express,2);
SizG = size(Given,2);
Ignore = setdiff([1:SizF], [Express, Given]);


muE = mu(Express); 
muG = mu(Given);

K = K([Express, Given, Ignore],:);
K = K(:,[Express, Given, Ignore]);
K = K(1:SizE+SizG,1:SizE+SizG);

B = K(SizE+1:end,SizE+1:end) + jitter*K(1,1)*eye(SizG);
C = K(1:SizE, (SizE+1):end);

KInvK = B\C';
cmu = muE - muG*KInvK;
cSigma = K(1:SizE,1:SizE) - C*KInvK;
                                                                                                                                                                                                                                                                          toolbox/gaussianFastSample.m                                                                        0000700 0003466 0000024 00000001071 11305234710 016212  0                                                                                                    ustar   mtitsias                        games                                                                                                                                                                                                                  function X = gaussianFastSample(N, mu, L)
% Function: Generates N samples from a d-dimension
%           Gaussian distribution.
%
% Inputs:  
%       -- N  mumber of samples to generate
%       -- mu 1 x d vector that is the mean of the Gaussian
%       -- L  Choleski decomposition of the covariance matrix (or the square root of 
%          this matrix in case the original covariance was positive semidefinite)
% Outputs: 
%        -- X (N x d) the outputs random vectors  

d = size(L,1);
X = randn(N,d); 
X = real(X*L);
X = X+mu(ones(N, 1), :); %repmat(mu,[N,1]);
                                                                                                                                                                                                                                                                                                                                                                                                                                                                       toolbox/gaussianSample.m                                                                            0000700 0003466 0000024 00000001022 11275533157 015405  0                                                                                                    ustar   mtitsias                        games                                                                                                                                                                                                                  function X = gaussianSample(N, mu, Sigma)
%Description:  It generates N samples from a d-dimensional
%              Gaussian distribution.
%
%N  : Number of samples to generate
%mu : mean
%Sigma : covariance matrix for the samples, should be positive definite

d = size(Sigma,1);

X=randn(N,d); 

% Cholesky decomposition for a positive definite Covariance
% If the covariance is not positive definite uses the square root 
[L,r]=chol(Sigma);
if r>0 
    X = real(X*sqrtm(Sigma));
else
    X = X*L;
end
X=X+repmat(mu,[N,1]);

                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  toolbox/jitterChol.m                                                                                0000700 0003466 0000024 00000000552 11275533160 014541  0                                                                                                    ustar   mtitsias                        games                                                                                                                                                                                                                  function [L er] = jitterChol(K)
% function [L er] = jitterChol(K)
%
% Description:  Computing Choleski decomposition by adding jitter  
%              when the matrix is semipositive definite  
%

jitter = 1e-8;
m = size(K,1); 
[L er] = chol(K);
if er > 0 % add jitter
   warning('Jitter added'); 
   K = K + jitter*mean(diag(K))*eye(m);
   [L er] = chol(K);
end                                                                                                                                                      toolbox/kernCompute.m                                                                               0000700 0003466 0000024 00000001367 11306047120 014722  0                                                                                                    ustar   mtitsias                        games                                                                                                                                                                                                                  function Knm = covfuncCompute(GPprior, X, Xu)
%function Knm = covfuncCompute(logtheta, X, Xu)
%
%Description:  It computes the covariance function between 
%              two set of inputs points: X and Xu.  
%

%
%Supported covariance functions:  RBF and ARD kernel.   

jitter = exp(2*GPprior.logtheta(end));

[n D] = size(X);
logtheta = GPprior.logtheta(:);
sigmaf = exp(logtheta(D+1));
X = X ./ repmat(exp(logtheta(1:D))',n,1);

if nargin == 3
   [m,D] = size(Xu);   
   Xu = Xu ./ repmat(exp(logtheta(1:D))',m,1);
   %
   Knm = -2*Xu*X' + repmat(sum(X.*X,2)',m,1) + repmat(sum(Xu.*Xu,2),1,n);
   Knm = sigmaf*exp(-0.5*Knm');
else
   Knm = -2*X*X' + repmat(sum(X.*X,2)',n,1) + repmat(sum(X.*X,2),1,n);
   Knm = sigmaf*exp(-0.5*Knm') + jitter*eye(n); 
end
                                                                                                                                                                                                                                                                         toolbox/kl.m                                                                                        0000700 0003466 0000024 00000000730 11275533161 013037  0                                                                                                    ustar   mtitsias                        games                                                                                                                                                                                                                  function out = KL(mu0, Sigma0, mu1, Sigma1)
% function out = KL(mu0, Sigma0, mu1, Sigma1)
%
% Description :  Computes the  KL divergence between two Gaussians distribution 
%                
%

mu0 = mu0(:);
mu1 = mu1(:);

N = size(Sigma0,1);

L0 = jitterChol(Sigma0)'; 
L1 = jitterChol(Sigma1)'; 
invSigma1 = L1'\(L1\eye(N));

out = sum(log(diag(L1))) - sum(log(diag(L0)))...
      + 0.5*trace(invSigma1*Sigma0) + 0.5*((mu1-mu0)'*(invSigma1*(mu1-mu0)))...
      - 0.5*N;
                                        toolbox/lin.m                                                                                       0000700 0003466 0000024 00000000040 11275533161 013205  0                                                                                                    ustar   mtitsias                        games                                                                                                                                                                                                                  function fx = lin(x)
%

fx = x;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                toolbox/lnLaplacepdf.m                                                                              0000700 0003466 0000024 00000000205 11275533162 015014  0                                                                                                    ustar   mtitsias                        games                                                                                                                                                                                                                  function logpdf = lnLaplacepdf(x,mu,beta)
% natural logarithm of the normal density
%
%

logpdf = -log(2*beta) - (1/beta).*abs(x-mu);                                                                                                                                                                                                                                                                                                                                                                                           toolbox/lngammapdf.m                                                                                0000700 0003466 0000024 00000000303 11312467750 014535  0                                                                                                    ustar   mtitsias                        games                                                                                                                                                                                                                  function logpdf = lngammapdf(x,prior)
%
%

%z = x ./ b;
%u = (a - 1) .* log(z) - z - gammaln(a);
%y = exp(u) ./ b;

a = prior.a;
b = prior.b;

logpdf = (a-1)*log(x) - b*x + a*log(b) - gammaln(a);                                                                                                                                                                                                                                                                                                                             toolbox/lnnormalpdf.m                                                                               0000700 0003466 0000024 00000000244 11312467740 014746  0                                                                                                    ustar   mtitsias                        games                                                                                                                                                                                                                  function logpdf = lnnormalpdf(x, prior)
% natural logarithm of the normal density
%
%

logpdf = -0.5*log(2*pi*prior.sigma2) - (0.5/prior.sigma2).*((x-prior.mu).^2);                                                                                                                                                                                                                                                                                                                                                            toolbox/lnspikeNormalpdf.m                                                                          0000600 0003466 0000024 00000000663 11320454422 015736  0                                                                                                    ustar   mtitsias                        games                                                                                                                                                                                                                  function logpdf = lnspikeNormalpdf(x, prior)
% natural logarithm of the normal density
%
%

%logpdf = -0.5*log(2*pi*prior.sigma2) - (0.5/prior.sigma2).*((x-prior.mu).^2);
%ok = exp(logpdf);
g1 = (1/sqrt(2*pi*prior.sigma2))*exp(- (0.5/prior.sigma2).*((x-prior.mu).^2));
g2 = (1/sqrt(2*pi*prior.spikeSigma2))*exp(- (0.5/prior.spikeSigma2).*((x-prior.spikeMu).^2));

%[ok; g1; g2]
%pause

logpdf = log((1-prior.pis).*g1 + prior.pis.*g2);                                                                              toolbox/lnspikeTruncNormalpdf.m                                                                     0000600 0003466 0000024 00000001062 11320643470 016747  0                                                                                                    ustar   mtitsias                        games                                                                                                                                                                                                                  function logpdf = lnspikeTruncNormalpdf(x, prior)
%  A mixture model with half truncated Gaussian (x>=mu) and a spike Gaussian 
%

% !!!! This will only work if prior.mu = prior.spikeMu = 0 !!!!


if (min(x) >= prior.mu) &  (min(x) >= prior.spikeMu) 
   g1 = (1/sqrt(2*pi*prior.sigma2))*exp(- (0.5/prior.sigma2).*((x-prior.mu).^2));
   g1 = 2*g1;

   g2 = (1/sqrt(2*pi*prior.spikeSigma2))*exp(- (0.5/prior.spikeSigma2).*((x-prior.spikeMu).^2));

   g2 = 2*g2;
   logpdf = log((1-prior.pis).*g1 + prior.pis.*g2);
else
    disp('x cannot be smaller than zero')
end                                                                                                                                                                                                                                                                                                                                                                                                                                                                              toolbox/lntruncNormalpdf.m                                                                          0000700 0003466 0000024 00000000444 11320462654 015762  0                                                                                                    ustar   mtitsias                        games                                                                                                                                                                                                                  function logpdf = lntruncNormalpdf(x, prior)
% natural logarithm of half truncated  normal density (x > mu)
%
%

sigma  = sqrt(prior.sigma2);

z = (x-prior.mu)./sigma;

logpdf = -0.5*log(2*pi*prior.sigma2) - 0.5.*(z.^2); 

den = 1 - normcdf(0, prior.mu, sigma); 

logpdf = logpdf - log(den);                                                                                                                                                                                                                             toolbox/lntwoMixNormalpdf.m                                                                         0000600 0003466 0000024 00000000664 11320474371 016120  0                                                                                                    ustar   mtitsias                        games                                                                                                                                                                                                                  function logpdf = lntwoMixNormalpdf(x, prior)
% natural logarithm of the normal density
%
%

%logpdf = -0.5*log(2*pi*prior.sigma2) - (0.5/prior.sigma2).*((x-prior.mu).^2);
%ok = exp(logpdf);
g1 = (1/sqrt(2*pi*prior.sigma2))*exp(- (0.5/prior.sigma2).*((x-prior.mu).^2));
g2 = (1/sqrt(2*pi*prior.spikeSigma2))*exp(- (0.5/prior.spikeSigma2).*((x-prior.spikeMu).^2));

%[ok; g1; g2]
%pause

logpdf = log((1-prior.pis).*g1 + prior.pis.*g2);                                                                             toolbox/logLGaussian.m                                                                              0000700 0003466 0000024 00000000443 11275533164 015025  0                                                                                                    ustar   mtitsias                        games                                                                                                                                                                                                                  function out = logLGaussian(lik, Y, F)
%function out = logLGaussian(Y,F,sigma2)
%
%Description: Log likelihood for the one dimensional 
%             Gaussian distribution   
%

sigma2 = exp(lik.logtheta); 

n = size(Y(:),1);
out = -0.5*n*log(2*pi*sigma2) - (0.5/sigma2)*sum((Y(:)-F(:)).^2);                                                                                                                                                                                                                             toolbox/logLPoisson.m                                                                               0000700 0003466 0000024 00000000301 11275533164 014676  0                                                                                                    ustar   mtitsias                        games                                                                                                                                                                                                                  function out = logLPoisson(lik, Y, F)
%function out = logLPoisson(lik, Y,F)
%
%Description: Poisson GP log-likelihood useful for counts data 
%

out = Y(:).*F(:) - exp(F(:)) - gammaln(Y(:)+1);                                                                                                                                                                                                                                                                                                                                toolbox/logLProbit.m                                                                                0000700 0003466 0000024 00000000171 11275533164 014510  0                                                                                                    ustar   mtitsias                        games                                                                                                                                                                                                                  function out = logLProbit(lik, Y, F)
%
%

yf = Y(:).*F(:); 
     
out = (1+erf(yf/sqrt(2)))/2;         
out = log(out);

                                                                                                                                                                                                                                                                                                                                                                                                       toolbox/logLSigmoid.m                                                                               0000700 0003466 0000024 00000000306 11275533164 014644  0                                                                                                    ustar   mtitsias                        games                                                                                                                                                                                                                  function out = logLSigmoid(lik, Y, F)
%function out = logLPoisson(lik, Y, F)
%
%Description: Sigmoid GP log-likelihood useful for binary classification
%

YF= Y(:).*F(:); 
out = -log(1 + exp(-YF));
                                                                                                                                                                                                                                                                                                                          toolbox/mcmcOptions.m                                                                               0000700 0003466 0000024 00000002440 11275533165 014730  0                                                                                                    ustar   mtitsias                        games                                                                                                                                                                                                                  function mcmcoptions = mcmcOptions(algorithm) 
%
%

% MCMC OPTIONS (you can change these options if you like)   
mcmcoptions.train.StoreEvery = 100; % keep samples every StoreEvery MCMC iterations
mcmcoptions.train.Burnin = 5000;  % burn in time
mcmcoptions.train.T = 50000; % sampling time
mcmcoptions.train.Store = 0;  % store the results regularly in a file         
mcmcoptions.train.disp = 0;  % display information during sampling (if applicable)

% define also options required when perform the adaptive MCMC phase 
switch algorithm 
    case 'controlPnts'
        %        
        mcmcoptions.adapt.T = 200;          
        mcmcoptions.adapt.Burnin = 100;
        mcmcoptions.adapt.StoreEvery = 10; 
        mcmcoptions.adapt.disp = 1;
        mcmcoptions.adapt.initialNumContrPnts = 3; 
        mcmcoptions.adapt.incrNumContrBy = 1;
        %
    case 'localRegion'
        %
        mcmcoptions.adapt.T = 200;          
        mcmcoptions.adapt.Burnin = 100;
        mcmcoptions.adapt.StoreEvery = 10; 
        mcmcoptions.adapt.disp = 1;
        mcmcoptions.adapt.initialNumRegions = 3; 
        mcmcoptions.adapt.incrNumRegionBy = 1;
        %
    case 'gibbs' 
        % for Gibbs sampling or Gibbs-like algorthm no adaptive MCMC is
        % needed
        mcmcoptions.adapt = [];
        %
end                                                                                                                                                                                                                                toolbox/metropolisHastings.m                                                                        0000700 0003466 0000024 00000000433 11275533166 016334  0                                                                                                    ustar   mtitsias                        games                                                                                                                                                                                                                  function [accept, A] = metropolisHastings(newLogPx, oldLogPx, newLogProp, oldLogProp)
%
%Desctiption:  The general Metropolis-Hastings step 
%

A = newLogPx + oldLogProp - oldLogPx - newLogProp;
A = exp(A);

accept = 0;
u = rand;
if u < A
   accept = 1;
end

                                                                                                                                                                                                                                               toolbox/minimize.m                                                                                  0000700 0003466 0000024 00000021445 11275533166 014265  0                                                                                                    ustar   mtitsias                        games                                                                                                                                                                                                                  function [X, fX, i] = minimize(X, f, length, varargin)

% Minimize a differentiable multivariate function. 
%
% Usage: [X, fX, i] = minimize(X, f, length, P1, P2, P3, ... )
%
% where the starting point is given by "X" (D by 1), and the function named in
% the string "f", must return a function value and a vector of partial
% derivatives of f wrt X, the "length" gives the length of the run: if it is
% positive, it gives the maximum number of line searches, if negative its
% absolute gives the maximum allowed number of function evaluations. You can
% (optionally) give "length" a second component, which will indicate the
% reduction in function value to be expected in the first line-search (defaults
% to 1.0). The parameters P1, P2, P3, ... are passed on to the function f.
%
% The function returns when either its length is up, or if no further progress
% can be made (ie, we are at a (local) minimum, or so close that due to
% numerical problems, we cannot get any closer). NOTE: If the function
% terminates within a few iterations, it could be an indication that the
% function values and derivatives are not consistent (ie, there may be a bug in
% the implementation of your "f" function). The function returns the found
% solution "X", a vector of function values "fX" indicating the progress made
% and "i" the number of iterations (line searches or function evaluations,
% depending on the sign of "length") used.
%
% The Polack-Ribiere flavour of conjugate gradients is used to compute search
% directions, and a line search using quadratic and cubic polynomial
% approximations and the Wolfe-Powell stopping criteria is used together with
% the slope ratio method for guessing initial step sizes. Additionally a bunch
% of checks are made to make sure that exploration is taking place and that
% extrapolation will not be unboundedly large.
%
% See also: checkgrad 
%
% Copyright (C) 2001 - 2006 by Carl Edward Rasmussen (2006-09-08).

INT = 0.1;    % don't reevaluate within 0.1 of the limit of the current bracket
EXT = 3.0;                  % extrapolate maximum 3 times the current step-size
MAX = 20;                         % max 20 function evaluations per line search
RATIO = 10;                                       % maximum allowed slope ratio
SIG = 0.1; RHO = SIG/2; % SIG and RHO are the constants controlling the Wolfe-
% Powell conditions. SIG is the maximum allowed absolute ratio between
% previous and new slopes (derivatives in the search direction), thus setting
% SIG to low (positive) values forces higher precision in the line-searches.
% RHO is the minimum allowed fraction of the expected (from the slope at the
% initial point in the linesearch). Constants must satisfy 0 < RHO < SIG < 1.
% Tuning of SIG (depending on the nature of the function to be optimized) may
% speed up the minimization; it is probably not worth playing much with RHO.

% The code falls naturally into 3 parts, after the initial line search is
% started in the direction of steepest descent. 1) we first enter a while loop
% which uses point 1 (p1) and (p2) to compute an extrapolation (p3), until we
% have extrapolated far enough (Wolfe-Powell conditions). 2) if necessary, we
% enter the second loop which takes p2, p3 and p4 chooses the subinterval
% containing a (local) minimum, and interpolates it, unil an acceptable point
% is found (Wolfe-Powell conditions). Note, that points are always maintained
% in order p0 <= p1 <= p2 < p3 < p4. 3) compute a new search direction using
% conjugate gradients (Polack-Ribiere flavour), or revert to steepest if there
% was a problem in the previous line-search. Return the best value so far, if
% two consecutive line-searches fail, or whenever we run out of function
% evaluations or line-searches. During extrapolation, the "f" function may fail
% either with an error or returning Nan or Inf, and minimize should handle this
% gracefully.

if max(size(length)) == 2, red=length(2); length=length(1); else red=1; end
if length>0, S='Linesearch'; else S='Function evaluation'; end 

i = 0;                                            % zero the run length counter
ls_failed = 0;                             % no previous line search has failed
[f0 df0] = feval(f, X, varargin{:});          % get function value and gradient
fX = f0;
i = i + (length<0);                                            % count epochs?!
s = -df0; d0 = -s'*s;           % initial search direction (steepest) and slope
x3 = red/(1-d0);                                  % initial step is red/(|s|+1)

while i < abs(length)                                      % while not finished
  i = i + (length>0);                                      % count iterations?!

  X0 = X; F0 = f0; dF0 = df0;                   % make a copy of current values
  if length>0, M = MAX; else M = min(MAX, -length-i); end

  while 1                             % keep extrapolating as long as necessary
    x2 = 0; f2 = f0; d2 = d0; f3 = f0; df3 = df0;
    success = 0;
    while ~success && M > 0
      try
        M = M - 1; i = i + (length<0);                         % count epochs?!
        [f3 df3] = feval(f, X+x3*s, varargin{:});
        if isnan(f3) || isinf(f3) || any(isnan(df3)+isinf(df3)), error(''), end
        success = 1;
      catch                                % catch any error which occured in f
        x3 = (x2+x3)/2;                                  % bisect and try again
      end
    end
    if f3 < F0, X0 = X+x3*s; F0 = f3; dF0 = df3; end         % keep best values
    d3 = df3'*s;                                                    % new slope
    if d3 > SIG*d0 || f3 > f0+x3*RHO*d0 || M == 0  % are we done extrapolating?
      break
    end
    x1 = x2; f1 = f2; d1 = d2;                        % move point 2 to point 1
    x2 = x3; f2 = f3; d2 = d3;                        % move point 3 to point 2
    A = 6*(f1-f2)+3*(d2+d1)*(x2-x1);                 % make cubic extrapolation
    B = 3*(f2-f1)-(2*d1+d2)*(x2-x1);
    x3 = x1-d1*(x2-x1)^2/(B+sqrt(B*B-A*d1*(x2-x1))); % num. error possible, ok!
    if ~isreal(x3) || isnan(x3) || isinf(x3) || x3 < 0 % num prob | wrong sign?
      x3 = x2*EXT;                                 % extrapolate maximum amount
    elseif x3 > x2*EXT                  % new point beyond extrapolation limit?
      x3 = x2*EXT;                                 % extrapolate maximum amount
    elseif x3 < x2+INT*(x2-x1)         % new point too close to previous point?
      x3 = x2+INT*(x2-x1);
    end
  end                                                       % end extrapolation

  while (abs(d3) > -SIG*d0 || f3 > f0+x3*RHO*d0) && M > 0  % keep interpolating
    if d3 > 0 || f3 > f0+x3*RHO*d0                         % choose subinterval
      x4 = x3; f4 = f3; d4 = d3;                      % move point 3 to point 4
    else
      x2 = x3; f2 = f3; d2 = d3;                      % move point 3 to point 2
    end
    if f4 > f0           
      x3 = x2-(0.5*d2*(x4-x2)^2)/(f4-f2-d2*(x4-x2));  % quadratic interpolation
    else
      A = 6*(f2-f4)/(x4-x2)+3*(d4+d2);                    % cubic interpolation
      B = 3*(f4-f2)-(2*d2+d4)*(x4-x2);
      x3 = x2+(sqrt(B*B-A*d2*(x4-x2)^2)-B)/A;        % num. error possible, ok!
    end
    if isnan(x3) || isinf(x3)
      x3 = (x2+x4)/2;               % if we had a numerical problem then bisect
    end
    x3 = max(min(x3, x4-INT*(x4-x2)),x2+INT*(x4-x2));  % don't accept too close
    [f3 df3] = feval(f, X+x3*s, varargin{:});
    if f3 < F0, X0 = X+x3*s; F0 = f3; dF0 = df3; end         % keep best values
    M = M - 1; i = i + (length<0);                             % count epochs?!
    d3 = df3'*s;                                                    % new slope
  end                                                       % end interpolation

  if abs(d3) < -SIG*d0 && f3 < f0+x3*RHO*d0          % if line search succeeded
    X = X+x3*s; f0 = f3; fX = [fX' f0]';                     % update variables
    %fprintf('%s %6i;  Value %4.6e\r', S, i, f0);
    s = (df3'*df3-df0'*df3)/(df0'*df0)*s - df3;   % Polack-Ribiere CG direction
    df0 = df3;                                               % swap derivatives
    d3 = d0; d0 = df0'*s;
    if d0 > 0                                      % new slope must be negative
      s = -df0; d0 = -s'*s;                  % otherwise use steepest direction
    end
    x3 = x3 * min(RATIO, d3/(d0-realmin));          % slope ratio but max RATIO
    ls_failed = 0;                              % this line search did not fail
  else
    X = X0; f0 = F0; df0 = dF0;                     % restore best point so far
    if ls_failed || i > abs(length)         % line search failed twice in a row
      break;                             % or we ran out of time, so we give up
    end
    s = -df0; d0 = -s'*s;                                        % try steepest
    x3 = 1/(1-d0);                     
    ls_failed = 1;                                    % this line search failed
  end
end
%fprintf('\n');
                                                                                                                                                                                                                           toolbox/trace_CondCov.m                                                                             0000700 0003466 0000024 00000002735 11306221062 015136  0                                                                                                    ustar   mtitsias                        games                                                                                                                                                                                                                  function [f, df] = trace_CondCov(Xu, X, logtheta, FixedFirst)
% Description:  Function and gradient evalaution with respect to  pseudo-inputs 
%               of the trace of the conditional prior (given inducing variables)
%               for the expenential kernel with varied length-scale (ARD kernel)
%

% number of examples and dimension of input space
[n, D] = size(X);
m = round(size(Xu,1)/D);
jitter = 1e-8;

Xu = reshape(Xu,m,D);

sigmaf = exp(2*logtheta(D+1));
X = X ./ repmat(exp(logtheta(1:D)),n,1);
Xu = Xu ./ repmat(exp(logtheta(1:D)),m,1);

Kmm = Xu*Xu';
Kmm = repmat(diag(Kmm),1,m) + repmat(diag(Kmm)',m,1) - 2*Kmm;
Kmm = sigmaf*exp(-0.5*Kmm);
Knm = -2*Xu*X' + repmat(sum(X.*X,2)',m,1) + repmat(sum(Xu.*Xu,2),1,n);
Knm = sigmaf*exp(-0.5*Knm');


[Lkmm er] = jitterChol(Kmm);
%if er > 0 % add jitter
%   Kmm = Kmm + 1e-7*sigmaf*eye(m);
%   disp('jitter added')
%   Lkmm = chol(Kmm);
%end
%

Cnm1 = Knm/Lkmm;
Cmnmn = Cnm1'*Cnm1;

% value of the objective function
f = n*sigmaf - sum(diag(Cmnmn));

% compute derivatives
Pmnmn = (Lkmm\Cmnmn)';
BB1 = Lkmm\Pmnmn; 
BB1 = Kmm.*BB1;

Cnm1 = (Lkmm\Cnm1')';
Cnm1 = Cnm1.*Knm;
%
for d=1:D
    %
    % pseudo inputs derivatives
    Knm = -((repmat(Xu(:,d)',n,1)-repmat(X(:,d),1,m))/exp(logtheta(d)));
    Kmm = -((repmat(Xu(:,d)',m,1)-repmat(Xu(:,d),1,m))/exp(logtheta(d)));       
    
    dXu(:,d) = -(sum(Knm.*Cnm1,1) - sum(Kmm.*BB1,1))'; 
    %
end

%
dXu = 2*dXu;
if FixedFirst == 1
    dXu(1,:) = zeros(1,D);
end
df = reshape(dXu, m*D, 1);
                                   activFuncts/jointactFunc.m                                                                          0000600 0003466 0000024 00000002263 11312233754 015662  0                                                                                                    ustar   mtitsias                        games                                                                                                                                                                                                                  function fx = jointactFunc(LikParams,ff,J)
%
%
%
%

SizF = size(ff,2);

%NumGenes = size(J(:),1);
W = LikParams.W(J,:);
W0 = LikParams.W0(J); 

%
switch LikParams.jointAct
    case 'sigmoid'
       %
       fx = W*ff + W0(:, ones(1, SizF)); %repmat(W0,[1 SizF]);
       % pass it through the sigmoid function 
       fx = sigmoid(fx);
    case 'genHill'
       % generalized hill function function 
       %ff(ff==0)=eps;
       %fx = W*log(ff) + W0(:, ones(1, SizF)); %repmat(W0,[1 SizF]);
       fx = W*log(ff+1e-100) + W0(:, ones(1, SizF)); %repmat(W0,[1 SizF]);
       %W0 = sum(W,2).*log(Par.Gammas(J)); 
       %
       fx = 1 ./ (1 + exp(-fx));
       %fx = sigmoid(fx); 
    case 'lin'
       % 
       fx = W*ff + W0(:, ones(1, SizF)); %repmat(W0,[1 SizF]); 
       %
    case 'michMenten'
       % 
       fx = michMenten(ff, W, LikParams.Net_X(J,:));
       %
    %case 'michMentenAct'
    %   xp = repmat(ff,[NumGenes 1]); 
    %   fx = xp./(xp + repmat(exp(-W0),[1 SizF]));
    %case 'michMentenRepres'
    %   xp = repmat(ff,[NumGenes 1]); 
    %   fx = 1./(xp + repmat(exp(-W0),[1 SizF])); 
end

% binarize the outputs if necessary 
if LikParams.jointActBin == 1
    fx = round(fx);
end

                                                                                                                                                                                                                                                                                                                                             activFuncts/jointactFunc2.m                                                                         0000700 0003466 0000024 00000003023 11267363056 015750  0                                                                                                    ustar   mtitsias                        games                                                                                                                                                                                                                  function fx = jointactFunc2(Par,ff,J)
%
%
%
%

NumGenes = size(J(:),1);

W = Par.W(J,:);
Tausindex = Par.Tausindex(J); 

W0 = Par.W0(J); 
if strcmp(Par.jointAct,'michMenten')
Net_X = Par.Net_X(J,:);
end
switch Par.jointAct
    case 'lin'
       % 
       % this account for delays in gene expression
       for m=1:size(J,2)
         fx(m,:) = W(m,:)*ff(:,Tausindex(m):Tausindex(m)+Par.sizTime-1);
       end
       fx = fx + repmat(W0,[1 Par.sizTime]);
       %
    case 'sigmoid'
       %
       % this account for delayes in gene expression
       for m=1:size(J,2)
         fx(m,:) = W(m,:)*ff(:,Tausindex(m):Tausindex(m)+Par.sizTime-1);
       end
       fx = fx + repmat(W0,[1 Par.sizTime]); 
       % pass it through the sigmoid function 
       fx = sigmoid(fx);
    case 'michMenten'
       ops.Net_X = Net_X; 
       ops.sizTime = Par.sizTime; 
       ops.Tausindex = Tausindex;
       fx = michMenten(ff, W, ops);
    case 'genHill'
       % generalized hill function function  
       for m=1:size(J,2)
         fx(m,:) = W(m,:)*log(ff(:,Tausindex(m):Tausindex(m)+Par.sizTime-1));
       end
       %W0 = sum(W,2).*log(Par.Gammas(J)); 
       %
       fx = fx + repmat(W0,[1 Par.sizTime]); 
       fx = sigmoid(fx); 
    %case 'michMentenAct'
    %   xp = repmat(ff,[NumGenes 1]); 
    %   fx = xp./(xp + repmat(exp(-W0),[1 SizF]));
    %case 'michMentenRepres'
    %   xp = repmat(ff,[NumGenes 1]); 
    %   fx = 1./(xp + repmat(exp(-W0),[1 SizF])); 
end

% binarize the outputs if necessary 
if Par.jointActBin == 1
    fx = round(fx);
end

                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             activFuncts/lin.m                                                                                   0000700 0003466 0000024 00000000040 11267363056 014015  0                                                                                                    ustar   mtitsias                        games                                                                                                                                                                                                                  function fx = lin(x)
%

fx = x;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                activFuncts/logOnePlusExp.m                                                                         0000700 0003466 0000024 00000000070 11267363056 016002  0                                                                                                    ustar   mtitsias                        games                                                                                                                                                                                                                  function fx = logOnePlusExp(x)
%

fx = log(1 + exp(x));
                                                                                                                                                                                                                                                                                                                                                                                                                                                                        activFuncts/loglogOnePlusExp.m                                                                      0000700 0003466 0000024 00000000100 11267363056 016476  0                                                                                                    ustar   mtitsias                        games                                                                                                                                                                                                                  function fx = loglogOnePlusExp(x)
%

fx = log(log(1 + exp(x)));
                                                                                                                                                                                                                                                                                                                                                                                                                                                                activFuncts/michMenten.m                                                                            0000700 0003466 0000024 00000003740 11301016752 015317  0                                                                                                    ustar   mtitsias                        games                                                                                                                                                                                                                  function fx = michMenten(ff, W, NetX)
% function fx = michMenten(x, W, netX)
% 
% Description: Impements the Michaelis-Menten multiple TF activation 
%              The single TF Michaelis-Menten equations are obtained as special cases
%
% Inputs: 
%    -- ff:This is NumoFTs x NumOfTimes size matrix that contains  
%           that stores TFs values.
%    -- W: This a NumOfGenes x NumOfTFs matrix that contains 
%           connectivity weights that take NON-NEGATIVE values. 
%           Each element W(i,j) stores the interaction weight 
%           between the i TF and the j gene. 
%    -- netX: This a NumOfGenes x NumOfTFs matrix that describes the 
%           network connectivity. Each element netX(i,j) can take three
%           values, (-1, 0, 1). When netX(j,i) = -1, then the i TF acts 
%           as repressor for the j gene. When netX(j,i) = 0, the i TF does not 
%           interact with the j gene, and when netX(j,i) = 1, the i TF 
%           activates the j gene. 
%
% Output: fx: This a  NumOfGenes x NumOfTimes matrix that computes the 
%          mulitple Tf Michaelis-Menten function:
%            
%             fx(j,t) =  (\sum_{i in R} w_ji+\sum_{i in A} w_ji f_i(t))/...
%                         (1 + \sum_{i in R} w_ji f(t)  + \sum_{i in A} w_ji f_i(t))
%
%          where 'A' is the set of TF activator for the j gene and
%          'R' is the set of repressors for the j gene. 
%                        
% Notes: When we use only one TF, then we obtain the M-M single TF equation. 
%        In that case the netX is NumOfGenes x 1 vector with -1s and 1s.
%        E.g. for one gene j netX(j) = -1, the equation becomes the repressor
%        M-M given by
%             fx(j,t) =  w_j /(1 + w_j f(t))
%        where you can extract the M-M constant form gamma_j =1/w_j. 


SizF = size(ff,2);

% find all repressors 
R = (NetX==-1);
% find all activators 
A = (NetX==1);
xp = (W.*A)*ff;
xpR = W.*R;

% combine the TFs
fx = (xp + repmat(sum(xpR,2),[1 SizF]))./(1 + xp + xpR*ff);


                                activFuncts/michMentenAct.m                                                                         0000700 0003466 0000024 00000000076 11267363056 015763  0                                                                                                    ustar   mtitsias                        games                                                                                                                                                                                                                  function fx = michMentenAct(x,gam)
%
%
%

fx = x./(x + gam);

                                                                                                                                                                                                                                                                                                                                                                                                                                                                  activFuncts/michMentenRepres.m                                                                      0000700 0003466 0000024 00000000100 11267363056 016500  0                                                                                                    ustar   mtitsias                        games                                                                                                                                                                                                                  function fx = michMentenRepres(x,gam)
%
%
%

fx = 1./(x + gam);
                                                                                                                                                                                                                                                                                                                                                                                                                                                                activFuncts/sigmoid.m                                                                               0000700 0003466 0000024 00000000064 11267363056 014674  0                                                                                                    ustar   mtitsias                        games                                                                                                                                                                                                                  function fx = sigmoid(x)
%
%

fx = 1./(1 + exp(-x));                                                                                                                                                                                                                                                                                                                                                                                                                                                                            activFuncts/singleactFunc.m                                                                         0000700 0003466 0000024 00000000112 11267363056 016020  0                                                                                                    ustar   mtitsias                        games                                                                                                                                                                                                                  function F = singleactFunc(singleAct,F)
%
%
%
%

F = feval(singleAct,F);

                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      