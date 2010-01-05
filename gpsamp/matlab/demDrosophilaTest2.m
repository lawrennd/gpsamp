
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



