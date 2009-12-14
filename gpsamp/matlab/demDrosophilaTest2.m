

dataName = 'drosophila_dataTest';
expNo = 1;
storeRes = 0;
printPlot = 0;


% This should go inside the loadDatasets function later
%%%%%%%%%%%%%%  Load the test genes data  %%%%%%%%%%%%%%%% 
%load datasets/drosophila_data.mat;
% load traning samples 
load /usr/local/michalis/mlprojects/gpsamp/matlab/drosTrainTotal;
load /usr/local/michalis/mlprojects/gpsamp/matlab/trainGenes;

TestGenes = Genes;
TestGenesVar = GenesVar;
TimesG = 0:11;
numTFs = 5;
%%%%%%%%%%%%%%  Load data  %%%%%%%%%%%%%%%% 

% model options
options = gpmtfOptions(TestGenes,numTFs); 
%options.constraints.replicas = 'coupled'; 

options.tauMax = 0; % no delays
% define the dense discretized grid in the time axis for the TF latent functions 
[options, TimesF] = gpmtfDiscretize(TimesG, options); 

% CREATE the model
modelTest = gpmtfCreate(TestGenes, TestGenesVar, [], [], TimesG, TimesF, options);

mcmcoptions = mcmcOptions('controlPnts'); 
mcmcoptions.adapt.T = 50; 
mcmcoptions.adapt.Burnin = 50;
mcmcoptions.train.StoreEvery = 20;

% precompute the TFs
if strcmp(modelTest.constraints.replicas,'free')
    for cnt=1:size(samples.F,2)
        modelTest.Likelihood.kineticsTF = samples.kineticsTF(:,:,cnt);
        for r=1:modelTest.Likelihood.numReplicas
            TFs{cnt}(:,:,r) = gpmtfComputeTF(modelTest.Likelihood, samples.F{cnt}(:,:,r), 1:numTFs);
        end
    end
else 
    % replicas are coupled 
    for cnt=1:size(samples.F,2)
        modelTest.Likelihood.kineticsTF = samples.kineticsTF(:,:,cnt);
        TFs{cnt} = gpmtfComputeTF(modelTest.Likelihood, samples.F{cnt}, 1:numTFs);
    end
    %
end

   
% compute all distances between samples (useful for sampling using the
% geometric distribution)
if strcmp(modelTest.constraints.replicas,'free')
    for j=1:modelTest.Likelihood.numTFs
    for r=1:modelTest.Likelihood.numReplicas     
        % collect all TFs  of the r replica in matrix 
        MatTF = zeros(size(samples.F,2), size(modelTest.Likelihood.TimesF,2)); 
        for cnt=1:size(samples.F,2)
            MatTF(cnt,:)  = TFs{cnt}(j,:,r); 
        end

        % compute the similarity matrix
        dd = dist2(MatTF,MatTF);
        for cnt=1:size(samples.F,2)
           [ok Ind] = sort(dd(cnt,:));
           simMat{j,r}(cnt,:) = uint16(Ind(2:end));
        end
    end
    end
    %
else 
    %
    for j=1:modelTest.Likelihood.numTFs
    MatTF = zeros(size(samples.F,2), size(modelTest.Likelihood.TimesF,2)); 
    for cnt=1:size(samples.F,2)
        MatTF(cnt,:)  = TFs{cnt}(j,:); 
    end

    % compute the similarity matrix
    dd = dist2(MatTF,MatTF);
    for cnt=1:size(samples.F,2)
       [ok Ind] = sort(dd(cnt,:));
       simMat{j}(cnt,:) = uint16(Ind(2:end));
    end
    end
    %
end


% parameter in the geometric distribution
modelTest.geo = 0.1;

[modelTest PropDist samplesTest accRates] = gpmtfTestGenesAdapt2(modelTest, TFs, simMat, mcmcoptions.adapt);
% training/sampling phase
[modelTest PropDist samplesTest accRates] = gpmtfTestGenesSample2(modelTest, TFs, simMat, PropDist, mcmcoptions.train);


