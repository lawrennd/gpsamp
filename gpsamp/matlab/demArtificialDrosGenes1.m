
numGenes = 1030; 
numConds = 2;

% more than one conditions simply create more replicas
numReplicas = numConds;

%%%%%%%%%%%%%%  Load the test genes data  %%%%%%%%%%%%%%%% 
load datasets/drosophila_data;
load datasets/trainScaleDros;
load datasets/testset;
load drosTrainTotal;
load CollectKin.mat;
    
% generate random genes
numTFs = 4;
% the last TF is an outlier
TFset = [1 2 3 4];
%TimesG = [0 1 2 3 4 6 10 11];
TimesG = [0 1 2 3 5 7 9 11 14 18]
sigma2 = 0.025;

% model options 
options = gpmtfOptions(ones(1,12,3),numTFs);
options.jointAct = 'sigmoid';
options.spikePriorW = 'yes';
options.noiseModel = {'pumaWhite'};
options.constraints.spaceW = 'positive';

% prior probablity for each interaction weight to be around zero 
% for each TF
%TFpis = [0.25 0.25 0.25 0.25 0.25];
TFpis = [0.5 0.5 0.5 0.5];
options.spikepriors = 1 - TFpis;

options.tauMax = 0; % no delays
% define the dense discretized grid in the time axis for the TF latent functions 
[options, TimesF] = gpmtfDiscretize(TimesG, options);  
modelTest = gpmtfCreate(ones(1, size(TimesG,2), numReplicas), ones(1, size(TimesG,2), numReplicas), [], [], TimesG, TimesF, options);

% ground GP function for the different conditions 
load mRNAToyGPF.mat;

%douplicate the outlier for both conditions
outlF = reshape([outlierF; outlierF]', [1 size(outlierF,2) 2]);
F = [F; outlF];
% precompute the TFs
modelTest.Likelihood.kineticsTF = [0.99461 1; 0.9457 1; 0.640 1;  1.2 1];
for r=1:numConds
     TFs(:,:,r) = gpmtfComputeTFODE(modelTest.Likelihood, F(TFset,:,r), 1:numTFs);  
end
TFs = log(TFs + 1e-100);

%% precompute the TFs
%for cnt=1:size(samples.F,2)
%    modelTest.Likelihood.kineticsTF = samples.kineticsTF(TFset,:,cnt);
%    if cnt == 1
%       for j=1:numTFs 
%          tmpTF = log(1+ exp(samples.F{cnt}(j,:,:)) ); 
%          scTF(j) = max(tmpTF(:));
%       end
%    end
%    %
%    for j=1:numTFs
%       tmpTF = log( 1 + exp(samples.F{cnt}(j,:,:)) ); 
%       tmpTF = (10/scTF(j))*tmpTF; 
%       % inverse logonePlusexp
%       samples.F{cnt}(j,:,:) = log( exp(tmpTF)  - 1);
%    end
%    for r=1:numReplicas 
%       TFs{cnt}(:,:,r) = gpmtfComputeTFODE(modelTest.Likelihood, samples.F{cnt}(TFset,:,r), 1:numTFs);
%    end
%    TFs{cnt} = log(TFs{cnt} + 1e-100); 
%end
    
% randomly select TFs profiles for all conditions 
%ch = randperm(size(samples.F,2));
%for c=1:numConds 
%    GPF{c} = samples.F{ch(c)};
%    TFind(c) = ch(c);
%end


% generate mRNA for the TF-Genes
for r=1:numReplicas
    GnsTF(:,:,r) = singleactFunc(modelTest.Likelihood.singleAct, F(:, :, r));
end
GenesTF = GnsTF(:,modelTest.Likelihood.comInds,:);
Vars = sigma2*GenesTF;
GenesTF = GenesTF + sqrt(Vars).*randn(size( GenesTF ));
GenesTFVar = Vars;
     
Genes = zeros(numGenes, size(TimesG,2), numReplicas);
GenesFree = zeros(numGenes, size(TimesG,2), numReplicas);
GenesVar = ones(numGenes, size(TimesG,2), numReplicas);
    
LikParams = modelTest.Likelihood; 
LikParams.TF = TFs;
% 
for n=1:numGenes

    while 1
       ch = rand(1,numTFs); 
       ch(ch<=TFpis)=1;
       ch(ch~=1)=0;
       if sum(ch(:)) == 0
          % at least one TF should be a regulator 
          ok = randperm(numTFs);
          ch(ok(1)) = 1; 
       end

       Net(n,:) = ch;
       % in the first 500  genes do not allow the outlier TF
       if n <= 530
          Net(n,end) = 0;
       end
              
       % do not allowe many genes with very high initial value
       while 1 
         ch1 = randperm(6177);
         Kinetics = Kin(ch1(1), :);
         if Kinetics(4) < 3
             break; 
         elseif (rand < 0.2)
             break;
         end
         
       end

       LikParams.kinetics = Kinetics;
       
       while 1
         k = (randn(1,numTFs) + 0.5);
         W(n,:) = k.*Net(n,:); 
         % do not allow genes expression to have very small interections
         % weights
         if min(abs(k)) > 0.1
             break; 
         end
       end
       W0(n) = options.constraints.W0.*randn;
       LikParams.W = W(n,:);
       LikParams.W0 = W0(n);
       
       % do not allow very noisy genes having samll expression or very
       % uniform expressions
       for r=1:numReplicas
           Gns(:,:,r) = gpmtfComputeGeneODE(LikParams, zeros(numTFs,111), r, 1);
       end   
       if (max(Gns(:)) > 2 & max(Gns(:)) < 10)  & (max(Gns(:))-min(Gns(:)))>2 
           break;
       end
    end
    
    Genes(n,:,:) = Gns(1, LikParams.comInds, :); 
    Vars = sigma2*Genes(n,:,:); 
    GenesFree(n,:,:) = Genes(n,:,:);
    Genes(n,:,:) = Genes(n,:,:) + sqrt(Vars).*randn(size(Genes(n,:,:)));
    GenesVar(n,:,:) = Vars(1,:,:);
    %
    GrTruth{n}.kinetics = LikParams.kinetics;
    GrTruth{n}.indParam = ch1(1);
    GrTruth{n}.W = LikParams.W;
    GrTruth{n}.W0 = LikParams.W0;
%
end

    
for n=1:numGenes 
   subplot(1,2,1);  
   plot(TimesG, Genes(n, :, 1), 'b');
   hold on;
   plot(TimesG, Genes(n, :, 1), '+b');
   hold on;
   plot(TimesG, Genes(n, :, 2), 'r');
   plot(TimesG, Genes(n, :, 2), '+r');
   hold off;
   
   n
   
   subplot(1,2,2);  
   plot(TimesG, GenesFree(n, :, 1), 'b');
   hold on;
   plot(TimesG, GenesFree(n, :, 1), '+b');
   hold on;
   plot(TimesG, GenesFree(n, :, 2), 'r');
   plot(TimesG, GenesFree(n, :, 2), '+r');
   hold off;
   
   Net(n,:)
   GrTruth{n}
   n
   pause;
end
