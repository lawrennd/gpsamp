function plotCreate(model, samples, Genes, TimesG, TimesF, demdata, printResults, GeneVars)
%function plotCreate(model, samples, Genes, TimesG, TimesF, demdata, printResults, GeneVars)
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

if strcmp(model.Likelihood.TFjointAct,'sigmoid') 
   model.Likelihood.TFsingleAct = 'exp';
end

ok = date;
fileName = [demdata 'MCMC' ok model.Likelihood.TFsingleAct model.Likelihood.TFjointAct]; 

NumOfTFs = model.Likelihood.NumOfTFs;

if strcmp(model.Constraints.sigma2,'free') 
  % plots sigma2s and the lengghscales
  for j=1:model.Likelihood.NumOfGenes
  sigma2j = squeeze(samples.sigma2(j,1,:));
  figure;
  hist(sigma2j,100);   
  if printResults
       print('-depsc', ['./results/' fileName 'Sigma2']);
  end
  titlestring = 'Observation variance: ';
  titlestring = [titlestring, num2str(j)]; 
  titlestring = [titlestring, ' gene'];
  title(titlestring,'fontsize', 20);
  end
end    

% plot the lengthscales
for j=1:NumOfTFs
    figure;
    hist(squeeze(exp(2*samples.logthetas(j,1,:))),100); 
    %title('Lengthscale','fontsize', 20);
    if printResults
      print('-depsc', ['./results/' fileName 'LengthSc' 'TF' num2str(j)]);
    end
    titlestring = 'Lengthscale: ';
    titlestring = [titlestring, num2str(j)]; 
    titlestring = [titlestring, ' TF'];
    title(titlestring,'fontsize', 20);
end

NumOfGenes = model.Likelihood.NumOfGenes;
order = 1:NumOfGenes;
NumOfReplicas = model.Likelihood.NumOfReplicas;
NumOfSamples = size(samples.F,2);
SizF = size(samples.F{1},2);
TimesF = TimesF(:); 
for r=1:NumOfReplicas
  % 
  for j=1:NumOfTFs
     % 
     FF = zeros(NumOfSamples,SizF);    
     for t=1:NumOfSamples
         FF(t,:) = samples.F{t}(j,:,r);
     end
    
     FF = feval(model.Likelihood.TFsingleAct,FF);
     mu = mean(FF)';
     stds1 = sqrt(var(FF))';
     stds2 = sqrt(var(FF))';
     if strcmp(model.Likelihood.TFsingleAct,'exp')==1
       mu = median(FF)';
       stds1 = (prctile(FF,95,1)'-mu)/2;
       stds2 = (mu-prctile(FF,5,1)')/2;
     end
     figure
     plot(TimesF,mu,'b','lineWidth',3);
     hold on;
     fillColor = [0.7 0.7 0.7];
     %fillColor = [0.8 0.8 0.8];  % for the paper
     fill([TimesF; TimesF(end:-1:1)], [mu; mu(end:-1:1)]...
            + 2*[stds1; -stds2(end:-1:1)], fillColor,'EdgeColor',fillColor);
     plot(TimesF,mu,'b','lineWidth',3);
     
     axis([TimesG(1) TimesG(end)+0.1 0 max(mu+2*stds1)+0.1]);
     
     
     
     % plot the ground truth if exist
     if strcmp(model.GroundTr,'yes')==1
     FFgt = feval(model.Likelihood.TFsingleAct,model.FF(j,:,r));
     plot(TimesF,FFgt,'r','lineWidth',3);
     end
     
     if printResults
      print('-depsc', ['./results/' fileName 'Replica' num2str(r) 'TFprof' num2str(j)]);
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


if 0

% plot statistics of the joint activation 
for r=1:NumOfReplicas
  % 
  for j=1:NumOfGenes
     % 
     FF = zeros(NumOfSamples,SizF);    
     for t=1:NumOfSamples
         FF(t,:) = samples.F{t}(1,:,r);
         LikParams.NumOfGenes = 1;
         %% Gene - TFs interaction weights 
         LikParams.W = samples.Weights(j,1,t);
         %% gene bias term in the regulatory part
         LikParams.W0 = samples.Weights0(j,t);
         LikParams.TFjointAct = model.Likelihood.TFjointAct;
         LikParams.TFsingleAct = model.Likelihood.TFsingleAct;
         LikParams.TFjointActBin = model.Likelihood.TFjointActBin;
         FF(t,:) = TFactivFun(LikParams,FF(t,:));
     end
     mu = mean(FF)';
     mu
     stds1 = sqrt(var(FF))';
     stds2 = sqrt(var(FF))';
     figure
     plot(TimesF,mu,'b','lineWidth',3);
     hold on;
     fillColor = [0.7 0.7 0.7];
     %fillColor = [0.8 0.8 0.8];  % for the paper
     fill([TimesF; TimesF(end:-1:1)], [mu; mu(end:-1:1)]...
            + 2*[stds1; -stds2(end:-1:1)], fillColor,'EdgeColor',fillColor);
     plot(TimesF,mu,'b','lineWidth',3);
     
     % plot the ground truth if exist
     if strcmp(model.GroundTr,'yes')==1
     FFgt = feval(model.Likelihood.TFsingleAct,model.FF(j,:,r));
     plot(TimesF,FFgt,'r','lineWidth',3);
     end
     
     %
  end
  %
end

end

%TimesF = TimesF(21:end);

% plot predicted gene expressions 
for r=1:NumOfReplicas
  %  
  for j=1:NumOfGenes
     % 
     GG = zeros(NumOfSamples,SizF);    
     for t=1:NumOfSamples
         GG(t,:) = samples.predGenes{t}(j,:,r);
     end
     
     mu = mean(GG)';
     %mu = mu(21:end);
     mu = mu(1:2:end);
     stds = sqrt(var(GG))';
     %stds = stds(21:end);
     stds = stds(1:2:end);
    
     TF = TimesF(1:2:end);
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
     if nargin == 8
     errorbar(TimesG,  Genes(j,:,r), 2*sqrt(GeneVars(j,:,r)), 'rx','lineWidth', 1.5);
     end
     axis([min(TimesG(:))-0.1 max(TimesG(:))+0.1 0.95*min(min(Genes(j,:,r))) 1.05*max(max(Genes(j,:,r)))]);
     
     if printResults
      print('-depsc', ['./results/' fileName 'Replica' num2str(r) 'GeneExp' num2str(j)]);
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
if strcmp(model.GroundTr,'yes')==1
bar([modelB(order); model.GtrKinetics(:,1)']', 0.7); colormap([0.9 0.9 0.9; 0 0 0]);
else
bar(modelB(order)', 0.7); colormap([0.9 0.9 0.9]);
end
hold on;
errorbar([1:NumOfGenes]-0.14, modelB(order), modelB(order)-stdBB1(order), stdBB2(order)-modelB(order),'.'); 
%errorbar([1:NumOfGenes], modelB(order), modelB(order)-stdBB1(order), stdBB2(order)-modelB(order),'.'); 
title('Basal rates','fontsize', 20);

if printResults
      print('-depsc', ['./results/' fileName 'Basal']);
end
     

% Plot the sensitivities.
figure;
if strcmp(model.GroundTr,'yes')==1
bar([modelS(order); model.GtrKinetics(:,3)']', 0.7); colormap([0.9 0.9 0.9; 0 0 0]);
else
bar(modelS(order)', 0.7); colormap([0.9 0.9 0.9]);
end
hold on;
errorbar([1:NumOfGenes]-0.14, modelS(order), modelS(order)-stdSS1(order), stdSS2(order)-modelS(order),'.');
%errorbar([1:NumOfGenes], modelS(order), modelS(order)-stdSS1(order), stdSS2(order)-modelS(order),'.');
title('Sensitivities','fontsize', 20);

if printResults
      print('-depsc', ['./results/' fileName 'Sensitivity']);
end

figure;
% plot degradation rates
if strcmp(model.GroundTr,'yes')==1
bar([modelD(order); model.GtrKinetics(:,2)']', 0.7); colormap([0.9 0.9 0.9; 0 0 0]);
else
bar(modelD(order)', 0.7); colormap([0.9 0.9 0.9]);
end
hold on;
errorbar([1:NumOfGenes]-0.14, modelD(order), modelD(order)-stdDD1(order), stdDD2(order)-modelD(order),'.');
%errorbar([1:NumOfGenes], modelD(order), modelD(order)-stdDD1(order), stdDD2(order)-modelD(order),'.');
title('Decays','fontsize', 20);

if printResults
      print('-depsc', ['./results/' fileName 'Decay']);
end

figure;
% plot initial conditions 
if strcmp(model.GroundTr,'yes')==1
bar([modelA(order); model.GtrKinetics(:,4)']', 0.7); colormap([0.9 0.9 0.9; 0 0 0]);
else
bar(modelA(order)', 0.7); colormap([0.9 0.9 0.9]);
end
hold on;
errorbar([1:NumOfGenes]-0.14, modelA(order), modelA(order)-stdAA1(order), stdAA2(order)-modelA(order),'.');
%errorbar([1:NumOfGenes], modelA(order), modelA(order)-stdAA1(order), stdAA2(order)-modelA(order),'.');
title('Initial conditions','fontsize', 20);

if printResults
      print('-depsc', ['./results/' fileName 'InitCond']);
end

for j=1:NumOfTFs
W1 = squeeze(samples.Weights(:,j,:))';
modelW1 = mean(W1,1);
stdW1_1 = sqrt(var(W1));
stdW1_2 = sqrt(var(W1));
% Plot first basal transcription rates.
figure;
if strcmp(model.GroundTr,'yes')==1
bar([modelW1(order); model.GtrW(:,j)']', 0.7); colormap([0.9 0.9 0.9; 0 0 0]);
else 
bar(modelW1(order)', 0.7); colormap([0.9 0.9 0.9]);
end
hold on;
errorbar([1:NumOfGenes]-0.14, modelW1(order), modelW1(order)-stdW1_1(order), stdW1_2(order)-modelW1(order),'.'); 
%errorbar([1:NumOfGenes], modelB(order), modelB(order)-stdBB1(order), stdBB2(order)-modelB(order),'.'); 
%title(j,'fontsize', 20);
if printResults
      print('-depsc', ['./results/' fileName 'IntWeights' 'TF' num2str(j)]);
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
if strcmp(model.GroundTr,'yes')==1
bar([modelW0(order); model.GtrW0']', 0.7); colormap([0.9 0.9 0.9; 0 0 0]);
else
bar(modelW0(order)', 0.7); colormap([0.9 0.9 0.9]);
end
hold on;
errorbar([1:NumOfGenes]-0.14, modelW0(order), modelW0(order)-stdW0_1(order), stdW0_2(order)-modelW0(order),'.');
%errorbar([1:NumOfGenes], modelA(order), modelA(order)-stdAA1(order), stdAA2(order)-modelA(order),'.');
%title('W0','fontsize', 20);
if printResults
      print('-depsc', ['./results/' fileName 'IntBias']);
end
title('Interaction biases','fontsize', 20);

