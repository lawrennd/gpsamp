function gpmtfPlot(model, samples, demdata, TFname, printResults)
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

FONTSIZE=8;

%if strcmp(model.Likelihood.TFjointAct,'sigmoid') 
%   model.Likelihood.TFsingleAct = 'exp';
%end
dirr = '/usr/local/michalis/mlprojects/gpsamp/tex/diagrams/';

Genes = model.Likelihood.Genes;

if model.Likelihood.noiseModel.active(1) == 1
    GeneVars = model.Likelihood.noiseModel.pumaSigma2;
end
 
TimesG = model.Likelihood.TimesG; 
TimesF = model.Likelihood.TimesF;
if isfield(model.Likelihood,'GenesTF')
    GenesTF = model.Likelihood.GenesTF;
    if model.Likelihood.noiseModel.active(1) == 1
        GeneTFVars = model.Likelihood.noiseModel.pumaSigma2_TF;
    end
end

TimesFF = TimesF(model.Likelihood.startTime:end);
ok = date;
fileName = [demdata 'MCMC' ok model.Likelihood.singleAct model.Likelihood.jointAct]; 

NumOfTFs = model.Likelihood.numTFs;

% if 0

% if model.Likelihood.noiseModel.active(1) == 1
%   % plots sigma2s and the lengghscales
%   for j=1:model.Likelihood.numGenes
%   sigma2j = squeeze(samples.sigma2(j,:));
%   figure;
%   hist(sigma2j,100);   
%   if printResults
%        print('-depsc', [dirr fileName 'Sigma2']);
%   end
%   titlestring = 'Observation variance: ';
%   titlestring = [titlestring, num2str(j)]; 
%   titlestring = [titlestring, ' gene'];
%   title(titlestring,'fontsize', FONTSIZE);
%   end
% end    

% % if isfield(model.Likelihood,'sigmasTF')
% % if 0
% %   % plots sigma2s and the lengghscales
% %   for j=1:model.Likelihood.numTFs
% %   sigma2j = squeeze(samples.sigmasTF(j,1,:));
% %   figure;
% %   hist(sigma2j,100);   
% %   if printResults
% %        print('-depsc', [dirr fileName 'SigmasTF']);
% %   end
% %   titlestring = 'Observation variance-TF Genes: ';
% %   titlestring = [titlestring, num2str(j)]; 
% %   titlestring = [titlestring, ' TF'];
% %   title(titlestring,'fontsize', FONTSIZE);
% %   end
% % end  
% % end

% % plot the lengthscales
% for j=1:NumOfTFs
%     figure;
%     hist(squeeze(exp(2*samples.lengthScale(j,:))),30); 
%     %title('Lengthscale','fontsize', FONTSIZE);
%     if printResults
%       print('-depsc', [dirr fileName 'LengthSc' 'TF' num2str(j)]);
%     end
%     titlestring = 'Lengthscale: ';
%     titlestring = [titlestring, num2str(j)]; 
%     titlestring = [titlestring, ' TF'];
%     title(titlestring,'fontsize', FONTSIZE);
% end


% end

NumOfGenes = model.Likelihood.numGenes;
order = 1:NumOfGenes;
NumOfReplicas = model.Likelihood.numReplicas;
NumOfSamples = size(samples.F,2);
SizF = size(samples.F{1},2);
TimesF = TimesF(:); 
figure;
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
     
     %figure
     subplot(NumOfTFs, NumOfReplicas, r + (j-1)*NumOfReplicas);
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
     
     titlestring = 'Profile: ';
     titlestring = [titlestring, num2str(r)]; 
     titlestring = [titlestring, ' replica, '];
     titlestring = [titlestring, TFname{j}];
     %titlestring = [titlestring, ' TF'];
     title(titlestring,'fontsize', FONTSIZE);
     %
  end
  %
end
if printResults
  print('-depsc', [dirr fileName 'Replica' num2str(r) 'TF' TFname(j)]);
end 

figure;
% predicted TF-Gene expressions 
if isfield(model.Likelihood,'GenesTF')
for r=1:NumOfReplicas
  %  
  for j=1:NumOfTFs
     % 
     GG = zeros(NumOfSamples, size(model.Likelihood.TimesF,2));    
     for t=1:NumOfSamples
         PredGenesTF = singleactFunc(model.Likelihood.singleAct, samples.F{t}(j,:,r));
         GG(t,:) = PredGenesTF;
     end
     
     mu = mean(GG)';
     %mu = mu(1:2:end);
     stds = sqrt(var(GG))';
     %stds = stds(1:2:end);
    
     TF = TimesF; % TimesFF(1:2:end)';
     %figure
     subplot(NumOfTFs, NumOfReplicas, r + (j-1)*NumOfReplicas);
     plot(TF,mu,'b','lineWidth',r);
     hold on;
     fillColor = [0.7 0.7 0.7];
     %fillColor = [0.8 0.8 0.8];  % for the paper
     fill([TF; TF(end:-1:1)], [mu; mu(end:-1:1)]...
            + 2*[stds; -stds(end:-1:1)], fillColor,'EdgeColor',fillColor);
     plot(TF,mu,'b','lineWidth',3);
   
     plot(TimesG,GenesTF(j,:,r),'rx','markersize', 14','lineWidth', 2);
     if model.Likelihood.noiseModel.active(1) == 1
     errorbar(TimesG,  GenesTF(j,:,r), 2*sqrt(GeneTFVars(j,:,r)), 'rx','lineWidth', 1.5);
     end
     
     axis([min(TimesG(:))-0.1 max(TimesG(:))+0.1  0.95*min([GenesTF(j,:,r), mu' - stds'])  1.05*max([GenesTF(j,:,r), mu' - stds'])]);
     
     titlestring = 'Expressions: ';
     titlestring = [titlestring, num2str(r)]; 
     titlestring = [titlestring, ' replica, '];
     %titlestring = [titlestring, num2str(j)];
     titlestring = [titlestring, TFname{j}];
     title(titlestring,'fontsize', FONTSIZE);
     %
  end
  %
end    
if printResults
  print('-depsc', [dirr fileName 'Replica' num2str(r) 'GeneTFExp' num2str(j)]);
end
end


% plot predicted gene expressions 
for j=1:NumOfGenes
  %  
  if mod(j,5) == 1,
    figure;
  end
  for r=1:NumOfReplicas
     % 
     GG = zeros(NumOfSamples,model.Likelihood.sizTime);    
     for t=1:NumOfSamples
         %
         LikParams = model.Likelihood; 
         LikParams.kinetics = samples.kinetics(:,:,t);
         LikParams.kineticsTF = samples.kineticsTF(:,:,t);
         LikParams.W = samples.W(:,:,t);
         LikParams.W0 = samples.W0(:,t);
         %LikParams.TF = TFs{samples.TFindex(t)};
         predgen = gpmtfComputeGeneODE(LikParams, samples.F{t}(:,:,r), r, j);
         GG(t,:) = predgen;
         %predgen(model.Likelihood.comInds)
         %Genes(j,:,r)
         %pause
         %
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
     %figure
     subplot(5, NumOfReplicas, r + mod(j-1, 5)*NumOfReplicas);
     plot(TF,mu,'b','lineWidth',r);
     hold on;
     fillColor = [0.7 0.7 0.7];
     %fillColor = [0.8 0.8 0.8];  % for the paper
     fill([TF; TF(end:-1:1)], [mu; mu(end:-1:1)]...
            + 2*[stds; -stds(end:-1:1)], fillColor,'EdgeColor',fillColor);
     plot(TF,mu,'b','lineWidth',3);
   
     plot(TimesG,Genes(j,:,r),'rx','markersize', 14','lineWidth', 2);
     if model.Likelihood.noiseModel.active(1) == 1
     errorbar(TimesG,  Genes(j,:,r), 2*sqrt(GeneVars(j,:,r)), 'rx','lineWidth', 1.5);
     end
     axis([min(TimesG(:))-0.1 max(TimesG(:))+0.1  0.95*min([Genes(j,:,r), mu' - stds'])  1.05*max([Genes(j,:,r), mu' - stds'])  ]);
     
     if printResults && r==NumOfReplicas && (mod(j,5)==0 || j==NumOfGenes),
      print('-depsc', [dirr fileName 'Replica' num2str(r) 'GeneExp' num2str(j)]);
     end
     titlestring = 'Expressions: ';
     titlestring = [titlestring, num2str(r)]; 
     titlestring = [titlestring, ' replica, '];
     titlestring = [titlestring, num2str(j)];
     titlestring = [titlestring, ' gene'];
     title(titlestring,'fontsize', FONTSIZE);
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
if isfield(model,'groundtr') == 1
bar([modelB(order); model.groundtr.kinetics(:,1)']', 0.7); colormap([0.9 0.9 0.9; 0 0 0]);
else
bar(modelB(order)', 0.7); colormap([0.9 0.9 0.9]);
end
hold on;
errorbar([1:NumOfGenes]-0.14, modelB(order), modelB(order)-stdBB1(order), stdBB2(order)-modelB(order),'.'); 
%errorbar([1:NumOfGenes], modelB(order), modelB(order)-stdBB1(order), stdBB2(order)-modelB(order),'.'); 
title('Basal rates','fontsize', FONTSIZE);

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
title('Sensitivities','fontsize', FONTSIZE);

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
title('Decays','fontsize', FONTSIZE);

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
title('Initial conditions','fontsize', FONTSIZE);

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
title('TF-Genes Decays','fontsize', FONTSIZE);

if printResults
      print('-depsc', [dirr fileName 'TFGeneDecay']);
end

% Plot the sensitivities.
% figure;
% if isfield(model,'groundtr') == 1
% bar([modelS(orderTF); model.groundtr.kineticsTF(:,2)']', 0.7); colormap([0.9 0.9 0.9; 0 0 0]);
% else
% bar(modelS(orderTF)', 0.7); colormap([0.9 0.9 0.9]);
% end
% hold on;
% errorbar([1:NumOfTFs]-0.14, modelS(orderTF), modelS(orderTF)-stdSS1(orderTF), stdSS2(orderTF)-modelS(orderTF),'.');
% %errorbar([1:NumOfGenes], modelS(order), modelS(order)-stdSS1(order), stdSS2(order)-modelS(order),'.');
% title('TF-Genes Sensitivities','fontsize', FONTSIZE);

% if printResults
%       print('-depsc', [dirr fileName 'TFGeneSensitivity']);
% end
%
end



for j=1:NumOfTFs
W1 = squeeze(samples.W(:,j,:))';
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
%title(j,'fontsize', FONTSIZE);
if printResults
      print('-depsc', [dirr fileName 'IntWeights' 'TF' num2str(j)]);
end
titlestring = 'Interaction weights: '; 
titlestring = [titlestring, num2str(j)];
titlestring = [titlestring, ' TF'];
title(titlestring,'fontsize', FONTSIZE);
end

W0 = samples.W0';
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
%title('W0','fontsize', FONTSIZE);
if printResults
      print('-depsc', [dirr fileName 'IntBias']);
end
title('Interaction biases','fontsize', FONTSIZE);


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
%title('W0','fontsize', FONTSIZE);
if printResults
      print('-depsc', [dirr fileName 'Delays']);
end
title('Delays','fontsize', FONTSIZE);
end
