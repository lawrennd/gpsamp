function gpmtfPlotMASAMB2(model, samples, genesAndChip, demdata)
%function gpmtfPlotMASAMB2(model, samples, genesAndChip, demdata)
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
%      -- GeneVars: if you know the Gene varariances (from PUMA)  give them here       
%
% Notes:
%       1. The confidence intervals are 95% (usually obtained with percentiles) 
%       2. In the current version the error bars of the TF profiles do not
%       contain likelihoood noise variances. If GeneVars are given then those
%       are plotted together with the observed gene expressions.    
%

FONTSIZE=8;
MARKERSIZE=6;
tfnames = {'twi', 'mef2'};

[a, b, c] = unique(genesAndChip.data, 'rows');
I = find(c <= 3);

%if strcmp(model.Likelihood.TFjointAct,'sigmoid') 
%   model.Likelihood.TFsingleAct = 'exp';
%end
%dirr = '/usr/local/michalis/mlprojects/gpsamp/tex/diagrams/';
dirr = 'figures/';

Genes = model.Likelihood.Genes;

if model.Likelihood.noiseModel.active(1) == 1
    GeneVars = model.Likelihood.noiseModel.pumaSigma2;
end
 
TimesG = model.Likelihood.TimesG + 1; 
TimesF = model.Likelihood.TimesF + 1;
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


r = 3;
J = [3, 5];

NumOfGenes = model.Likelihood.numGenes;
order = 1:NumOfGenes;
NumOfReplicas = model.Likelihood.numReplicas;
NumOfSamples = size(samples.F,2);
SizF = size(samples.F{1},2);
TimesF = TimesF(:); 
figure;
%for r=1:NumOfReplicas
  % 
  for k=1:2
    
    j = J(k);
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
     %subplot(NumOfTFs, NumOfReplicas, r + (j-1)*NumOfReplicas);
     subplot(3, 6, k+8)
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
     
     set(gca, 'FontSize', FONTSIZE);
     title(sprintf('%s TF protein', tfnames{k}),'fontsize', FONTSIZE);
     %
  end
  %
%end

%figure;
% predicted TF-Gene expressions 
if isfield(model.Likelihood,'GenesTF')
%for r=1:NumOfReplicas
  %  
  for k=1:2
    j = J(k);
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
     subplot(3, 6, k+2)
     %subplot(NumOfTFs, NumOfReplicas, r + (j-1)*NumOfReplicas);
     plot(TF,mu,'b','lineWidth',r);
     hold on;
     fillColor = [0.7 0.7 0.7];
     %fillColor = [0.8 0.8 0.8];  % for the paper
     fill([TF; TF(end:-1:1)], [mu; mu(end:-1:1)]...
            + 2*[stds; -stds(end:-1:1)], fillColor,'EdgeColor',fillColor);
     plot(TF,mu,'b','lineWidth',3);
   
     plot(TimesG,GenesTF(j,:,r),'rx','markersize', MARKERSIZE,'lineWidth', 2);
     if model.Likelihood.noiseModel.active(1) == 1
     errorbar(TimesG,  GenesTF(j,:,r), 2*sqrt(GeneTFVars(j,:,r)), 'rx','lineWidth', 1.5);
     end
     
     axis([min(TimesG(:))-0.1 max(TimesG(:))+0.1  0.95*min([GenesTF(j,:,r), mu' - stds'])  1.05*max([GenesTF(j,:,r), mu' - stds'])]);
     
     set(gca, 'FontSize', FONTSIZE);
     title(sprintf('%s TF mRNA', tfnames{k}),'fontsize', FONTSIZE);
     %
  end
  %
%end    
end

I = [2, 7, 23, 46, 53, 68];


% plot predicted gene expressions 
for j=1:length(I)
  %  
  %if mod(j,5) == 1,
  %  figure;
  %end
  %for r=1:NumOfReplicas
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
         predgen = gpmtfComputeGeneODE(LikParams, samples.F{t}(:,:,r), r, I(j));
         GG(t,:) = predgen;
         %predgen(model.Likelihood.comInds)
         %Genes(I(j),:,r)
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
     subplot(3, 6, 12+j);
     plot(TF,mu,'b','lineWidth',r);
     hold on;
     fillColor = [0.7 0.7 0.7];
     %fillColor = [0.8 0.8 0.8];  % for the paper
     fill([TF; TF(end:-1:1)], [mu; mu(end:-1:1)]...
            + 2*[stds; -stds(end:-1:1)], fillColor,'EdgeColor',fillColor);
     plot(TF,mu,'b','lineWidth',3);
   
     plot(TimesG,Genes(I(j),:,r),'rx','markersize', MARKERSIZE,'lineWidth', 2);
     if model.Likelihood.noiseModel.active(1) == 1
     errorbar(TimesG,  Genes(I(j),:,r), 2*sqrt(GeneVars(I(j),:,r)), 'rx','lineWidth', 1.5);
     end
     axis([min(TimesG(:))-0.1 max(TimesG(:))+0.1  0.95*min([Genes(I(j),:,r), mu' - stds'])  1.05*max([Genes(I(j),:,r), mu' - stds'])  ]);
     
     set(gca, 'FontSize', FONTSIZE);
     title(genesAndChip.textdata{I(j)+1, 1},'fontsize', FONTSIZE);
     %
  %end
  %
end

set(gcf, 'PaperUnits', 'centimeters')
set(gcf, 'PaperPosition', [0, 0, 18, 6])
