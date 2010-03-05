function gpmtfTestPlot(model, Genes, GeneVars, fbgn, samples, TFs, demdata, printResults, Grtruth)

NumPlotRows = 4;
NumPlotCols = 3;
FONTSIZE=8;

dirr = '/usr/local/michalis/mlprojects/gpsamp/tex/diagrams/';
dirrhtml = '/usr/local/michalis/mlprojects/gpsamp/html/';

%Genes = model.Likelihood.Genes;
gg = 1;
warning off;

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
D = 30;

if isfield(samples, 'sigma2f') & isfield(samples, 'lengthScale')
  % plot the variacne of the rbf kernel and its lengthscale
  hist(samples.sigma2f, D);   
  if printResults
       print('-depsc', [dirr fileName 'Sigma2f']);
       print('-dpng', [dirrhtml fileName 'Sigma2f']);
  end
  titlestring = 'Rbf kernel variance';
  title(titlestring,'fontsize', 20);
  
  figure;
  hist(samples.lengthScale, D);   
  if printResults
       print('-depsc', [dirr fileName 'LengthScale']);
       print('-dpng', [dirrhtml fileName 'LengthScale']);
  end
  titlestring = 'Rbf kernel lenthscale';
  title(titlestring,'fontsize', 20);
  %
end   

if isfield(samples, 'sigma2')
  % plot the variance of the added white noise
  figure;
  h = hist(samples.sigma2, D);
  hist(samples.sigma2, D);   
  if printResults
       print('-depsc', [dirr fileName 'Sigma2']);
       print('-dpng', [dirrhtml fileName 'Sigma2']);
  end
  titlestring = 'White noise variance';
  title(titlestring,'fontsize', 20);
  
  if nargin == 9
    %
    hold on;
    
    plot([Grtruth.sigma2 Grtruth.sigma2], [0 max(h)], 'LineWidth', 5,'Color', 'r');
    %
  end
  %
end   


meanRbf = zeros(model.Likelihood.numTimes,1);
covRbf = zeros(model.Likelihood.numTimes, model.Likelihood.numTimes);  
GGT = zeros(NumOfSamples, model.Likelihood.sizTime, NumOfReplicas);   
% plot predicted gene expressions 
for t=1:NumOfSamples
  % 
  LikParams = model.Likelihood; 
  LikParams.kinetics = samples.kinetics(:,t)';
  LikParams.W = samples.W(:,t)';
  LikParams.W0 = samples.W0(t);
  LikParams.TF = TFs{samples.TFindex(t)};
 
  if model.Likelihood.noiseModel.active(3) == 1
  mRbf = zeros(model.Likelihood.numTimes, 1);
  sigma2f = samples.sigma2f(t);
  lengthscale = samples.lengthScale(t);  
  X2 = model.Likelihood.noiseModel.X2;
  Kern = sigma2f*exp(-(0.5/(lengthscale^2)).*X2);
  Sigma = zeros(size(Kern)); 
  end
  
  for r=1:NumOfReplicas
     % 
     predgen = gpmtfComputeGeneODE(LikParams, zeros(NumOfTFs, SizF), r, 1);
     GGT(t,:,r) = predgen;
     
     % precomputation for the mean of the rbf noise model
     
    if model.Likelihood.noiseModel.active(3) == 1
     mRbf = mRbf + (1./(model.Likelihood.noiseModel.pumaSigma2(1,:,r)')).*(Genes(1,:,r)' - predgen(model.Likelihood.comInds)');

     Sigma = Sigma + diag(1./model.Likelihood.noiseModel.pumaSigma2(1,:,r));   
    end
     % 
  end
  if model.Likelihood.noiseModel.active(3) == 1
  mmrbf(t,:) = Kern*inv(Kern + inv(Sigma))*(inv(Sigma)*mRbf);
  meanRbf = meanRbf + mmrbf(t,:)';
  covRbf = covRbf + Kern - Kern*inv(Kern + inv(Sigma))*Kern; 
  end
  %
end


for r=1:NumOfReplicas
  %  
  GG = GGT(:,:,r); 
  
  mu = mean(GG)';
  stds = sqrt(var(GG))';
    
  TF = TimesFF'; % TimesFF(1:2:end)';
  figure;
  plot(TF, mu, 'b','lineWidth',r);
  hold on;
  fillColor = [0.7 0.7 0.7];
  %fillColor = [0.8 0.8 0.8];  % for the paper
  fill([TF; TF(end:-1:1)], [mu; mu(end:-1:1)]...
        + 2*[stds; -stds(end:-1:1)], fillColor,'EdgeColor',fillColor);
  plot(TF,mu,'b','lineWidth',3);
   
  plot(TimesG,Genes(1,:,r),'rx','markersize', 6,'lineWidth', 2);
  if model.Likelihood.noiseModel.active(1) == 1 & sum(model.Likelihood.noiseModel.active(2:3)==0)
     errorbar(TimesG,  Genes(1,:,r), 2*sqrt(GeneVars(1,:,r)), 'rx','lineWidth', 1.5);
  end
 
  % plot rbf GP posterior
  if model.Likelihood.noiseModel.active(3) == 1
  mu = meanRbf./NumOfSamples;
  stds = sqrt(diag(covRbf./NumOfSamples) + mean((mmrbf - repmat(mu', NumOfSamples, 1)).^2)');
 
  TG = TimesG';
  plot(TimesG, mu, 'g', 'lineWidth',3);
  hold on;
  fillColor = [0.8 0.8 0.8];
  %fillColor = [0.8 0.8 0.8];  % for the paper
  fill([TG; TG(end:-1:1)], [mu; mu(end:-1:1)]...
        + 2*[stds; -stds(end:-1:1)], fillColor,'EdgeColor',fillColor);
  plot(TG, mu, 'g','lineWidth',3);
  end
  
  V = axis; 
  axis([min(TimesG(:))-0.1 max(TimesG(:))+0.1 V(3) V(4)]);
     
  if printResults
     print('-depsc', [dirr fileName 'Replica' num2str(r) 'GeneExp' num2str(gg)]);
     print('-dpng', [dirrhtml fileName 'Replica' num2str(r) 'GeneExp' num2str(gg)]);
  end
  titlestring = 'Expressions: ';
  titlestring = [titlestring, num2str(r)]; 
  titlestring = [titlestring, ' replica, '];
  titlestring = [titlestring, num2str(gg)];
  titlestring = [titlestring, ' gene'];
  title(titlestring,'fontsize', FONTSIZE);
  %
  %
end
plotIndex = r;


%figure;
%stds = sqrt(diag(covRbf./NumOfSamples));
%mu = meanRbf./NumOfSamples;
%
%TG = TimesG';
%plot(TimesG, mu, 'b', 'lineWidth',2);
%hold on;
%fillColor = [0.7 0.7 0.7];
%%fillColor = [0.8 0.8 0.8];  % for the paper
%fill([TG; TG(end:-1:1)], [mu; mu(end:-1:1)]...
%        + 2*[stds; -stds(end:-1:1)], fillColor,'EdgeColor',fillColor);
%plot(TG, mu, 'b','lineWidth',3);
%
%for r=1:NumOfReplicas
%   plot(TimesG, Genes(1,:,r),'rx','markersize', 14','lineWidth', 2);
%   if model.Likelihood.noiseModel.active(1) == 1 & sum(model.Likelihood.noiseModel.active(2:3)==0)
%      errorbar(TimesG,  Genes(1,:,r), 2*sqrt(GeneVars(1,:,r)), 'rx','lineWidth', 1.5);
%   end
%   %axis([min(TimesG(:))-0.1 max(TimesG(:))+0.1 0.95*min(min(Genes(1,:,r))) 1.05*max(max(Genes(1,:,r)))]);
%end
%if printResults
%   print('-depsc', [dirr fileName 'Replica' num2str(r) 'GeneExp' num2str(gg)]);
%   print('-dpng', [dirrhtml fileName 'Replica' num2str(r) 'GeneExp' num2str(gg)]);
%end
%titlestring = 'Rbf GP posterior';
%title(titlestring,'fontsize', 20);


figure;
h = hist(samples.kinetics(1,:),D);
hist(samples.kinetics(1,:),D);
title('Basal rates','fontsize', FONTSIZE);
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
plotIndex = plotIndex+1;
subplot(NumPlotRows, NumPlotCols, plotIndex);
h = hist(samples.kinetics(3,:), D);
hist(samples.kinetics(3,:), D);
title('Sensitivities','fontsize', FONTSIZE);
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
plotIndex = plotIndex+1;
subplot(NumPlotRows, NumPlotCols, plotIndex);
h = hist(samples.kinetics(2,:), D);
hist(samples.kinetics(2,:), D);
title('Decays','fontsize', FONTSIZE);
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
plotIndex = plotIndex+1;
subplot(NumPlotRows, NumPlotCols, plotIndex);
h = hist(samples.kinetics(4,:), D);
hist(samples.kinetics(4,:), D);
title('Initial conditions','fontsize', FONTSIZE);
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
W1 = samples.W(j,:);
plotIndex = plotIndex+1;
subplot(NumPlotRows, NumPlotCols, plotIndex);
h = hist(W1, D);
hist(W1, D);
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
title(titlestring,'fontsize', FONTSIZE);
end

W0 = samples.W0';

plotIndex = plotIndex+1;
subplot(NumPlotRows, NumPlotCols, plotIndex);
h = hist(W0, D);
hist(W0, D);
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
title('Interaction biases','fontsize', FONTSIZE);
