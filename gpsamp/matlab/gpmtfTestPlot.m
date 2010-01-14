function gpmtfTestPlot(model, Genes, GeneVars, fbgn, samples, TFs, demdata, printResults, Grtruth)

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
      LikParams.W = samples.W(:,t)';
      LikParams.W0 = samples.W0(t);
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
  if model.Likelihood.noiseModel.active(1) == 1 & sum(model.Likelihood.noiseModel.active(2:3)==0)
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
W1 = samples.W(j,:);
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

W0 = samples.W0';

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
