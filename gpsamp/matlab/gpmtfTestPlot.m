function gpmtfTestPlot(model, Genes, GeneVars, samples, demdata, printResults)

dirr = '/usr/local/michalis/mlprojects/gpsamp/tex/diagrams/';

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
fileName = [demdata 'Test' 'MCMC' ok model.Likelihood.singleAct model.Likelihood.jointAct]; 

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
      GG(t,:) = samples.predGenes(r,:,t);
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
hist(samples.kinetics(1,:));
title('Basal rates','fontsize', 20);

if printResults
      print('-depsc', [dirr fileName 'Basal']);
end
  figure;   
hist(samples.kinetics(3,:));
title('Sensitivities','fontsize', 20);

if printResults
      print('-depsc', [dirr fileName 'Sensitivity']);
end
figure;
hist(samples.kinetics(2,:));
title('Decays','fontsize', 20);

if printResults
      print('-depsc', [dirr fileName 'Decay']);
end
figure;
hist(samples.kinetics(4,:));
title('Initial conditions','fontsize', 20);

if printResults
      print('-depsc', [dirr fileName 'InitCond']);
end


for j=1:NumOfTFs
W1 = samples.Weights(j,:);
figure;
hist(W1);
if printResults
     print('-depsc', [dirr fileName 'IntWeights' 'TF' num2str(j)]);
end
titlestring = 'Interaction weights: '; 
titlestring = [titlestring, num2str(j)];
titlestring = [titlestring, ' TF'];
title(titlestring,'fontsize', 20);
end

W0 = samples.Weights0';
figure;
hist(W0);
if printResults
      print('-depsc', [dirr fileName 'IntBias']);
end
title('Interaction biases','fontsize', 20);
