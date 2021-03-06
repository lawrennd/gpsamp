function gpmtfTestPlot(testGenes, Genes, GeneVars, TFs, model, fbgn, demdata, printResults, Grtruth)
% function gpmtfTestPlot(testGenes, Genes, GeneVars, TFs, model, fbgn, demdata, printResults, Grtruth)

NumPlotRows = 4;
NumPlotCols = 3;
FONTSIZE=8;

PNGSIZE=[400 300];

%% precompute the TFs
%TFset = model{end}.Likelihood.TFset; 
%numTFs = model{end}.Likelihood.numTFs; 
%load drosTrainTotal;
%for cnt=1:size(samples.F,2)
%%
%    model{end}.Likelihood.kineticsTF = samples.kineticsTF(TFset,:,cnt);
%    for r=1:model{end}.Likelihood.numReplicas
%        TFs{cnt}(:,:,r) = gpmtfComputeTFODE(model{end}.Likelihood, samples.F{cnt}(TFset,:,r), 1:numTFs);
%    end    
%    TFs{cnt} = log(TFs{cnt} + 1e-100);
%%
%end

numModels = size(testGenes,2);
NumRowsParams = 2;
NumColsParams = 4;

dirr = '~/mlprojects/gpsamp/tex/diagrams/';
dirrhtml = '~/mlprojects/gpsamp/html/';

warning off;
TimesG = model{numModels}.Likelihood.TimesG; 
TimesF = model{numModels}.Likelihood.TimesF;
%if isfield(model{numModels}.Likelihood,'GenesTF')
%    GenesTF = model{numModels}.Likelihood.GenesTF;
%    if strcmp(model{numModels}.constraints.sigmasTF,'fixed')
%        GeneTFVars = model.Likelihood.sigmasTF;
%    end
%end

        
TimesFF = TimesF(model{numModels}.Likelihood.startTime:end);
ok = date;
fileName = [demdata 'Test' 'MCMC' ok]; 

NumOfTFs = size(TFs{1},1);


NumOfGenes = model{1}.Likelihood.numGenes;
order = 1:NumOfGenes;
NumOfReplicas = model{1}.Likelihood.numReplicas;
NumOfSamples = size(testGenes{1}.LogL,2);
TimesF = TimesF(:);
SizF = size(TimesF,1);
D = 30;

%if isfield(samples, 'sigma2f') & isfield(samples, 'lengthScale')
%  % plot the variacne of the rbf kernel and its lengthscale
%  hist(samples.sigma2f, D);   
%  if printResults
%       print('-depsc', [dirr fileName 'Sigma2f']);
%       print('-dpng', [dirrhtml fileName 'Sigma2f']);
%  end
%  titlestring = 'Rbf kernel variance';
%  title(titlestring,'fontsize', 20);
%  
%  figure;
%  hist(samples.lengthScale, D);   
%  if printResults
%       print('-depsc', [dirr fileName 'LengthScale']);
%       print('-dpng', [dirrhtml fileName 'LengthScale']);
%  end
%  titlestring = 'Rbf kernel lenthscale';
%  title(titlestring,'fontsize', 20);
%  %
%end   


%meanRbf = zeros(model.Likelihood.numTimes,1);
%covRbf = zeros(model.Likelihood.numTimes, model.Likelihood.numTimes);  
% plot predicted gene expressions 



for m=1:numModels
if model{m}.Likelihood.numTFs == 0  
   GGT = zeros(NumOfSamples, model{m}.Likelihood.numTimes, NumOfReplicas);      
else
   GGT = zeros(NumOfSamples, model{m}.Likelihood.sizTime, NumOfReplicas);    
end
for t=1:NumOfSamples
  %  
  LikParams = model{m}.Likelihood; 
  if min(size(testGenes{m}.kinetics)) == 1
      LikParams.kinetics = testGenes{m}.kinetics(:)';
  else  
      LikParams.kinetics = testGenes{m}.kinetics(:,t)';
  end
  if model{m}.Likelihood.numTFs == 0
  %    
      B = LikParams.kinetics(1);
      DD = LikParams.kinetics(2);
      A = LikParams.kinetics(3);
      GGT(t,:,:) = repmat(B/DD  + (A - B/DD)*exp(-TimesG*DD)', 1, NumOfReplicas);
  else
  % 
      LikParams.W = testGenes{m}.W(:,t)';
      LikParams.W0 = testGenes{m}.W0(t);
      if isfield(testGenes{m},'TFindex')   
         LikParams.TF = TFs{testGenes{m}.TFindex(t)}(model{m}.Likelihood.TFcomb==1,:, :);
      else
         LikParams.TF = TFs{1}(model{m}.Likelihood.TFcomb==1,:, :);
      end
      
      %if model.Likelihood.noiseModel.active(3) == 1
      %mRbf = zeros(model.Likelihood.numTimes, 1);
      %sigma2f = samples.sigma2f(t);
      %lengthscale = samples.lengthScale(t);  
      %X2 = model.Likelihood.noiseModel.X2;
      %Kern = sigma2f*exp(-(0.5/(lengthscale^2)).*X2);
      %Sigma = zeros(size(Kern)); 
      %end
      for r=1:NumOfReplicas
      % 
          predgen = gpmtfComputeGeneODE(LikParams, zeros(NumOfTFs, SizF), r, 1);
          GGT(t,:,r) = predgen;
     
          % precomputation for the mean of the rbf noise model     
          %if model.Likelihood.noiseModel.active(3) == 1
          %mRbf = mRbf + (1./(model.Likelihood.noiseModel.pumaSigma2(1,:,r)')).*(Genes(1,:,r)' - predgen(model.Likelihood.comInds)');
          %Sigma = Sigma + diag(1./model.Likelihood.noiseModel.pumaSigma2(1,:,r));   
          %end
      % 
      end
      %if model.Likelihood.noiseModel.active(3) == 1
      %mmrbf(t,:) = Kern*inv(Kern + inv(Sigma))*(inv(Sigma)*mRbf);
      %meanRbf = meanRbf + mmrbf(t,:)';
      %covRbf = covRbf + Kern - Kern*inv(Kern + inv(Sigma))*Kern; 
      %end
  %    
  end
  %
end
  
  % one plot all replicas (rows) and all models (columns)
  for r=1:NumOfReplicas
      %
      subplot(NumOfReplicas, numModels, (r-1)*numModels + m);  
    
      GG = GGT(:,:,r); 
      mu = mean(GG, 1)';
      stds = sqrt(var(GG,0,1))';
      TF = TimesFF'; % TimesFF(1:2:end)';
      %figure;
      if size(mu,1) == size(TimesG(:),1)
          TF = TimesG';
      else
         TimesFF';
      end

      % max/min Gene values with offest
      maxG = max(Genes(1,:,r)) + 0.1*max(Genes(1,:,r));
      minG = min(Genes(1,:,r)) - 0.1*min(Genes(1,:,r));
      minG(minG<0)=0; 
      

      plot(TF, mu, 'b','lineWidth',r);
      hold on;
      fillColor = [0.7 0.7 0.7];
      %fillColor = [0.8 0.8 0.8];  % for the paper
      fill([TF; TF(end:-1:1)], [mu; mu(end:-1:1)]...
               + 2*[stds; -stds(end:-1:1)], fillColor,'EdgeColor',fillColor);
      plot(TF,mu,'b','lineWidth',3);


      
      plot(TimesG, Genes(1,:,r),'rx','markersize', 6,'lineWidth', 2);
      if model{m}.Likelihood.noiseModel.active(1) == 1 & sum(model{m}.Likelihood.noiseModel.active(2:3)==0)
          errorbar(TimesG,  Genes(1,:,r), 2*sqrt(GeneVars(1,:,r)), 'rx','lineWidth', 1.5);
      end
      axis([min(TimesG) max(TimesG) minG maxG]);
 
      % plot rbf GP posterior
      %if model.Likelihood.noiseModel.active(3) == 1
      %   mu = meanRbf./NumOfSamples;
      %   stds = sqrt(diag(covRbf./NumOfSamples) + mean((mmrbf - repmat(mu', NumOfSamples, 1)).^2)');
      %
      %   TG = TimesG';
      %   plot(TimesG, mu, 'g', 'lineWidth',3);
      %   hold on;
      %   fillColor = [0.8 0.8 0.8];
      %   %fillColor = [0.8 0.8 0.8];  % for the paper
      %   fill([TG; TG(end:-1:1)], [mu; mu(end:-1:1)]...
      %        + 2*[stds; -stds(end:-1:1)], fillColor,'EdgeColor',fillColor);
      %   plot(TG, mu, 'g','lineWidth',3);
      %end
      V = axis; 
      axis([min(TimesG(:))-0.1 max(TimesG(:))+0.1 V(3) V(4)]);
      %if printResults
      %   print('-depsc', [dirr fileName 'Replica' num2str(r) 'GeneExp' num2str(gg)]);
      %   print('-dpng', [dirrhtml fileName 'Replica' num2str(r) 'GeneExp' num2str(gg)]);
      %end
      titlestring = ['m:', num2str(m)];
      titlestring = [titlestring, ' r:', num2str(r)];
      titlestring = [titlestring, ' g:' fbgn];
      title(titlestring,'fontsize', FONTSIZE);
  %
  end
%
end

if printResults
    %print('-depsc', [dirr fileName 'GeneExp' fbgn]);
    set(gcf, 'PaperPosition', [0 0 PNGSIZE(1)/72 PNGSIZE(2)/72]);
    print('-dpng', '-r72', [dirrhtml fileName 'GeneExp' fbgn]);
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

modelB=zeros(numModels, 1);  modelD=zeros(numModels, 1); modelS=zeros(numModels, 1);  modelA=zeros(numModels, 1);
stdSS1=zeros(numModels, 1); stdBB1=zeros(numModels, 1); stdDD1=zeros(numModels, 1); stdAA1=zeros(numModels, 1);
stdSS2=zeros(numModels, 1); stdBB2=zeros(numModels, 1); stdDD2=zeros(numModels, 1); stdAA2=zeros(numModels, 1); 
if isfield(testGenes{1}, 'sigma2')
    modelSigma2 = zeros(numModels, 1);
    stdSigma2_1 = zeros(numModels, 1);
    stdSigma2_2 = zeros(numModels, 1);
end
modelW = zeros(numModels, model{m}.Likelihood.numTFs);
stdW_1 = zeros(numModels, model{m}.Likelihood.numTFs);
stdW_2 = zeros(numModels, model{m}.Likelihood.numTFs);
modelW0 = zeros(numModels, 1);
stdW0_1 = zeros(numModels, 1);
stdW0_2 = zeros(numModels, 1);
for m=1:numModels
% 
    if model{m}.Likelihood.numTFs > 0
        ind = find(model{m}.Likelihood.TFcomb==1);
        ok = testGenes{m}.W';
        stdok_1 = sqrt(var(ok));
        stdok_2 = sqrt(var(ok));
        modelW(m, ind) = mean(ok); 
        stdW_1(m, ind) = stdok_1;
        stdW_2(m, ind) = stdok_2;
        
        modelW0(m) = mean(testGenes{m}.W0, 2);
        stdW0_1(m) = sqrt(var(testGenes{m}.W0)); 
        stdW0_2(m) = sqrt(var(testGenes{m}.W0)); 
        
    end 
    
    if min(size(testGenes{m}.kinetics)) == 1
       BB = testGenes{m}.kinetics(1);
       DD = testGenes{m}.kinetics(2);
    else  
       BB = testGenes{m}.kinetics(1,:)';
       DD = testGenes{m}.kinetics(2,:)';
    end
    
    modelB(m) = median(BB);
    modelD(m) = median(DD);
    if model{m}.Likelihood.numTFs > 0 
       if min(size(testGenes{m}.kinetics)) == 1
          SS = testGenes{m}.kinetics(3);
          AA = testGenes{m}.kinetics(4);
       else  
          SS = testGenes{m}.kinetics(3,:)';
          AA = testGenes{m}.kinetics(4,:)';
       end
       modelS(m) = median(SS);
       stdSS1(m) = prctile(SS,5);
       stdSS2(m) = prctile(SS,95);
    else
       AA = testGenes{m}.kinetics(3,:)';  
    end
    modelA(m) = median(AA);
    stdBB1(m) = prctile(BB,5);
    stdDD1(m) = prctile(DD,5);
    stdBB2(m) = prctile(BB,95);
    stdDD2(m) = prctile(DD,95);
    stdAA1(m) = prctile(AA,5);
    stdAA2(m) = prctile(AA,95);
     
    modelSigma2(m) = median(testGenes{m}.sigma2);
    stdSigma2_1(m) = prctile(testGenes{m}.sigma2, 5);
    stdSigma2_2(m) = prctile(testGenes{m}.sigma2, 95);
end

figure;
plotIndex=1;
subplot(NumRowsParams, NumColsParams, plotIndex);
bar(modelB', 0.7); colormap([0.9 0.9 0.9]);
hold on;
errorbar([1:numModels]-0.14, modelB, modelB-stdBB1, stdBB2-modelB,'.'); 
%errorbar([1:NumOfGenes], modelB(order), modelB(order)-stdBB1(order), stdBB2(order)-modelB(order),'.'); 
title('Basal rates','fontsize', FONTSIZE);

%if printResults
%      print('-depsc', [dirr fileName 'Basal']);
%      print('-dpng', [dirrhtml fileName 'Basal']);
%end

plotIndex = plotIndex+1;
subplot(NumRowsParams, NumColsParams, plotIndex);
bar(modelD', 0.7); colormap([0.9 0.9 0.9]);
hold on;
errorbar([1:numModels]-0.14, modelD, modelD-stdDD1, stdDD2-modelD,'.'); 
title('Decay rates','fontsize', FONTSIZE);

%if printResults
%      print('-depsc', [dirr fileName 'Decay']);
%      print('-dpng', [dirrhtml fileName 'Decay']);
%end

plotIndex = plotIndex+1;
subplot(NumRowsParams, NumColsParams, plotIndex);
bar(modelS', 0.7); colormap([0.9 0.9 0.9]);
hold on;
errorbar([1:numModels]-0.14, modelS, modelS-stdSS1, stdSS2-modelS,'.'); 
title('Sensitivities','fontsize', FONTSIZE);

%if printResults
%      print('-depsc', [dirr fileName 'Sensitivity']);
%      print('-dpng', [dirrhtml fileName 'SEnsitivity']);
%end


plotIndex = plotIndex+1;
subplot(NumRowsParams, NumColsParams, plotIndex);
bar(modelA', 0.7); colormap([0.9 0.9 0.9]);
hold on;
errorbar([1:numModels]-0.14, modelA, modelA-stdAA1, stdAA2-modelA,'.'); 
title('Initial Conds','fontsize', FONTSIZE);

%if printResults
%      print('-depsc', [dirr fileName 'InitCond']);
%      print('-dpng', [dirrhtml fileName 'InitCond']);
%end


for j=1:NumOfTFs
plotIndex = plotIndex+1;
subplot(NumRowsParams, NumColsParams, plotIndex); 
bar(modelW(:, j)', 0.7); colormap([0.9 0.9 0.9]);
hold on;
errorbar([1:numModels]-0.14, modelW(:, j), stdW_1(:, j),'.'); 
%errorbar([1:NumOfGenes], modelB(order), modelB(order)-stdBB1(order), stdBB2(order)-modelB(order),'.'); 
%title(j,'fontsize', FONTSIZE);
titlestring = 'Inter. weights: '; 
titlestring = [titlestring, ' TF '];
titlestring = [titlestring, num2str(j)];
%titlestring = [titlestring, ' TF'];
title(titlestring,'fontsize', FONTSIZE);
end


plotIndex = plotIndex+1;
subplot(NumRowsParams, NumColsParams, plotIndex); 
bar(modelW0', 0.7); colormap([0.9 0.9 0.9]);
hold on;
errorbar([1:numModels]-0.14, modelW0, 2*stdW0_1,'.'); 
title('Inter. biases','fontsize', FONTSIZE);

 
% plot the variance of the added white noise
if isfield(testGenes{1}, 'sigma2')
  % plot the variance of the added white noise
  
  plotIndex = plotIndex+1;
  subplot(NumRowsParams, NumColsParams, plotIndex);
  bar(modelSigma2', 0.7); colormap([0.9 0.9 0.9]);
  hold on;
  errorbar([1:numModels]-0.14, modelSigma2, modelSigma2-stdSigma2_1, stdSigma2_2-modelSigma2,'.'); 
  title('Observ. noise','fontsize', FONTSIZE);
  
  %if nargin == 9
  %  %
  %  hold on;
  %  
  %  plot([Grtruth.sigma2 Grtruth.sigma2], [0 max(h)], 'LineWidth', 5,'Color', 'r');
  %  %
  %end
  %
end


if printResults
    %print('-depsc', [dirr fileName 'Params' fbgn]);
    set(gcf, 'PaperPosition', [0 0 PNGSIZE(1)/72 PNGSIZE(2)/72]);
    print('-dpng', '-r72', [dirrhtml fileName 'Params' fbgn]);
end
