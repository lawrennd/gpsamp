function gpmtfTestPlot16Models(testGenes, Genes, GeneVars, TFs, model, fbgn, demdata, printResults, Grtruth)
% function gpmtfTestPlot(testGenes, Genes, GeneVars, TFs, model, fbgn, demdata, printResults, Grtruth)

NumPlotRows = 4;
NumPlotCols = 3;
FONTSIZE=8;
nMod = 4; 

PNGSIZE=[400 300];


numModels = size(testGenes,2);
NumRowsParams = 3;
NumColsParams = 4;

dirr = '~/mlprojects/gpsamp/tex/diagrams/';
dirrhtml = '~/mlprojects/gpsamp/html/';

warning off;
TimesG = model{numModels}.Likelihood.TimesG; 
TimesF = model{numModels}.Likelihood.TimesF;


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

% create rhe figure identifiers  for the different replical and data fit
% plots
for r=1:NumOfReplicas
      h(r) = figure; 
end

for m=1:numModels
   if model{m}.Likelihood.numTFs == 0  
     GGT = zeros(NumOfSamples, model{m}.Likelihood.numTimes, NumOfReplicas);      
   else
     GGT = zeros(NumOfSamples, model{m}.Likelihood.sizTime, NumOfReplicas);    
   end
   for t=1:NumOfSamples
   %  
     LikParams = model{m}.Likelihood; 
     if prod(size(testGenes{m}.kinetics)) <= 4
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
        if prod(size(testGenes{m}.W)) <= NumOfTFs
           LikParams.W = testGenes{m}.W(:)';
        else
           LikParams.W = testGenes{m}.W(:,t)';
        end
        LikParams.W0 = testGenes{m}.W0(t);
        if isfield(testGenes{m},'TFindex')   
           LikParams.TF = TFs{testGenes{m}.TFindex(t)}(model{m}.Likelihood.TFcomb==1, :, :);
        else
           LikParams.TF = TFs{1}(model{m}.Likelihood.TFcomb==1, :, :);
        end
        %if model{m}.Likelihood.noiseModel.active(2) == 1
        %  LikParams.noiseModel.sigma2 = testGenes{m}.sigma2;
        %end
        %LikParams.Genes = Genes;
        %LikParams.noiseModel.pumaSigma2 = GeneVars;
        %for r=1:NumOfReplicas
        %    % call the function only with j gene expressions  
        %    newLogLik(r,:) = gpmtfLogLikelihoodGene(LikParams, model{m}.F(:,:,r), r, 1);
        %end
        for r=1:NumOfReplicas
        % 
           predgen = gpmtfComputeGeneODE(LikParams, zeros(NumOfTFs, SizF), r, 1);
           GGT(t,:,r) = predgen;
        % 
        end
     %    
     end
    %
   end
    
   % separate plots for each replica
   for r=1:NumOfReplicas
      %
      figure(h(r)); 
      hold on;
      subplot(nMod, nMod, m);  
    
      GG = GGT(:,:,r); 
      mu = mean(GG,1)';
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
      fill([TF; TF(end:-1:1)], [mu; mu(end:-1:1)]...
               + 2*[stds; -stds(end:-1:1)], fillColor,'EdgeColor',fillColor);
      plot(TF,mu,'b','lineWidth',3);
      
      plot(TimesG, Genes(1,:,r),'rx','markersize', 6,'lineWidth', 2);
      if model{m}.Likelihood.noiseModel.active(1) == 1 & sum(model{m}.Likelihood.noiseModel.active(2:3)==0)
          errorbar(TimesG,  Genes(1,:,r), 2*sqrt(GeneVars(1,:,r)), 'rx','lineWidth', 1.5);
      end
      axis([min(TimesG) max(TimesG) minG maxG]);
  
      V = axis; 
      axis([min(TimesG(:))-0.1 max(TimesG(:))+0.1 V(3) V(4)]);    
      %titlestring = [num2str(m)];
      %titlestring = [titlestring, ' r:', num2str(r)];
      %titlestring = [titlestring, ' g:' fbgn];
      title(num2str(model{m}.Likelihood.TFcomb),'fontsize', FONTSIZE);
   %
   end
%
end

if printResults
    %print('-depsc', [dirr fileName 'GeneExp' fbgn]);
    for r=1:3,
      set(h(r), 'PaperPosition', [0 0 PNGSIZE(1)/72 PNGSIZE(2)/72]);
      print(h(r), '-dpng', '-r72', [dirrhtml fileName 'GeneExp' fbgn 'Rep' num2str(r)]);
    end
end

plotIndex = r;

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
        if prod(size(testGenes{m}.W)) <= NumOfTFs
            ok = testGenes{m}.W(:)';
        else 
            ok = testGenes{m}.W';
        end
        stdok_1 = sqrt(var(ok,0,1));
        stdok_2 = sqrt(var(ok,0,1));
        
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

plotIndex = plotIndex+1;
subplot(NumRowsParams, NumColsParams, plotIndex);
bar(modelD', 0.7); colormap([0.9 0.9 0.9]);
hold on;
errorbar([1:numModels]-0.14, modelD, modelD-stdDD1, stdDD2-modelD,'.'); 
title('Decay rates','fontsize', FONTSIZE);


plotIndex = plotIndex+1;
subplot(NumRowsParams, NumColsParams, plotIndex);
bar(modelS', 0.7); colormap([0.9 0.9 0.9]);
hold on;
errorbar([1:numModels]-0.14, modelS, modelS-stdSS1, stdSS2-modelS,'.'); 
title('Sensitivities','fontsize', FONTSIZE);


plotIndex = plotIndex+1;
subplot(NumRowsParams, NumColsParams, plotIndex);
bar(modelA', 0.7); colormap([0.9 0.9 0.9]);
hold on;
errorbar([1:numModels]-0.14, modelA, modelA-stdAA1, stdAA2-modelA,'.'); 
title('Initial Conds','fontsize', FONTSIZE);




for j=1:NumOfTFs
plotIndex = plotIndex+1;
subplot(NumRowsParams, NumColsParams, plotIndex); 
bar(modelW(:, j)', 0.7); colormap([0.9 0.9 0.9]);
hold on;
errorbar([1:numModels]-0.14, modelW(:, j), 2*stdW_1(:, j),'.'); 
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
end


if printResults
    %print('-depsc', [dirr fileName 'Params' fbgn]);
    set(gcf, 'PaperPosition', [0 0 PNGSIZE(1)/72 PNGSIZE(2)/72]);
    print('-dpng', '-r72', [dirrhtml fileName 'Params' fbgn]);
end
