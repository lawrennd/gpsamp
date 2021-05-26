function gpmtfTestPlotSeparatemRNA(testGenes, Genes, GeneVars, TFs, model, fbgn, demdata, printResults, dirr)


% USER DEFINED PARAMETERS 
FONTSIZE=8;
PNGSIZE=[400 300];
timeshift = 1;
XLabel = 'time (h)';
NumRowsParams = 2;
NumColsParams = 4;
% END OF USER DEFINED PARAMETERS

dirr = '~/mlprojects/gpsamp/tex/diagrams/';
dirrhtml = '~/mlprojects/gpsamp/html/';
% dirrhtml =  dirr; 

warning off;
numModels = size(testGenes,2);
TimesG = model{numModels}.Likelihood.TimesG; 
TimesF = model{numModels}.Likelihood.TimesF;

        
TimesFF = TimesF(model{numModels}.Likelihood.startTime:end);
ok = datestr(now, 29); % '2010-11-10';
fileName = [demdata 'Test' 'MCMC' ok]; 

NumOfTFs = size(TFs{1},1);


NumOfGenes = model{1}.Likelihood.numGenes;
order = 1:NumOfGenes;
NumOfReplicas = model{1}.Likelihood.numReplicas;
NumOfSamples = size(testGenes{1}.LogL,2);
TimesF = TimesF(:);
SizF = size(TimesF,1);
D = 30;


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
  
  % one plot all replicas (rows) and all models (columns)
  for r=1:NumOfReplicas
      %
      figure; 
   
      GG = GGT(:,:,r); 
      mu = mean(GG, 1)';
      stds = sqrt(var(GG,0,1))';
      TF = TimesFF' + timeshift; 
      %figure;
      if size(mu,1) == size(TimesG(:),1)
          TF = TimesG' + timeshift;
      else
         TimesFF';
      end

      % max/min Gene values with offest
      maxG = max(Genes(1,:,r) + 2*sqrt(GeneVars(1,:,r))) + 0.1*max(Genes(1,:,r));
      minG = 0; %min(Genes(1,:,r)) - 0.1*min(Genes(1,:,r)); 
      

      plot(TF, mu, 'b','lineWidth',r);
      hold on;
      fillColor = [0.7 0.7 0.7];
      %fillColor = [0.8 0.8 0.8];  % for the paper
      fill([TF; TF(end:-1:1)], [mu; mu(end:-1:1)]...
               + 2*[stds; -stds(end:-1:1)], fillColor,'EdgeColor',fillColor);
      plot(TF,mu,'b','lineWidth',3);

      plot(TimesG+timeshift, Genes(1,:,r),'rx','markersize', 6,'lineWidth', 2);
      if model{m}.Likelihood.noiseModel.active(1) == 1 & sum(model{m}.Likelihood.noiseModel.active(2:3)==0)
          errorbar(TimesG+timeshift,  Genes(1,:,r), 2*sqrt(GeneVars(1,:,r)), 'rx','lineWidth', 1.5);
      end
      axis([min(TimesG)+timeshift max(TimesG)+timeshift minG maxG]);
 
      V = axis; 
      axis([min(TimesG(:))-0.1+timeshift max(TimesG(:))+0.1+timeshift V(3) V(4)]);

      set(gca, 'FontSize', FONTSIZE)
      set(gca, 'YTickLabel', [])
      xlabel('time (h)')
      set(gcf, 'PaperUnits', 'centimeters')
      set(gcf, 'PaperPosition', [0 0 4.5 3.5])

      if printResults
         print('-depsc2',[dirr fileName 'Model'  num2str(model{m}.Likelihood.TFcomb,'%1g') 'Replica' num2str(r) 'GeneExp' fbgn]);
         print('-dpng', [dirrhtml fileName 'Model'  num2str(model{m}.Likelihood.TFcomb,'%1g')  'Replica' num2str(r) 'GeneExp' fbgn]);
      end
      %titlestring = ['m:', num2str(m)];
      %titlestring = [titlestring, ' r:', num2str(r)];
      %titlestring = [titlestring, ' g:' fbgn];
      %title(titlestring,'fontsize', FONTSIZE);
  %
  end
%
end

