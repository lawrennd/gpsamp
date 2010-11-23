function gpmtfPlot(model, samples, demdata, TFname, Genename, printResults, dirr)
%function gpmtfPlot(model, samples, demdata, TFname, printResults, ddir)
%
% Description: Creates plots to display the outcome of MCMC in the traning phase 
%
% Inputs: 
%      -- model: Contains the structure of the GP model 
%      -- samples: A structure that contains the samples
%      -- demdata: a string that characterizes the experiments, e.g. 'p53' 
%      -- printResults: if 0 then plots will not be printed to files 
%                       if 1 the plots will be printed to the directory
%                       dirr
%      -- dirr: directory where the plots  are going to be printed 
%
% Things to check:
%       1. The confidence intervals are 95% (usually obtained with percentiles) 
%       2. In the current version the error bars of the TF profiles do not
%       contain likelihoood noise variances. If GeneVars are given then those
%       are plotted together with the observed gene expressions.    
%

% USER DEFINED PARAMETERS
ok = '15-Nov-2010';% date; % '24-Oct-2010';  %date;
fileName = [demdata 'MCMC' ok model.Likelihood.singleAct model.Likelihood.jointAct]; 
FONTSIZE=10;
FS2 = 16;
timeshift = 1;
XLabel = 'time (h)';
% how gene per plot
ManyGenes = 5;
% separate plots or grouped 
SepPlot = 0;
% allow different colours fot the TF profiles 
%colSet = {'b','r','g'}; 
colSet = {'b','b','b', 'b','b'}; 
%scTFAxis = [0 1.2];
scTFAxis = [];
% write replicas or conditions 
rep = 'rep';
%rep = 'cond'; 
% write the wrod rep or cond or not
Wrtrep = 1;
% which replicas are plotted in for the mRNA data  (default is all)
sRep = ones(1,model.Likelihood.numReplicas); 
%sRep = [1 1 0];

%demdata = 'demEcoli';
%order = [1 5 3 4 2];
%order = 1:NumOfGenes;
%order = [1 5 3 4 2]; % for the Barenco data 
%order = 1:1:size(Y,1);  % for ecoli or other data
%DataID = 1; % 1 for Barenco data, 2 for Ecoli data
%bbar = 1;  % draw or not erroor when ploting basal , decay and sensitivities
%dataIDvasknownPuma = 1;  % plot the errors bard of the observed gene expressions if they have been
                         % computed by some package e.g. PUMA
%%%%%%  END OF USER DEFINED PARAMETERS


if nargin < 6
 dirr = '/usr/local/michalis/mlprojects/gpsamp/tex/diagrams/';
end

dirrhtml =  dirr;  

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


NumOfTFs = model.Likelihood.numTFs;
NumOfGenes = model.Likelihood.numGenes;
order = 1:NumOfGenes;
NumOfReplicas = model.Likelihood.numReplicas;
NumOfSamples = size(samples.LogL,2);
if isfield(samples, 'F')
    SizF = size(samples.F{1},2);
else
    SizF = size(TimesF,2);
end
TimesF = TimesF(:); 


% PLOT PREDICTED TFS
%
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
         if isfield(samples, 'F')
           PredTF(t,:) = gpmtfComputeTF(lik, samples.F{t}(j,:,r), j);
         else
           % the baseline model is used  
           PredTF(t,:) = gpmtfComputeTF(lik, model.F(j,:,r), j);  
         end
     end
     mu = median(PredTF,1)';
     stds1 = (prctile(PredTF,95,1)'-mu)/2;
     stds2 = (mu-prctile(PredTF,5,1)')/2;
     
    if SepPlot==0      
         subplot(NumOfReplicas, NumOfTFs, (r-1)*NumOfTFs + j);
     elseif (r*j) > 1
         figure;
     end
     TimesFP = TimesF + timeshift;
     plot(TimesFP,mu, colSet{j}, 'lineWidth',3);
     hold on;
     fillColor = [0.7 0.7 0.7];
     %fillColor = [0.8 0.8 0.8];  % for the paper
     fill([TimesFP; TimesFP(end:-1:1)], [mu; mu(end:-1:1)]...
            + 2*[stds1; -stds2(end:-1:1)], fillColor,'EdgeColor',fillColor);
     plot(TimesFP,mu,colSet{j},'lineWidth',3);
     
     %if isempty(scTFAxis) 
     %   axis([TimesF(1) TimesF(end)+0.1 0 max(mu+2*stds1)+0.1]);
     %else
     %   axis([TimesF(1) TimesF(end)+0.1 scTFAxis(1) scTFAxis(2)]);
     %end
     
     % plot the ground truth if exist
     if isfield(model,'groundtr') == 1
        FFgt = model.groundtr.TF(j,:,r);
        %FFgt = feval(model.Likelihood.TFsingleAct,model.GroundTruth.F(j,:,r));
        plot(TimesFP,FFgt,'r','lineWidth',3);
     end
     
     axis tight;
     set(gca, 'FontSize', FONTSIZE);
     set(gca, 'YTickLabel', []);
     xlabel(XLabel);
     set(gcf, 'PaperUnits', 'centimeters');
     set(gcf, 'PaperPosition', [0 0 4.5 3.5]);
     if printResults & SepPlot==1
         print('-depsc', [dirr fileName 'TFproteins_Rep' num2str(r) 'TF' num2str(j)]); 
         print('-dpng', [dirrhtml fileName 'TFproteins_Rep' num2str(r) 'TF' num2str(j)]);
     end
     %titlestring = 'TF protein: '; 
     titlestring = [];
     titlestring = [titlestring, TFname{j}];
     % if NumOfReplicas > 1
     %    if Wrtrep == 1 
     %    titlestring = [titlestring, ', ', rep, ' ', num2str(r)];   
     %    end    
     % end
     %titlestring = [titlestring, ' TF'];
     title(titlestring,'fontsize', FONTSIZE);
     %
  end
  %
end
 
if printResults & SepPlot==0
  print('-depsc', [dirr fileName 'TFproteins']);
end 

if SepPlot==0 
   figure;
end
%
% END OF PLOTING PREDTICTED TFs


% PLOT PREDICTED TF mRNA EXPRESSION PROFILES  
%
if isfield(model.Likelihood,'GenesTF')
for r=1:NumOfReplicas
  %  
  for j=1:NumOfTFs
     % 
     GG = zeros(NumOfSamples, size(model.Likelihood.TimesF,2));    
     for t=1:NumOfSamples
         if isfield(samples, 'F')
            PredGenesTF = singleactFunc(model.Likelihood.singleAct, samples.F{t}(j,:,r));
         else
            PredGenesTF = model.F(j,:,r);
         end
         GG(t,:) = PredGenesTF;
     end
     
     mu = mean(GG,1)';
     stds = sqrt(var(GG,0,1))';
   
     if SepPlot==0 
       subplot(NumOfReplicas, NumOfTFs, (r-1)*NumOfTFs + j);
     else
       figure;
     end
    
     TimesFP = TimesF + timeshift; 
     TimesGP = TimesG + timeshift; 
     plot(TimesFP,mu,'b','lineWidth',r);
     hold on;
     fillColor = [0.7 0.7 0.7];
     %fillColor = [0.8 0.8 0.8];  % for the paper
     fill([TimesFP; TimesFP(end:-1:1)], [mu; mu(end:-1:1)]...
            + 2*[stds; -stds(end:-1:1)], fillColor,'EdgeColor',fillColor);
     plot(TimesFP,mu,'b','lineWidth',3);
   
     plot(TimesGP,GenesTF(j,:,r),'rx','markersize', 14','lineWidth', 2);
     if model.Likelihood.noiseModel.active(1) == 1
     errorbar(TimesGP,  GenesTF(j,:,r), 2*sqrt(GeneTFVars(j,:,r)), 'rx','lineWidth', 1.5);
     end
     
     axis tight;
     set(gca, 'FontSize', FONTSIZE);
     %set(gca, 'YTickLabel', []);
     xlabel(XLabel);
     set(gcf, 'PaperUnits', 'centimeters');
     set(gcf, 'PaperPosition', [0 0 4.5 3.5]);
     %axis([min(TimesG(:))-0.1 max(TimesG(:))+0.1  0.95*min([GenesTF(j,:,r), mu' - stds'])  1.05*max([GenesTF(j,:,r), mu' - stds'])]);
     
     if printResults & SepPlot==1
         print('-depsc', [dirr fileName 'TFmRNAs_Rep' num2str(r) 'Gene' num2str(j)]);
     end
     
     %titlestring = 'TF mRNA: ';
     titlestring = [];
     titlestring = [titlestring, TFname{j}];
     if NumOfReplicas > 1
        if Wrtrep == 1 
        titlestring = [titlestring, ', ' rep, ' ', num2str(r)];
        end
     end
     title(titlestring,'fontsize', FONTSIZE);
     %
  end
  %
end    

if printResults & SepPlot==0
  print('-depsc', [dirr fileName  'TFmRNAs']);
end
end
%
% END OF PLOTING THE PREDICTED TF mRNAs


% PLOT PREDICTED mRNA PROFILES FOR THE TARGET GENES  
%
for j=1:NumOfGenes
  % 
  if SepPlot==0 
    if mod(j,ManyGenes) == 1,
      figure;
    end
  end
  
  rrc = 0;
  for r=1:NumOfReplicas 
     if sRep(r) == 1
     rrc = rrc + 1;    
     % 
     GG = zeros(NumOfSamples,model.Likelihood.sizTime);    
     for t=1:NumOfSamples
         %
         LikParams = model.Likelihood; 
         LikParams.kinetics = samples.kinetics(:,:,t);
         LikParams.kineticsTF = samples.kineticsTF(:,:,t);
         LikParams.W = samples.W(:,:,t);
         LikParams.W0 = samples.W0(:,t);
         if isfield(samples, 'F')
            predgen = gpmtfComputeGeneODE(LikParams, samples.F{t}(:,:,r), r, j);
         else
            predgen = gpmtfComputeGeneODE(LikParams, model.F(:,:,r), r, j);
         end
         GG(t,:) = predgen;
         %
     end
     
     mu = mean(GG,1)';
     stds = sqrt(var(GG,0,1))';
    
     TF = TimesFF'; 
     TFP = TF + timeshift; 
     TimesGP = TimesG + timeshift; 
     if SepPlot==0 
        subplot(size(sRep,2), ManyGenes, (rrc-1)*ManyGenes + 1 + mod(j-1, ManyGenes));
     else
        figure;
     end
     
     %subplot(5, NumOfReplicas, r + mod(j-1, 5)*NumOfReplicas);
     plot(TFP,mu,'b','lineWidth',r);
     hold on;
     fillColor = [0.7 0.7 0.7];
     %fillColor = [0.8 0.8 0.8];  % for the paper
     fill([TFP; TFP(end:-1:1)], [mu; mu(end:-1:1)]...
            + 2*[stds; -stds(end:-1:1)], fillColor,'EdgeColor',fillColor);
     plot(TFP,mu,'b','lineWidth',3);
   
     plot(TimesGP,Genes(j,:,r),'rx','markersize', 14','lineWidth', 2);
     if model.Likelihood.noiseModel.active(1) == 1
        errorbar(TimesGP,  Genes(j,:,r), 2*sqrt(GeneVars(j,:,r)), 'rx','lineWidth', 1.5);
     end
     
     axis tight;
     set(gca, 'FontSize', FONTSIZE);
     %set(gca, 'YTickLabel', []);
     xlabel(XLabel);
     set(gcf, 'PaperUnits', 'centimeters');
     set(gcf, 'PaperPosition', [0 0 4.5 3.5]);
     %axis([min(TimesGP(:))-0.1 max(TimesG(:))+0.1  0.95*min([Genes(j,:,r), mu' - stds'])  1.05*max([Genes(j,:,r), mu' - stds'])  ]);
     
     
     %titlestring = 'mRNA: '; 
     titlestring = [];
     if isnumeric(Genename(j))
        titlestring = [titlestring, num2str(j)];
     else
        titlestring = [titlestring, Genename{j}]; 
     end
     if NumOfReplicas > 1
        if Wrtrep == 1 
        titlestring = [titlestring, ', ' rep, ' ', num2str(r)];
        end
     end
     title(titlestring,'fontsize', FONTSIZE);
          
     if SepPlot==0 
        if printResults && r==sum(sRep) && (mod(j,ManyGenes)==0 || j==NumOfGenes),
        if NumOfReplicas > 1
           print('-depsc', [dirr fileName 'Replica' num2str(sRep,'%1g') 'GENEmRNA' num2str(j)]);
        else
           print('-depsc', [dirr fileName 'GENEmRNA' num2str(j)]); 
        end
        end
     else
        print('-depsc', [dirr fileName 'Replica' num2str(sRep,'%1g') 'GENEmRNA' num2str(j)]);
     end
     end
     %
  end
  %
end
%
% END PLOTING PREDICTED mRNA PROFILES FOR THE TARGET GENES


ok = mean(samples.kinetics,3);  
BB = squeeze(samples.kinetics(:,1,:))';
DD = squeeze(samples.kinetics(:,2,:))';
SS = squeeze(samples.kinetics(:,3,:))';
AA = squeeze(samples.kinetics(:,4,:))';
modelB = median(BB,1);
modelS = median(SS,1);
modelD = median(DD,1);
modelA = median(AA,1);

stdBB1 = prctile(BB,5,1);
stdDD1 = prctile(DD,5,1);
stdSS1 = prctile(SS,5,1);
stdBB2 = prctile(BB,95,1);
stdDD2 = prctile(DD,95,1);
stdSS2 = prctile(SS,95,1);
stdAA1 = prctile(AA,5,1);
stdAA2 = prctile(AA,95,1);

% Plot first basal transcription rates.
figure;
if isfield(model,'groundtr') == 1
bar([modelB(order); model.groundtr.kinetics(:,1)']', 0.7); colormap([0.9 0.9 0.9; 0 0 0]);
else
bar(modelB(order)', 0.7); colormap([0.9 0.9 0.9]);
end
hold on;
errorbar([1:NumOfGenes]-0.14, modelB(order), modelB(order)-stdBB1(order), stdBB2(order)-modelB(order),'.'); 
axis tight;
set(gca,'fontsize',FS2);
%errorbar([1:NumOfGenes], modelB(order), modelB(order)-stdBB1(order), stdBB2(order)-modelB(order),'.'); 
title('Basal rates (b_j)','fontsize', FS2);

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
axis tight;
set(gca,'fontsize',FS2);
%errorbar([1:NumOfGenes], modelS(order), modelS(order)-stdSS1(order), stdSS2(order)-modelS(order),'.');
title('Sensit. (s_j)','fontsize', FS2);

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
axis tight;
set(gca,'fontsize',FS2);
%errorbar([1:NumOfGenes], modelD(order), modelD(order)-stdDD1(order), stdDD2(order)-modelD(order),'.');
title('mRNA degr. (d_j)','fontsize', FS2);

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
axis tight;
set(gca,'fontsize',FS2);
%errorbar([1:NumOfGenes], modelA(order), modelA(order)-stdAA1(order), stdAA2(order)-modelA(order),'.');
title('Init. conds (a_j)','fontsize', FS2);
if printResults
      print('-depsc', [dirr fileName 'InitCond']);
end


if isfield(model.Likelihood,'GenesTF')
%
orderTF = 1:NumOfTFs;
ok = mean(samples.kineticsTF,3);  
DD = squeeze(samples.kineticsTF(:,1,:))';
SS = squeeze(samples.kineticsTF(:,2,:))';
AA = squeeze(samples.kineticsTF(:,3,:))';
modelS = median(SS,1);
modelD = median(DD,1);
modelA = median(AA,1);

stdDD1 = prctile(DD,5,1);
stdSS1 = prctile(SS,5,1);
stdDD2 = prctile(DD,95,1);
stdSS2 = prctile(SS,95,1);
stdAA1 = prctile(AA,5,1);
stdAA2 = prctile(AA,95,1);

figure;
% plot degradation rates
if isfield(model,'groundtr') == 1
bar([modelD(orderTF); model.groundtr.kineticsTF(:,1)']', 0.7); colormap([0.9 0.9 0.9; 0 0 0]);
else
bar(modelD(orderTF)', 0.7); colormap([0.9 0.9 0.9]);
end
hold on;
errorbar([1:NumOfTFs]-0.14, modelD(orderTF), modelD(orderTF)-stdDD1(orderTF), stdDD2(orderTF)-modelD(orderTF),'.');
set(gca, 'xticklabel', TFname);
axis tight;


set(gca,'fontsize',FS2);
title('TF mRNA degr. (\delta_i)','fontsize', FS2);
if printResults
      print('-depsc', [dirr fileName 'TFGeneDecay']);
end

figure;
% plot degradation rates
if isfield(model,'groundtr') == 1
bar([modelA(orderTF); model.groundtr.kineticsTF(:,3)']', 0.7); colormap([0.9 0.9 0.9; 0 0 0]);
else
bar(modelA(orderTF)', 0.7); colormap([0.9 0.9 0.9]);
end
hold on;
errorbar([1:NumOfTFs]-0.14, modelA(orderTF), modelA(orderTF)-stdAA1(orderTF), stdAA2(orderTF)-modelA(orderTF),'.');
set(gca, 'xticklabel', TFname);
axis tight;
set(gca,'fontsize',FS2);
title('TF mRNA Init. Cond. (c_i)','fontsize', FS2);
if printResults
      print('-depsc', [dirr fileName 'TFGeneInitCond']);
end
%
end


figure;
% plot lengthscales
modelLengsc = mean(samples.lengthScale');
stdLengsc_1 = prctile(samples.lengthScale',5,1);
stdLengsc_2 = prctile(samples.lengthScale',95,1);
%% plot initial conditions
%if isfield(model,'groundtr') == 1
%bar([modelTaus(order); model.groundtr.Taus]', 0.7); colormap([0.9 0.9 0.9; 0 0 0]);
%else
bar(modelLengsc', 0.7); colormap([0.9 0.9 0.9]);
%end
hold on;
errorbar([1:NumOfTFs]-0.14, modelLengsc, modelLengsc-stdLengsc_1, stdLengsc_2-modelLengsc,'.');
set(gca,'fontsize',FS2);
axis tight;
%errorbar([1:NumOfGenes], modelA(order), modelA(order)-stdAA1(order), stdAA2(order)-modelA(order),'.');
%title('W0','fontsize', FONTSIZE);
if printResults
      print('-depsc', [dirr fileName 'Lengthscales']);
end
title('Lengthscales','fontsize', FS2);




for j=1:NumOfTFs
W1 = squeeze(samples.W(:,j,:))';
modelW1 = mean(W1,1);
stdW1_1 = sqrt(var(W1,0,1));
stdW1_2 = sqrt(var(W1,0,1));
% Plot first basal transcription rates.
figure;
if isfield(model,'groundtr') == 1
bar([modelW1(order); model.groundtr.W(:,j)']', 0.7); colormap([0.9 0.9 0.9; 0 0 0]);
else 
bar(modelW1(order)', 0.7); colormap([0.9 0.9 0.9]);
end
hold on;
errorbar([1:NumOfGenes]-0.14, modelW1(order), 2*stdW1_1(order),'.'); 
axis tight;
set(gca,'fontsize',FS2);

%titlestring = 'Interaction weights: '; 
titlestring = []; 
titlestring = [titlestring, TFname{j}];
%titlestring = [titlestring, ' TF'];
title(titlestring,'fontsize', FS2);
if printResults
      print('-depsc', [dirr fileName 'IntWeights' 'TF' num2str(j)]);
end

end

W0 = samples.W0';
modelW0 = mean(W0,1);
stdW0_1 = sqrt(var(W0,0,1));
stdW0_2 = sqrt(var(W0,0,1));
figure;
% plot initial conditions
if isfield(model,'groundtr') == 1
bar([modelW0(order); model.groundtr.W0']', 0.7); colormap([0.9 0.9 0.9; 0 0 0]);
else
bar(modelW0(order)', 0.7); colormap([0.9 0.9 0.9]);
end
hold on;
errorbar([1:NumOfGenes]-0.14, modelW0(order), 2*stdW0_1(order),'.');
axis tight;
set(gca,'fontsize',FS2);
title('Biases','fontsize', FS2);
if printResults
      print('-depsc', [dirr fileName 'IntBias']);
end


if model.Likelihood.noiseModel.active(2) == 1
figure;
% plot Sigma2s
modelSigma2 = mean(samples.sigma2',1);
stdSigma2_1 = prctile(samples.sigma2',5,1);
stdSigma2_2 = prctile(samples.sigma2',95,1);
%% plot initial conditions
%if isfield(model,'groundtr') == 1
%bar([modelTaus(order); model.groundtr.Taus]', 0.7); colormap([0.9 0.9 0.9; 0 0 0]);
%else
bar(modelSigma2(order)', 0.7); colormap([0.9 0.9 0.9]);
%end
hold on;
errorbar([1:NumOfGenes]-0.14, modelSigma2(order), modelSigma2(order)-stdSigma2_1(order), stdSigma2_2(order)-modelSigma2(order),'.');
set(gca,'fontsize',FS2);
axis tight;
%errorbar([1:NumOfGenes], modelA(order), modelA(order)-stdAA1(order), stdAA2(order)-modelA(order),'.');
%title('W0','fontsize', FONTSIZE);
if printResults
      print('-depsc', [dirr fileName 'Sigma2']);
end
title('Sigma2','fontsize', FS2);
end


% plot delayes if they were inferred
if abs(model.Likelihood.tauMax) > 0
Taus = samples.Taus';
modelTaus = mean(Taus,1);
stdTaus_1 = prctile(Taus,5,1);%sqrt(var(Taus));
stdTaus_2 = prctile(Taus,95,1);%sqrt(var(Taus));
figure;
% plot initial conditions
if isfield(model,'groundtr') == 1
bar([modelTaus(order); model.groundtr.Taus]', 0.7); colormap([0.9 0.9 0.9; 0 0 0]);
else
bar(modelTaus(order)', 0.7); colormap([0.9 0.9 0.9]);
end
hold on;
errorbar([1:NumOfGenes]-0.14, modelTaus(order), modelTaus(order)-stdTaus_1(order), stdTaus_2(order)-modelTaus(order),'.');
set(gca,'fontsize',FS2);
%errorbar([1:NumOfGenes], modelA(order), modelA(order)-stdAA1(order), stdAA2(order)-modelA(order),'.');
%title('W0','fontsize', FONTSIZE);
if printResults
      print('-depsc', [dirr fileName 'Delays']);
end
title('Delays','fontsize', FS2);
end
