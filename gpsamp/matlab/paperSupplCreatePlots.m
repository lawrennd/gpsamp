% TOY DATA - ONE CONDITION
%
ddpath ../../disimrank/matlab/;
addpath toolbox/;
addpath activFuncts/;
clear;
demToyTrain_3TFs30GenesOneCondImData1;
% plot ground truth TF profiles
TimesFP = model.Likelihood.TimesF;
plot(TimesFP,exp(TFs(1,:,1)), 'b','lineWidth',2);
hold on; 
plot(TimesFP,exp(TFs(2,:,1)), 'r','lineWidth',2);
plot(TimesFP,exp(TFs(3,:,1)), 'g','lineWidth',2);
plot(TimesFP,exp(TFs(4,:,1)), 'k-.','lineWidth',2);
axis tight;
set(gca, 'FontSize',10);
set(gca, 'YTickLabel', []);
xlabel('time');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 4.5 3.5]);
dirr = 'figures/';
fileName = 'groudTruthTFs_Rep1';
print('-depsc',[dirr fileName]); 
print('-dpng', [dirr fileName]);
     
% % plot ground truth TF profiles
TimesFP = model.Likelihood.TimesF;
plot(TimesFP,exp(TFs(1,:,1)), 'b','lineWidth',2);
hold on; 
plot(TimesFP,exp(TFs(2,:,1)), 'r','lineWidth',2);
plot(TimesFP,exp(TFs(3,:,1)), 'g','lineWidth',2);
plot(TimesFP,exp(TFs(4,:,1)), 'k-.','lineWidth',2);
axis tight;
set(gca, 'FontSize',10);
set(gca, 'YTickLabel', []);
xlabel('time');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 4.5 3.5]);
dirr = 'figures/';
fileName = 'groudTruthTFs_Rep1';
print('-depsc',[dirr fileName]); 
print('-dpng', [dirr fileName]);


% plot ground TF mRNA functions
figure; 
load mRNAToyGPF.mat;
plot(TimesFP,exp(F(1,:,1)), 'b','lineWidth',2);
hold on; 
plot(TimesFP,exp(F(2,:,1)), 'r','lineWidth',2);
plot(TimesFP,exp(F(3,:,1)), 'g','lineWidth',2);
plot(TimesFP,exp(outlierF), 'k-.','lineWidth',2);
axis tight;
set(gca, 'FontSize',10);
set(gca, 'YTickLabel', []);
%set(gca, 'XTickLabel', TimesG);
xlabel('time');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 4.5 3.5]);
fileName = 'groudTruthTFmRNAs_Rep1';
print('-depsc',[dirr fileName]); 
print('-dpng', [dirr fileName]);
figure; 
plot(TimesFP,exp(F(1,:,2)), 'b','lineWidth',2);
hold on; 
plot(TimesFP,exp(F(2,:,2)), 'r','lineWidth',2);
plot(TimesFP,exp(F(3,:,2)), 'g','lineWidth',2);
plot(TimesFP,exp(outlierF), 'k-.','lineWidth',2);
axis tight;
set(gca, 'FontSize',10);
set(gca, 'YTickLabel', []);
%set(gca, 'XTickLabel', TimesG);
xlabel('time');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 4.5 3.5]);
fileName = 'groudTruthTFmRNAs_Rep2';
print('-depsc',[dirr fileName]); 
print('-dpng', [dirr fileName]);


%
clear;
close all;
load datasets/toy4TFs28_June_10.mat;
%load demtoy_dataOneCond108-Jul-2010.mat
load demtoy_dataOneCond127-Jun-2011.mat
model.groundtr.kineticsTF = [0.99461 1 0; 0.9457 1 0; 0.640 1 0];
tmp = zeros(3,3,1000);
tmp(:,1,:) = samples.kineticsTF(:,1,:);
tmp(:,2,:) = samples.kineticsTF(:,2,:);
samples.kineticsTF = tmp;
for j=1:size(model.Likelihood.Genes,1)
   model.groundtr.kinetics(j,:) = GrTruth{j}.kinetics;
   model.groundtr.W(j,:) = GrTruth{j}.W;
   model.groundtr.W0(j) = GrTruth{j}.W0;
end
model.groundtr.W0 = model.groundtr.W0';
%model.groundtr.TF = exp(TFs(1:3,:,:));
TFnames = {'ANT', 'BEE', 'CAR', 'UNK'};
%timeshift 
ops(1) = 0;
% separate plots or grouped (1:separate, 0:not separate) 
ops(2) = 1;
% allow different colours fot the TF profiles (!: not allow, 0:allow)
ops(3)=0;
% write replicas or conditions (1:rep , 0:cond)
ops(4)=0;
% write the wrod rep or cond or not (1: no print)
ops(5)=1;
ops(6) = 0;
% return after TF plots
ops(7) = 1;
% time (0:in nuetral units, 1: in hours, 2:in seconds)
ops(8)=0;
gpmtfPlot(model, samples, 'toyOneCond', TFnames, 1:30, 1, 'figures/', ops);
close all;
% Plot also genes parameters
ops(7)=0;
% bacth plots
ops(2)=0;
% time (0:in nuetral units, 1: in hours, 2:in seconds)
ops(8)=0;
gpmtfPlot(model, samples, 'toyOneCond', TFnames, 1:30, 1, 'figures/', ops);


% TOY DATA - TWO CONDITIONS
%
%
clear;
close all;
%demToyTrain_3TFs30GenesTwoCondsImData1;
demToyTrain_3TFs30GenesTwoCondsImData2;
load datasets/toy4TFs28_June_10.mat;
%load demtoy_dataTwoConds109-Jul-2010.mat
load demtoy_dataTwoConds127-Jun-2011.mat
model.groundtr.kineticsTF = [0.99461 1 0; 0.9457 1 0; 0.640 1 0];
tmp = zeros(3,3,1000);
tmp(:,1,:) = samples.kineticsTF(:,1,:);
tmp(:,2,:) = samples.kineticsTF(:,2,:);
samples.kineticsTF = tmp;
for j=1:size(model.Likelihood.Genes,1)
   model.groundtr.kinetics(j,:) = GrTruth{j}.kinetics;
   model.groundtr.W(j,:) = GrTruth{j}.W;
   model.groundtr.W0(j) = GrTruth{j}.W0;
end
model.groundtr.W0 = model.groundtr.W0';
%model.groundtr.TF = TFs(1:3,:,:);
TFnames = {'ANT', 'BEE', 'CAR', 'UNK'};
%timeshift 
ops(1) = 0;
% separate plots or grouped (1:separate, 0:not separate) 
ops(2) = 1;
% allow different colours fot the TF profiles (!: not allow, 0:allow)
ops(3)=0;
% write replicas or conditions (1:rep , 0:cond)
ops(4)=0;
% write the wrod rep or cond or not (1: no print)
ops(5)=1;
% which replicas are plotted in for the mRNA data (1:all, 0:one, 2: first 2)
ops(6) = 2;
% return after TF plots
ops(7) = 1;
% time (0:in nuetral units, 1: in hours, 2:in seconds)
ops(8)=0;
gpmtfPlot(model, samples, 'toyTwoCond', TFnames, 1:30, 1, 'figures/', ops);
close all;
% Plot also genes parameters
ops(7)=0;
% batch plots
ops(2)=0;
% time (0:in neutral units, 1: in hours, 2:in seconds)
ops(8)=0;
gpmtfPlot(model, samples, 'toyTwoCond', TFnames, 1:30, 1, 'figures/', ops);


% DROSOPHILA
%
%
clear;
close all;
% makes the run with  92  genes  
% demDrosophilaTrain_5TFs92GenesImData1;   
% Filters out 25 genes and re-runs using 25 genes 
% demDrosophilaTrain_5TFs92GenesImData6;  
%load drosTrainTotal.mat;
load drosTrainTotal_14DEC2010.mat
%samples.W = samples.Weights;
%samples.W0 = samples.Weights0;
%timeshift 
ops(1) = 1;
% separate plots or grouped (1:separate, 0:not separate) 
ops(2) = 1;
% allow different colours fot the TF profiles (!: not allow, 0:allow)
ops(3)=1;
% write replicas or conditions (1:rep , 0:cond)
ops(4)=1;
% write the wrod rep or cond or not (1: no print)
ops(5)=1;
% which replicas are plotted in for the mRNA data (1:all, 0:one) 
ops(6) = 0;
% return after TF plots
ops(7) = 1;
% time (0:in neutral units, 1: in hours, 2:in seconds)
ops(8)=1;
drosTF.names = {'TIN', 'BIN', 'TWI', 'BAP', 'MEF2'};
gpmtfPlot(model, samples, 'dros', drosTF.names, fbgns, 1, 'figures/', ops);
% return after TF mRNA plots and do not show vertical ticks in the plots
ops(7) = 2;
gpmtfPlot(model, samples, 'dros', drosTF.names, fbgns, 1, 'figures/', ops);
% Plot also genes parameters
ops(7)=0;
% batch plots
ops(2)=0;
gpmtfPlot(model, samples, 'dros', drosTF.names, fbgns, 1, 'figures/', ops); 
close all;

clear;
ops(1) = 1;ops(2) = 0; ops(3)=1; ops(4)=1; ops(5)= 1; ops(6) = 0; ops(8)=1; ops(7)=0;  ops(2)=0;
thr = 0.01;
% load a training phase with puma+white corresponding to 92 genes
% Plot only the genes that are not fitted well by the model (i.e the 67 genes)
load drosTrainTotal_28NOV2010.mat;
if ops(8) == 0, XLabel = 'time'; elseif ops(8) == 1, XLabel = 'time (h)'; else  XLabel = 'time (s)'; end; 
FONTSIZE=10; FS2 = 16; Genename = fbgns; printResults = 1;  demdata = 'dros'; ok= ''; rep = 'rep';
stdSigma2_2 = prctile(samples.sigma2',95,1);
selectGenes = (stdSigma2_2<thr); 
SepPlot = 0; ManyGenes = 5; NumOfSamples = size(samples.LogL,2); timeshift = ops(1);
Genes = model.Likelihood.Genes;  TimesF = model.Likelihood.TimesF; TimesG = model.Likelihood.TimesG; TimesFF = TimesF(model.Likelihood.startTime:end);
GeneVars = model.Likelihood.noiseModel.pumaSigma2;
sRep = [1 0 0];NumOfGenes=92; Wrtrep = 0; NumOfReplicas=3; jj = 0;
fileName = [demdata 'MCMC' ok model.Likelihood.singleAct model.Likelihood.jointAct]; dirr = 'figures/';
list = 1:92;
for j=1:NumOfGenes
  ok=0;  
  if  (selectGenes(j) ==0 & ismember(j,list)==1),  jj =  jj + 1;  ok=1; end;
  if  (ok ==1  & SepPlot==0), if mod(jj,ManyGenes) == 1,  figure; end; end;
  if  ok ==1  
  rrc = 0;
  for r=1:NumOfReplicas 
     if sRep(r) == 1
     rrc = rrc + 1;    
     GG = zeros(NumOfSamples,model.Likelihood.sizTime);    
     for t=1:NumOfSamples
         LikParams = model.Likelihood;    LikParams.kinetics = samples.kinetics(:,:,t);   LikParams.kineticsTF = samples.kineticsTF(:,:,t);
         LikParams.W = samples.W(:,:,t);   LikParams.W0 = samples.W0(:,t);
         predgen = gpmtfComputeGeneODE(LikParams, samples.F{t}(:,:,r), r, j);
         GG(t,:) = predgen;
     end
     mu = mean(GG,1)';   stds = sqrt(var(GG,0,1))';
     TF = TimesFF';      TFP = TF + timeshift; 
     TimesGP = TimesG + timeshift; 
     if SepPlot==0, subplot(3, ManyGenes, (rrc-1)*ManyGenes + 1 + mod(jj-1, ManyGenes)); else figure; end;
     plot(TFP,mu,'b','lineWidth',2);     hold on;
     fillColor = [0.7 0.7 0.7];
     fill([TFP; TFP(end:-1:1)], [mu; mu(end:-1:1)] + 2*[stds; -stds(end:-1:1)], fillColor,'EdgeColor',fillColor);
     plot(TFP,mu,'b','lineWidth',2);
     plot(TimesGP,Genes(j,:,r),'rx','markersize', 14','lineWidth', 2);
     if model.Likelihood.noiseModel.active(1) == 1, errorbar(TimesGP,  Genes(j,:,r), 2*sqrt(GeneVars(j,:,r)), 'rx','lineWidth', 1.5); end;
     axis tight;     set(gca, 'FontSize', FONTSIZE);
     %set(gca, 'YTickLabel', []);    
     xlabel(XLabel);
     jj
     if SepPlot==1,   set(gcf, 'PaperUnits', 'centimeters');      set(gcf, 'PaperPosition', [0 0 4.5 3.5]); end;
     titlestring = [];
     if isnumeric(Genename(j)), titlestring = [titlestring, num2str(j)]; else  titlestring = [titlestring, Genename{j}]; end;
     if NumOfReplicas > 1,  if Wrtrep == 1,   titlestring = [titlestring, ', ' rep, ' ', num2str(r)]; end;  end;
     title(titlestring,'fontsize', FONTSIZE);
     if SepPlot==0 
        if printResults && r==sum(sRep) && (mod(jj,ManyGenes)==0 || jj==NumOfGenes),
        if NumOfReplicas > 1
           print('-depsc', [dirr fileName 'Replica' num2str(sRep,'%1g') 'GENEmRNA_FilteredOut' num2str(jj)]);
        else
           print('-depsc', [dirr fileName 'GENEmRNA_FilteredOut' num2str(jj)]); 
        end
        end
     else
        print('-depsc', [dirr fileName 'Replica' num2str(sRep,'%1g') 'GENEmRNA_FilteredOut' num2str(jj)]);
     end
  end
  end
  end
end
clear;
close all; 




% produce drosophila bars 
addpath ../../disimrank/matlab/;
load datasets/drosophila_data.mat
analyseDrosophila_joint_2010_11_24;
analyseDrosophila_combinations_2010_11_24;
clear;
analyseToy_2010_07_29;


% Compute a table that summarizes the  values of the ground truth kinetic 
% parameters for the toy example 
clear;
close all;
load datasets/toy4TFs28_June_10.mat;
Kin = [];
for j=1:1030
    Kin = [Kin; GrTruth{j}.kinetics];
end
med = prctile(Kin,50);
pr5 = prctile(Kin,5);
pr95 = prctile(Kin, 95);


% Compute a table that summarizes the values of kinetics parameters found
% in the drosophila 25 training genes 
clear;
close all;
load drosTrainTotal_14DEC2010.mat
DD = squeeze(samples.kineticsTF(:,1,:))';
SS = squeeze(samples.kineticsTF(:,2,:))';
AA = squeeze(samples.kineticsTF(:,3,:))';
medS = median(SS,1);
medD = median(DD,1);
medA = median(AA,1);
prD5 = prctile(DD,5,1);
prS5 = prctile(SS,5,1);
prA5 = prctile(AA,5,1);
prD95 = prctile(DD,95,1);
prS95 = prctile(SS,95,1);
prA95 = prctile(AA,95,1);
DD = squeeze(log(2)./samples.kineticsTF(:,1,:))';
SS = squeeze(log(2)./samples.kineticsTF(:,2,:))';
AA = squeeze(log(2)./samples.kineticsTF(:,3,:))';
HmedS = median(SS,1);
HmedD = median(DD,1);
HmedA = median(AA,1);
HprD5 = prctile(DD,5,1);
HprS5 = prctile(SS,5,1);
HprA5 = prctile(AA,5,1);
HprD95 = prctile(DD,95,1);
HprS95 = prctile(SS,95,1);
HprA95 = prctile(AA,95,1);
disp('degradation estimates, eah row corresponds to protein')
log(2)./[prD5; medD; prD95]'
disp('degradation estimates, eah row corresponds to protein')
[HprD5; HmedD; HprD95]'





