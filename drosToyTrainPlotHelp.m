
clear;
demToyTrain_3TFs30GenesOneCondImData1;
load datasets/toy4TFs28_June_10.mat;
load demtoy_dataOneCond108-Jul-2010.mat
model.groundtr.kineticsTF = [0.99461 1; 0.9457 1; 0.640 1];
for j=1:size(model.Likelihood.Genes,1)
   model.groundtr.kinetics(j,:) = GrTruth{j}.kinetics;
   model.groundtr.W(j,:) = GrTruth{j}.W;
   model.groundtr.W0(j) = GrTruth{j}.W0;
end
model.groundtr.W0 = model.groundtr.W0';
model.groundtr.TF = TFs(1:3,:,:);
TFnames = {'ANT', 'BEE', 'CAR', 'UNK'};
gpmtfPlot(model, samples, 'toyOneCond', TFnames, 1:30, 1, 'figures/');
% important argumetns in plot fucntino must be as follows
%rep = 'cond'; 
%Wrtrep = 0; 
%sRep = [1 0 0];

clear;
demToyTrain_3TFs30GenesTwoCondsImData1;
load datasets/toy4TFs28_June_10.mat;
load demtoy_dataTwoConds109-Jul-2010.mat
model.groundtr.kineticsTF = [0.99461 1; 0.9457 1; 0.640 1];
for j=1:size(model.Likelihood.Genes,1)
   model.groundtr.kinetics(j,:) = GrTruth{j}.kinetics;
   model.groundtr.W(j,:) = GrTruth{j}.W;
   model.groundtr.W0(j) = GrTruth{j}.W0;
end
model.groundtr.W0 = model.groundtr.W0';
model.groundtr.TF = TFs(1:3,:,:);
TFnames = {'ANT', 'BEE', 'CAR', 'UNK'};
gpmtfPlot(model, samples, 'toyTwoCond', TFnames, 1:30, 1, 'figures/');
% important argumetns in plot fucntino must be as follows
%rep = 'cond'; 
%Wrtrep = 1; 
%sRep = [1 1 0];

demDrosophilaTrain_5TFs92GenesImData1;
load drosTrainTotal.mat;
samples.W = samples.Weights;
samples.W0 = samples.Weights0;
gpmtfPlot(model, samples, 'dros', drosTF.names, fbgns, 1, 'figures/');
% important argumetns in plot fucntino must be as follows
%rep = 'rep'; 
%Wrtrep = 0; 
%sRep = [1 0 0];