% TOY DATA - ONE CONDITION
%
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
     
figure; 
plot(TimesFP,exp(TFs(1,:,2)), 'b','lineWidth',2);
hold on; 
plot(TimesFP,exp(TFs(2,:,2)), 'r','lineWidth',2);
plot(TimesFP,exp(TFs(3,:,2)), 'g','lineWidth',2);
plot(TimesFP,exp(TFs(4,:,2)), 'k-.','lineWidth',2);
axis tight;
set(gca, 'FontSize',10);
set(gca, 'YTickLabel', []);
xlabel('time');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 4.5 3.5]);
fileName = 'groudTruthTFs_Rep2';
print('-depsc',[dirr fileName]); 
print('-dpng', [dirr fileName]);

%
clear;
close all;
load datasets/toy4TFs28_June_10.mat;
load demtoy_dataOneCond108-Jul-2010.mat
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
demToyTrain_3TFs30GenesTwoCondsImData1;
load datasets/toy4TFs28_June_10.mat;
load demtoy_dataTwoConds109-Jul-2010.mat
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
%demDrosophilaTrain_5TFs92GenesImData1;
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
gpmtfPlot(model, samples, 'dros', drosTF.names, fbgns, 1, 'figures/', ops);
% Plot also genes parameters
ops(7)=0;
% bacth plots
ops(2)=0;
gpmtfPlot(model, samples, 'dros', drosTF.names, fbgns, 1, 'figures/', ops); 
close all;

% produce drosophila bars 
addpath ../../disimrank/matlab/;
load datasets/drosophila_data.mat
analyseDrosophila_joint_2010_11_24;
analyseDrosophila_combinations_2010_11_24;
clear;
analyseToy_2010_07_29;