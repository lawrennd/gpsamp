% USER-specified: Do you want constrained (positive interactions weights)
% or unconstrained models? 
% !! Warning: Baseline models was run only for the constrained case !!
flag = 1; % "1" for  constrained; "anything else" for unconstrained 

% USER-specified: Sizes of the ranking sets for the first plot
T1 = [20 50 100 200, 400, 800, 1600 3200];
T2 = [10 20 100 200, 400, 800, 915];

% USER-specified:: Directory where you store figures
ddir = 'figures/';
printPlot = 1; % 0 means not printing

% plot MAP models
plotMAP = [1 1 1]; 

% USER-specified: whether to include empirical prior line
figSize = [7 5];
fontSize = 7;

% models (as stored  in the results variables; see below) 
% correspidng to 5 TFs being active/inactive 
combConstr = [0 0 0 0 0;
	1 0 0 0 0; 0 1 0 0 0; 0 0 1 0 0; 0 0 0 1 0; 0 0 0 0 1;
	1 1 0 0 0; 1 0 1 0 0; 1 0 0 1 0; 1 0 0 0 1;
	0 1 1 0 0; 0 1 0 1 0; 0 1 0 0 1;
	0 0 1 1 0; 0 0 1 0 1;
	0 0 0 1 1;
	1 1 1 0 0; 1 1 0 1 0; 1 1 0 0 1;
	1 0 1 1 0; 1 0 1 0 1;
	1 0 0 1 1;
	0 1 1 1 0; 0 1 1 0 1;
	0 1 0 1 1;
	0 0 1 1 1;
	1 1 1 1 0; 1 1 1 0 1; 1 1 0 1 1; 1 0 1 1 1; 0 1 1 1 1;
	1 1 1 1 1];

combUnconstr = [0 0 0 0 0;
	1 0 0 0 0; 0 1 0 0 0; 0 0 1 0 0; 0 0 0 1 0; 0 0 0 0 1;
	1 1 0 0 0; 1 0 1 0 0; 1 0 0 1 0; 1 0 0 0 1;
	0 1 1 0 0; 0 1 0 1 0; 0 1 0 0 1;
	0 0 1 1 0; 0 0 1 0 1;
	0 0 0 1 1;
	1 1 1 0 0; 1 1 0 1 0; 1 1 0 0 1;
	1 0 1 1 0; 1 0 1 0 1;
	1 0 0 1 1;
	0 1 1 1 0; 0 1 1 0 1;
	0 1 0 1 1;
	0 0 1 1 1;
	1 1 1 1 0; 1 1 1 0 1; 1 1 0 1 1; 1 0 1 1 1; 0 1 1 1 1;
	1 1 1 1 1];

baselinecomb = [0 0 0 0 0;
		1 0 0 0 0; 0 1 0 0 0; 0 0 1 0 0 ; 0 0 0 1 0; 0 0 0 0 1;
		1 1 0 0 0; 1 0 1 0 0; 1 0 0 1 0; 1 0 0 0 1;
		0 1 1 0 0; 0 1 0 1 0; 0 1 0 0 1;
		0 0 1 1 0; 0 0 1 0 1;
		0 0 0 1 1];

    
if flag == 1     
    % load the results
    %results_b = sortResults(load('results/multitf8b_2010-12-06_summary.mat'));
    results_b = sortResults(load('results/multitf8c_2010-12-14_summary.mat'));   
    % load the 15 baseline models (zeroth model is excluded)
    baseline_a = sortResults(load('results/multitf5a_baseline_2010-07-02_summary.mat'));
    % load the zeroth model and append
    baseline_a2 = sortResults(load('results/multitf5a_baseline2_2010-07-02_summary.mat'));
    baseline_a.marlls = [baseline_a2.marlls, baseline_a.marlls];
else 
    % load the unconstrained summaries 
    results_b = sortResults(load('results/multitf7a_2010-07-13_summary.mat'));    
    % !! For Baseline uses the constrained models (if we train the uncosntrained ones, modify the code below) !!
    % load the 15 baseline models (zerhoth model is excluded)
    baseline_a = sortResults(load('results/multitf5a_baseline_2010-07-02_summary.mat'));
    % load the zeroth model and append
    baseline_a2 = sortResults(load('results/multitf5a_baseline2_2010-07-02_summary.mat'));
    baseline_a.marlls = [baseline_a2.marlls, baseline_a.marlls];    
    % the code below refers to the combConstr variable, so re-define it
    combConstr = combUnconstr; 
end


load results/res_linearPositive.mat 
regression_a.marlls = -predErrors1;
%regression_a.marlls = -predErrors1norm;
regression_a.genes = results_b.genes; 


% You need to exclude the 92 training genes 
genesAndChip = importdata('datasets/eileen_nature_training_set.txt'); 
Trfbgns = genesAndChip.textdata(2:end,1);
mask = ones(size(results_b.genes,1),1); 
for i=1:size(results_b.genes,1)
    for j=1:size(Trfbgns,1)
         if strcmp(results_b.genes(i), Trfbgns(j)) 
             mask(i) = 0; 
         end
    end
end
regression_a.marlls = regression_a.marlls(mask==1,:);
regression_a.genes = regression_a.genes(mask==1); 
results_b.marlls = results_b.marlls(mask==1,:);
results_b.genes = results_b.genes(mask==1); 
baseline_a.marlls = baseline_a.marlls(mask==1,:);
baseline_a.genes = baseline_a.genes(mask==1);

% number of TFs
numTFs = size(combConstr,2);
M = drosMakeValidationMatrix(chipdistances, results_b.genes, 2000);
numGenes = size(M,1); 

% indices of all 32 models 
ind32 = 1:size(results_b.marlls, 2); 
posteriors = {};
% numerically stable computed posterior over each model from the 32 models 
for k=1:length(ind32)
    % posterior of each model in the 32-model 
    posteriors{1}(:,k) = results_b.marlls(:,k) - logsumexp(results_b.marlls(:, setdiff(ind32, k)), 2);
end

% baseline models likelihood ratios 
posteriors{end+1} = baseline_a.marlls - repmat(baseline_a.marlls(:, 1), ...
					   [1, size(baseline_a.marlls, 2)]);
posteriors{end+1} = regression_a.marlls;
        
% Compute global ranking performance using random prediction.
accRand = 0;
countRand = 0;
for gg =1:1000
% generate a random network
Mrand =  round(rand(size(M,1),5));
for i=1:numGenes,
    if ~isnan(sum( M(i,:) ))
      accRand = accRand + 1 - sum(abs(Mrand(i,:) - M(i,:)))/5;
      countRand = countRand + 1;
    end
end
end
prioraccs31 = accRand/countRand;

% compute focussed ranking performance using random prediction. Compute
[C, IA, IB] = intersect(drosinsitu.genes(any(drosinsitu.data, 2)), results_b.genes);
accRand = 0;
countRand = 0;
accPrior = 0;
countPrior = 0;
N = length(C);
for gg=1:1000
% generate a random network
Mrand =  round(rand(N,5));
for i=1:N
    if ~isnan(sum( M(i,:) ))
    accRand = accRand + sum(abs(Mrand(i,:) - M(IB(i),:)))/5;
    countRand = countRand + 1;
    end
end
end
focused_prioraccs31 =  accRand/countRand;

                               
% 1 PLOT ---------------------------------------------------------------
% GLOBAL RANKING BASED ON THE MAP MODEL 
% START  ---------------------------------------------------------------             
scores = {}; I={}; sscores={}; J={};
r1 = zeros(length(posteriors), length(T1));
pvals1 = r1;
for k=1:length(posteriors)   
  if ((k==1) | (k==3))
    mycomb = combConstr;
  elseif  (k==2)
    mycomb = baselinecomb;
  end
  [scores{k}, I{k}] = max(posteriors{k}, [], 2);
  [sscores{k}, J{k}] = sort(scores{k}, 'descend');
  for l=1:length(T1),
    acc = 0;
    count = 0;
    for m=1:T1(l)
        if ~isnan(sum(  M(J{k}(m), :)  ))   
	      acc = acc + 1 - sum(abs( M(J{k}(m), :) - mycomb(I{k}(J{k}(m)), :)))/5;   
	      count = count + 1;
        end
    end
    r1(k, l) = acc / count;
    pvals1(k, l) = 1 - binocdf(acc - 1, count, prioraccs31);
  end
end
h1 = figure;
set(gca, 'FontSize', fontSize);
% plot bars
h = bar(100*r1(plotMAP==1,:)');
set(gca, 'XTickLabel', T1);
hold on
v = axis;
v(3:4) = [0 100];
axis(v);
plot(v(1:2), 100*prioraccs31*[1 1], 'b');
hold off
legends = {'MAP-32', 'Baseline', 'Regression', 'Random', 'Location', 'EastOutside'};
legend(legends([1 2 3 4]), 'Location', 'SouthEast');
axis(v)
xlabel('# of global top genes')
ylabel('Accuracy (%)')
drosStarBars(h, pvals1(plotMAP==1,:)');
set(gca, 'FontSize', fontSize);
set(gcf, 'PaperUnits', 'centimeters')
set(gcf, 'PaperPosition', [0, 0, figSize])
% 1 PLOT ---------------------------------------------------------------
% GLOBAL RANKING BASED ON THE MAP MODEL 
% END    ---------------------------------------------------------------



% 1B PLOT ---------------------------------------------------------------
% GLOBAL RANKING BASED ON THE MAP MODEL, IN-SITU FILTERING
% START  ---------------------------------------------------------------
% Find the genes with positive in situ annotations
[C, IA, IB] = intersect(drosinsitu.genes(any(drosinsitu.data, 2)), results_b.genes);
% In situ genes
Minsitu = M(IB,:);
scores = {}; I={}; sscores={}; J={};
r1 = zeros(length(posteriors), length(T2));
pvals1 = r1;
for k=1:length(posteriors)   
  if ((k==1) | (k==3))
    mycomb = combConstr;
  elseif  (k==2)
    mycomb = baselinecomb;
  end
  [scores{k}, I{k}] = max(posteriors{k}(IB,:), [], 2);
  [sscores{k}, J{k}] = sort(scores{k}, 'descend');
  for l=1:length(T2),
    acc = 0; acc2 = 0;
    count = 0;
    for m=1:T2(l), 
      if ~isnan(sum( Minsitu(J{k}(m), :) ))
          acc = acc +  1 - sum(abs( Minsitu(J{k}(m), :) - mycomb(I{k}(J{k}(m)), :) ))/5;   
	      count = count + 1;
      end
    end
    r1(k, l) = acc / count;
    pvals1(k, l) = 1 - binocdf(acc - 1, count, focused_prioraccs31);
  end
end
h1b = figure;
set(gca, 'FontSize', fontSize);
% plot bars
h = bar(100*r1(plotMAP==1,:)');
set(gca, 'XTickLabel', T2);
hold on
v = axis;
v(3:4) = [0 100];
axis(v);
plot(v(1:2), 100*focused_prioraccs31*[1 1], 'b');
hold off
legends = {'MAP-32',  'Baseline', 'Regression', 'Random', 'Location', 'EastOutside'};
legend(legends([1 2 3 4]), 'Location', 'SouthEast');
axis(v)
xlabel('# of global top genes')
ylabel('Accuracy (%)');
drosStarBars(h, pvals1(plotMAP==1,:)');
set(gca, 'FontSize', fontSize);
set(gcf, 'PaperUnits', 'centimeters')
set(gcf, 'PaperPosition', [0, 0, figSize])
% 1B PLOT ---------------------------------------------------------------
% GLOBAL RANKING BASED ON THE MAP MODEL, IN-SITU FILTERING
% END    ---------------------------------------------------------------


% print Plots
property = 'Constrained';
if flag ~= 1
    property = 'Unconstrained'; 
end
if printPlot 
   print(h1, '-depsc2', [ddir 'drosophilaBars_' 'GlobalHammingMAP' num2str(2) '_', property '.eps']);
   print(h1b, '-depsc2', [ddir 'drosophilaBars_' 'GlobalHammingMAPIN_SITU' num2str(2) '_', property '.eps']);
end
