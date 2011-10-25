% USER-specified: Do you want constrained (positive interactions weights)
% or unconstrained models? 
% !! Warning: Baseline models was run only for the constrained case !!
flag = 1; % "1" for  constrained; "anything else" for unconstrained 

% USER-specified: Sizes of the ranking sets for the first plot
T1 = [200, 400, 800, 1600, 3200, 6000];
% Sizes of the ranking sets for the second plot
T2 = [100, 200, 400, 800, 1600, 3200];

% USER-specified:: Directory where you store figures
ddir = 'figures/';
printPlot = 1; % 0 means not printing

% plot MAP models
plotMAP = [1 0 0 0 1 1]; 
%plotMAP = [1 1 0 0 0]; 
% for the rest plots 
plotRest = [1 0 0 0 1 0 1 1]; 
%plotRest = [1 1 0 0 1 1 0]; 

% USER-specified: whether to include empirical prior line
incPrior = 0;
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
    
    % load the 15 baseline models (zerhoth model is excluded)
    baseline_a = sortResults(load('results/multitf5a_baseline_2010-07-02_summary.mat'));
    % load the zeroth model and append
    baseline_a2 = sortResults(load('results/multitf5a_baseline2_2010-07-02_summary.mat'));
    baseline_a.marlls = [baseline_a2.marlls, baseline_a.marlls];
else 
    % load the unconstrained summaries 
    results_b = sortResults(load('results/multitf7a_2010-07-13_summary.mat'));
    
    % !! For Baseline uses the constrained models (if we train the uncosntrained ones, modify the code below) !!
    % load the 15 base inferelatorline models (zerhoth model is excluded)
    baseline_a = sortResults(load('results/multitf5a_baseline_2010-07-02_summary.mat'));
    % load the zeroth model and append
    baseline_a2 = sortResults(load('results/multitf5a_baseline2_2010-07-02_summary.mat'));
    baseline_a.marlls = [baseline_a2.marlls, baseline_a.marlls];
    
    % the code below refers to the combConstr variable, so re-define it
    combConstr = combUnconstr; 
end

%baselinecomb = combConstr;

load results/res_linearPositive.mat 
regression_a.marlls = -predErrors1;
%regression_a.marlls = -predErrors1norm;
regression_a.genes = results_b.genes; 

% load inferelator results
load inferelator_models.mat; 



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
% mask out the nans
mask = ones(size(results_b.genes,1),1); 
for i=1:size(results_b.genes,1)
    if isnan(sum( M(i,:) ))
       mask(i) = 0;
    end
end
regression_a.marlls = regression_a.marlls(mask==1,:);
regression_a.genes = regression_a.genes(mask==1); 
results_b.marlls = results_b.marlls(mask==1,:);
results_b.genes = results_b.genes(mask==1); 
baseline_a.marlls = baseline_a.marlls(mask==1,:);
baseline_a.genes = baseline_a.genes(mask==1);
singles = singles(mask==1,:);
doubles = doubles(mask==1,:);
M = M(mask==1,:);
numGenes = size(M,1); 


% computation of the prior over 32 models
prior32 = zeros(1, size(combConstr, 1));
for k=1:size(combConstr, 1),
  indOne = find(combConstr(k, :)==1);
  %prior(k)=0;
  for j=1:numGenes
      if ~isnan(prod(M(j,:))) 
      if all(M(j, indOne),2)  & sum(M(j, setdiff(1:numTFs, indOne)),2)==0
          prior32(k) = prior32(k) + 1;
      end
      end
  end
end
% add smoothing/Laplace factor 
alpha = 1;
prior32 = prior32 + 1;
prior32 = prior32/sum(prior32); 
% end of computing the prior 

% prior over all models excluding the zeroth model 
prior31 = prior32(2:end)/sum(prior32(2:end));

% Bayes correct classification rate based on a classifier that  uses
% a unifrom prior (as classification rule)
%prioraccs31(1) = 1/size(prior31,2); 
% Bayes correct classification rate based on a classifier that uses 
% the empirical prior
%prioraccs32(2) = sum(prior32(2:end).*prior32(2:end),2)/sum(prior32(2:end)); 
%prioraccs31(2) = sum(prior31.*prior31,2); 


    
% indices of all 32 models 
ind32 = 1:size(results_b.marlls, 2); 
% indices that correspond to 16 models (at most 2 TFs) 
ind16 = find(sum(combConstr,2)<=2); 
% index of the zeroth model 
ind0 = find(sum(combConstr,2)==0);
linkMargPosteriors = {};
comb16 = combConstr(ind16,:);


posteriors = {};
% non-uniform prior over models 
withprior = results_b.marlls + repmat(log(prior32), [size(results_b.marlls, 1), 1]);
% numerically stable computed posterior over each model from the 32 models 
for k=1:length(ind32)
    % posterior of each model in the 32-model 
    posteriors{1}(:,k) = results_b.marlls(:,k) - logsumexp(results_b.marlls(:, setdiff(ind32, k)), 2);
    % posterior of each model in the 32-model case plus prior 
    posteriors{2}(:,k) = withprior(:,k) - logsumexp(withprior(:, setdiff(ind32, k)), 2);
end

% numerically stable computed posterior over each model from the 16 models 
for k=1:length(ind16)
    % posterior of each model in the 16-model 
    posteriors{3}(:,k) = results_b.marlls(:,ind16(k)) - logsumexp(results_b.marlls(:, setdiff(ind16, k)), 2);
    % posterior of each model in the 16-model case plus prior 
    posteriors{4}(:,k) = withprior(:,ind16(k)) - logsumexp(withprior(:,  setdiff(ind16, k)), 2);
end


prioraccinds = [1, 2, 1, 2, 1, 1, 1];
% baseline models likelihood ratios 
posteriors{end+1} = baseline_a.marlls - repmat(baseline_a.marlls(:, 1), ...
					   [1, size(baseline_a.marlls, 2)]);

posteriors{end+1} = regression_a.marlls;

          
% compute global ranking perforamcne using random prediction. Compute

accRand = 0;
countRand = 0;
accPrior = 0;
countPrior = 0;
for l=1:30000
    % random ranking and prediction 
    perm = randperm(size(results_b.genes,1)); 
    % random network
    Mrand = round(rand(1,5));
    Mprior = combConstr(sample_discrete(prior32),:);
    if sum(Mrand(1,:),2) > 0 
	   accRand = accRand + all( M(perm(1), Mrand==1), 2);
	   countRand = countRand + 1;
    end
    if sum(Mprior(1,:),2) > 0 
	   accPrior = accPrior + all( M(perm(1), Mprior==1), 2);
	   countPrior = countPrior + 1;
    end
end
prioraccs31(1) =  accRand/countRand;
prioraccs31(2) =  accPrior/countPrior;

% compute focussed ranking performance using random prediction. Compute
[C, IA, IB] = intersect(drosinsitu.genes(any(drosinsitu.data, 2)), results_b.genes);
accRand = 0;
countRand = 0;
accPrior = 0;
countPrior = 0;
N = length(C);
for l=1:30000
    % random ranking and prediction 
    perm = randperm(N);
    % random network
    Mrand = round(rand(1,5));
    Mprior = combConstr(sample_discrete(prior32),:);
    if sum(Mrand(1,:),2) > 0 
	   accRand = accRand + all( M(IB(perm(1)), Mrand==1), 2);
	   countRand = countRand + 1;
    end
    if sum(Mprior(1,:),2) > 0 
	   accPrior = accPrior + all( M(IB(perm(1)), Mprior==1), 2);
	   countPrior = countPrior + 1;
    end
end
focused_prioraccs31(1) =  accRand/countRand;
focused_prioraccs31(2) =  accPrior/countPrior;

                   
                   
% 1 PLOT ---------------------------------------------------------------
% GLOBAL RANKING BASED ON THE MAP MODEL 
% START  ---------------------------------------------------------------             
scores = {}; I={}; sscores={}; J={};
r1 = zeros(length(posteriors), length(T1));
pvals1 = r1;
for k=1:length(posteriors)   
  if ((k<=2) | (k==6)), % 32 models 
    mycomb = combConstr;
  elseif  (k>2 & k<=4), % 16 models 
    mycomb = comb16;
  elseif (k==5)
    mycomb = baselinecomb;
  end
  [scores{k}, I{k}] = max(posteriors{k}, [], 2);
  [sscores{k}, J{k}] = sort(scores{k}, 'descend');
  for l=1:length(T1),
    acc = 0; acc2 = 0;
    count = 0;
    for m=1:T1(l),
      if I{k}(J{k}(m)) ~= 1,
	      acc = acc + all( M(J{k}(m), mycomb( I{k}(J{k}(m)), :)==1), 2); 
          %acc2 = acc2 + all(M(J{k}(m), :) == mycomb(I{k}(J{k}(m)), :), 2);
	      count = count + 1;
      end
    end
    r1(k, l) = acc / count;
    %r2(k, l) = acc2 / count;
    pvals1(k, l) = 1 - binocdf(acc - 1, count, prioraccs31(prioraccinds(k)));
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
plot(v(1:2), 100*prioraccs31(1)*[1 1], 'b');
if incPrior,
  plot(v(1:2), 100*prioraccs31(2)*[1 1], 'g');
end
hold off
legends = {'MAP-32', 'MAP-32 + prior', 'MAP-16', 'MAP-16 + prior', 'Baseline', 'Regression', 'Uniform prior', 'Empirical prior', 'Location', 'EastOutside'};
legend(legends([plotMAP, 1, incPrior]==1));
axis(v)
xlabel('# of global top genes')
ylabel('Enrichment (%)')
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

scores = {}; I={}; sscores={}; J={};
r1 = zeros(length(posteriors), length(T1));
pvals1 = r1;
for k=1:length(posteriors)   
  if ((k <=2) | (k==6)), % 32 models 
    mycomb = combConstr;
  elseif  (k>2 & k<=4), % 16 models 
    mycomb = comb16;
  elseif (k==5) % 32 models
    mycomb = baselinecomb;
  end
  [scores{k}, I{k}] = max(posteriors{k}, [], 2);
  [sscores{k}, J{k}] = sort(scores{k}, 'descend');
  for l=1:length(T1),
    acc = 0; acc2 = 0;
    count = 0;
    for m=1:T1(l),
      if I{k}(J{k}(m)) ~= 1 && any(IB == J{k}(m)),
	      acc = acc + all( M(J{k}(m), mycomb( I{k}(J{k}(m)), :)==1), 2); 
          %acc2 = acc2 + all(M(J{k}(m), :) == mycomb(I{k}(J{k}(m)), :), 2);
	      count = count + 1;
      end
    end
    r1(k, l) = acc / count;
    %r2(k, l) = acc2 / count;
    pvals1(k, l) = 1 - binocdf(acc - 1, count, focused_prioraccs31(prioraccinds(k)));
  end
end
h1b = figure;
set(gca, 'FontSize', fontSize);
% plot bars
h = bar(100*r1(plotMAP==1,:)');
set(gca, 'XTickLabel', T1);
hold on
v = axis;
v(3:4) = [0 100];
axis(v);
plot(v(1:2), 100*focused_prioraccs31(1)*[1 1], 'b');
%plot(v(1:2), 100*prioraccs31(1)*[1 1], 'b--');
if incPrior,
  plot(v(1:2), 100*focused_prioraccs31(2)*[1 1], 'g');
  %plot(v(1:2), 100*prioraccs31(2)*[1 1], 'g--');
end
hold off
%legends = {'MAP-32', 'MAP-32 + prior', 'MAP-16', 'MAP-16 + prior', 'Baseline', 'Focused prior', 'Global prior', 'Focused Empirical prior', 'Global Empirical prior', 'Location', 'EastOutside'};
legends = {'MAP-32', 'MAP-32 + prior', 'MAP-16', 'MAP-16 + prior', 'Baseline', 'Regression', 'Uniform prior', 'Empirical prior', 'Location', 'EastOutside'};
%legend(legends([plotMAP, 1, 1, incPrior, incPrior]==1));
legend(legends([plotMAP, 1, incPrior]==1));
axis(v)
xlabel('# of global top genes')
ylabel('Enrichment (%)')
drosStarBars(h, pvals1(plotMAP==1,:)');
set(gca, 'FontSize', fontSize);
set(gcf, 'PaperUnits', 'centimeters')
set(gcf, 'PaperPosition', [0, 0, figSize])
% 1B PLOT ---------------------------------------------------------------
% GLOBAL RANKING BASED ON THE MAP MODEL, IN-SITU FILTERING
% END    ---------------------------------------------------------------


% Computation of marginal posterior probabilities over single links 
for k=1:numTFs
%
  % posterior of the TF-link being active using the 32 models
  indSingle = find(combConstr(:,k)==1);
  linkMargPosteriors{1}(:,k) = logsumexp(results_b.marlls(:, indSingle), 2) - ...
                          logsumexp(results_b.marlls(:, setdiff(ind32, indSingle)), 2);
        
  % posterior of the TF-link being active using the 32 models plus prior  
  linkMargPosteriors{2}(:,k) = logsumexp(withprior(:, indSingle), 2) - ...
                        logsumexp(withprior(:, setdiff(ind32, indSingle)), 2); 
                      
  % posterior of the TF-link being active using the 16 models   
  linkMargPosteriors{3}(:,k) = logsumexp(results_b.marlls(:, comb16(:, k)==1), 2) - ...
                          logsumexp(results_b.marlls(:, comb16(:, k)==0), 2);
               
  % posterior of the TF-link being active using the 16 models plus prior  
  linkMargPosteriors{4}(:,k) = logsumexp(withprior(:, comb16(:, k)==1), 2) - ...
                        logsumexp(withprior(:, comb16(:, k)==0), 2);
 
  % find the *single* index inside comb where only "k" is active 
  indIndiv = find(combConstr(:,k)==1 & sum(combConstr,2)==1 );     
  
  % posterior of the TF-link being active using the 2 models              
  linkMargPosteriors{5}(:,k) = results_b.marlls(:, indIndiv) - results_b.marlls(:, ind0);
  
  % posterior of the TF-link being active using the 2 models plus prior             
  linkMargPosteriors{6}(:,k) = withprior(:, indIndiv) - withprior(:, ind0);
  
  % posterior of the TF-link being active using maximum likelihood model             
  linkMargPosteriors{7}(:,k) = baseline_a.marlls(:, k+1) - baseline_a.marlls(:, 1);
  
  % posterior of the TF-link being active using inferelator           
  linkMargPosteriors{8}(:,k) = abs(singles(:,k));
%  
end


% 2 PLOT ---------------------------------------------------------------
% GLOBAL RANKING BASED ON PRESENCE OF SINGLE LINKS
% START  ---------------------------------------------------------------
%
% separate prior for eacn link being active 
%priorSingleTF = nanmean(M);
%prioraccsSingleTF(1) = mean(priorSingleTF);
%prioraccsSingleTF(2) = sum(priorSingleTF .* priorSingleTF) / sum(priorSingleTF).^2;
for k=1:numTFs
    indSingle = find(combConstr(:,k)==1);
    priorSingleTF(k) = sum(prior32(indSingle));
end

% Bayes correct classification rate based on a classifier that  uses
% a unifrom prior (as classification rule)
prioraccsSingleTF(1) = mean(priorSingleTF);
% Bayes correct classification rate based on a classifier that uses 
% the empirical prior
prioraccsSingleTF(2) = sum(priorSingleTF .* priorSingleTF) / sum(priorSingleTF);

%prioraccsSingleTF(1) = mean(priorSingleTF);
%prioraccsSingleTF(2) = sum(priorSingleTF .* priorSingleTF) / sum(priorSingleTF).^2;
baselines2 = [1, 2, 1, 2, 1. 2, 1, 1];
numGenes = size(linkMargPosteriors{1},1);
r2 = zeros(length(T2), length(linkMargPosteriors));
pvals2 = r2;
for k=1:length(linkMargPosteriors), 
    BestLink = [];
    BestPost = []; 
    % for loop over the number of genes 
    % that computes the best single-TF-link for each gene
    for n=1:size(linkMargPosteriors{k},1)
       [foo, BL] = sort(linkMargPosteriors{k}(n,:), 'descend');
       BestLink(n) = BL(1); 
       BestPost(n) = foo(1);
    end
    [foo, I] = sort(BestPost, 'descend');
    for l=1:length(T2),
       r2(l,k) = 0;  
       nanCnt =0;
       for j=1:T2(l)
          gene = I(j);   
          TF = BestLink(I(j));
          MM = M(gene, TF);
          if ~isnan(MM)
             r2(l,k) = r2(l,k) + MM;
          else
             nanCnt = nanCnt + 1;
          end
       end
       pvals2(l, k) = 1 - binocdf(r2(l,k) - 1, ...
			          T2(l)-nanCnt, ...
	  		          prioraccsSingleTF(baselines2(k)));
       r2(l,k) = r2(l,k)/(T2(l)-nanCnt); 
    end
    % OLD CODE THA WAS DOING SOMETHING NOT EXACTLY THE SAME
    %[foo, I] = sort(linkMargPosteriors{k}(:), 'descend');
    %for l=1:length(T2),
    %   r2(l,k) = 0;  
    %   nanCnt =0;
    %   for j=1:T2(l)
    %      gene = mod(I(j), numGenes);
    %      gene(gene==0)=numGenes;          
    %      TF = floor(I(j)/numGenes) + 1;
    %      TF(TF==6)=5; 
    %      MM = M(gene, TF);
    %      if ~isnan(MM)
    %         r2(l,k) = r2(l,k) + MM;
    %      else
    %         nanCnt = nanCnt + 1;
    %      end
    %   end
    %   %pvals3(l, k) = 1 - binocdf(nansum(M(I(1:T2(l)))) - 1, ...
	%   %		          sum(~isnan(M(I(1:T2(l))))), ...
	%   %		          prioraccs2(baselines2(k)));
    %   pvals2(l, k) = 1 - binocdf(r2(l,k) - 1, ...
	%		          T2(l)-nanCnt, ...
	%		          prioraccsSingleTF(baselines2(k)));
    %   r2(l,k) = r2(l,k)/(T2(l)-nanCnt); 
    %end
%
end
% plots bars 
h2 = figure;
set(gca, 'FontSize', fontSize);
h = bar(100*r2(:, plotRest==1));
set(gca, 'XTickLabel', T2);
hold on
v = axis;
v(3:4) = [0 100];
axis(v)
plot(v(1:2), 100*prioraccsSingleTF(1)*[1 1], 'b');
if incPrior,
  plot(v(1:2), 100*prioraccsSingleTF(2)*[1 1], 'g');
end
hold off
legends = {'Posterior-32', 'Posterior-32 + prior', 'Posterior-16', 'Posterior-16 + prior', ...
       'Posterior-2', 'Posterior-2 + prior',...
       'ML-Baseline', 'Inferelator', 'Uniform prior', 'Empirical prior'};
legend(legends([plotRest, 1, incPrior]==1));  
%legend('Posterior-32', 'Posterior-32 + prior', 'Posterior-16', 'Posterior-16 + prior', ...
%       'Posterior-2', 'Posterior-2 + prior',...
%       'Baseline', 'Uniform random', 'Random from prior', ...
%       'Location', 'EastOutside');
axis(v)
xlabel('# of top predictions')
ylabel('Enrichment (%)')
drosStarBars(h, pvals2(:, plotRest==1));
set(gca, 'FontSize', fontSize);
set(gcf, 'PaperUnits', 'centimeters')
set(gcf, 'PaperPosition', [0, 0, figSize])
% 2 PLOT ---------------------------------------------------------------
% GLOBAL RANKING BASED ON PRESENCE OF SINGLE LINKS
% END    ---------------------------------------------------------------


% Computation of marginal posterior probability over pairs of links 
linkPairPosteriors = {};
cnt = 0;
for k=1:numTFs
  for g=(k+1):numTFs
  %   
  cnt = cnt + 1;
  
  pairs(cnt,:) = [k g];
  % find all the indices inside comb where the TFs "k" and "g" are active 
  indPair = find(combConstr(:,k)==1 & combConstr(:,g)==1);
   
  % posterior probability of the TF-pair under the 32 hypotheses case
  linkPairPosteriors{1}(:, cnt) = logsumexp(results_b.marlls(:,  indPair), 2) - ...
			                      logsumexp(results_b.marlls(:, setdiff(ind32, indPair)), 2);
    
  % posterior probability of the TF-pair under the 32 hypotheses case
  linkPairPosteriors{2}(:, cnt) = logsumexp(withprior(:,  indPair), 2) - ...
			                      logsumexp(withprior(:, setdiff(ind32, indPair)), 2);
                                  
  % find the *single* index inside comb (restricted to 16 hypotheses) where 
  % the TFs "k" and "g" are active 
  indPairSingle = find(combConstr(:,k)==1 & combConstr(:,g)==1 & sum(combConstr,2)==2 );
 
  % posterior probability of the TF-pair under the 16 hypotheses case
  linkPairPosteriors{3}(:, cnt) = results_b.marlls(:,  indPairSingle) - ...
			              logsumexp(results_b.marlls(:, setdiff(ind16, indPairSingle)), 2);                   
                      
  % posterior probability of the TF-pair under the 16 hypotheses case
  linkPairPosteriors{4}(:, cnt) = withprior(:,  indPairSingle) - ...
			              logsumexp(withprior(:, setdiff(ind16, indPairSingle)), 2);
    
  % find the indices inside comb restricted to 4 models 
  %(0 0; k 0; 0 g; k g )
  ind4 = find(  sum(combConstr(:, setdiff(1:numTFs, [k g]) ), 2)==0  );
  
  % posterior probability of the TF-pair under the 4 hypotheses case
  linkPairPosteriors{5}(:, cnt) = results_b.marlls(:, indPairSingle) - ...
  			             logsumexp(results_b.marlls(:, setdiff(ind4, indPairSingle)), 2);
          
  % posterior probability of the TF-pair under the 4 hypotheses case
  linkPairPosteriors{6}(:, cnt) = withprior(:, indPairSingle) - ...
  			             logsumexp(withprior(:, setdiff(ind4, indPairSingle)), 2);
                     
  % baseline maximum likelihood model
  indPSinglebase = find(baselinecomb(:,k)==1 & baselinecomb(:,g)==1 & sum(baselinecomb,2)==2 ); 
  linkPairPosteriors{7}(:, cnt) = baseline_a.marlls(:, indPSinglebase) - baseline_a.marlls(:, 1);  
  
  %             
  end
end

% Inferelator predictions
linkPairPosteriors{8} = zeros(size(linkPairPosteriors{7}));
% twi & Mef2 
linkPairPosteriors{8}(:,9) = abs(doubles(:,1));
% bin & bap 
linkPairPosteriors{8}(:,6) = abs(doubles(:,2));
% tin & bap 
linkPairPosteriors{8}(:,3) = abs(doubles(:,3));
% bap & Mef2
linkPairPosteriors{8}(:,10) = abs(doubles(:,4));
% tin & Mef2
linkPairPosteriors{8}(:,4) = abs(doubles(:,5));

% compute also Antti's doubles 
doublesA = zeros(size(doubles)); 
for i=1:numGenes
   for k=1:10
       linkPairPosteriors{8}(i,k) = max(linkPairPosteriors{8}(i,k), min(abs(singles(i,pairs(k,1))),  abs(singles(i,pairs(k,2)))));
   end
end


% 3 PLOT ---------------------------------------------------------------
% GLOBAL RANKING BASED ON THE PRESENCE OF A PAIR OF  LINKS
% START  ---------------------------------------------------------------
%
% separate prior for each pairs of links being active 
cnt = 0;
for k=1:numTFs
  for g=(k+1):numTFs
      cnt = cnt + 1;
      indPair = find(combConstr(:,k)==1 & combConstr(:,g)==1);
      priorPairTF(cnt) = sum(prior32(indPair));
  end
end

% Bayes correct classification rate based on a classifier that  uses
% a unifrom prior (as classification rule)
prioraccsPairTF(1) = mean(priorPairTF);
% Bayes correct classification rate based on a classifier that uses 
% the empirical prior
prioraccsPairTF(2) = sum(priorPairTF.*priorPairTF) / sum(priorPairTF);


baselines2 = [1, 2, 1, 2, 1. 2, 1, 1];
r3 = zeros(length(T2), length(linkPairPosteriors));
pvals3 = r3;
for k=1:length(linkPairPosteriors),
    BestLink = [];
    BestPost = []; 
    % for loop over the number of genes 
    % that computes the besrt pair-TF-link for each gene
    for n=1:size(linkPairPosteriors{k},1)
       [foo, BL] = sort(linkPairPosteriors{k}(n,:), 'descend');
       BestLink(n) = BL(1); 
       BestPost(n) = foo(1);
    end
    [foo, I] = sort(BestPost, 'descend');
    for l=1:length(T2),
       r3(l,k) = 0;  
       nanCnt =0;
       for j=1:T2(l)
          gene = I(j);   
          TFpair = BestLink(I(j));
          MM =  prod(M(gene, pairs(TFpair,:)));
          if ~isnan(MM)
             r3(l,k) = r3(l,k) + MM;
          else
             nanCnt = nanCnt + 1;
          end
       end
       pvals3(l, k) = 1 - binocdf(r3(l,k) - 1, ...
			          T2(l)-nanCnt, ...
			          prioraccsPairTF(baselines2(1)));                  
       r3(l,k) = r3(l,k)/(T2(l)-nanCnt);
    end
    %[foo, I] = sort(linkPairPosteriors{k}(:), 'descend');
    %for l=1:length(T2),
    %   r3(l,k) = 0;  
    %   nanCnt = 0;
    %   for j=1:T2(l)
    %      gene = mod(I(j), numGenes);
    %      gene(gene==0)=numGenes;          
    %      TFpair = floor(I(j)/numGenes) + 1;
    %      TFpair(TFpair==11)=10;
    % 
    %      MM =  prod(M(gene, pairs(TFpair,:)));
    %      if ~isnan(MM)
    %         r3(l,k) = r3(l,k) + MM;
    %      else
    %         nanCnt = nanCnt + 1;
    %      end
    %   end
    %   pvals3(l, k) = 1 - binocdf(r3(l,k) - 1, ...
	%		          T2(l)-nanCnt, ...
	%		          prioraccsPairTF(baselines2(1)));                  
    %   r3(l,k) = r3(l,k)/(T2(l)-nanCnt);
    %end
%
end
% plots bars 
h3 = figure;
set(gca, 'FontSize', fontSize);
h = bar(100*r3(:, plotRest==1));
set(gca, 'XTickLabel', T2);
hold on
v = axis;
v(3:4) = [0 100];
axis(v)
plot(v(1:2), 100*prioraccsPairTF(1)*[1 1], 'b');
if incPrior,
  plot(v(1:2), 100*prioraccsPairTF(2)*[1 1], 'g');
end
hold off
legends = {'Posterior-32', 'Posterior-32 + prior', 'Posterior-16', 'Posterior-16 + prior', 'Posterior-4','Posterior-4 + prior',...
       'ML-Baseline', 'Inferelator', 'Uniform prior', 'Empirical prior'};
legend(legends([plotRest, 1, incPrior]==1));  
%legend('Posterior-32', 'Posterior-32 + prior', 'Posterior-16', 'Posterior-16 + prior', 'Posterior-4','Posterior-4 + prior',...
%       'Baseline', 'Uniform random', 'Random from prior', ...
%       'Location', 'EastOutside');
axis(v)
xlabel('# of top predictions')
ylabel('Enrichment (%)')
drosStarBars(h, pvals3(:, plotRest==1));
set(gca, 'FontSize', fontSize);
set(gcf, 'PaperUnits', 'centimeters')
set(gcf, 'PaperPosition', [0, 0, figSize])
% 3 PLOT ---------------------------------------------------------------
% GLOBAL RANKING BASED ON THE PRESENCE OF A PAIR OF  LINKS
% END    ---------------------------------------------------------------



% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %       BELOW  ARE PLOTS FOR PREDICTING ABSENCE OF LINKS                                 %                 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% % Computation of marginal posterior probabilities over single links of being
% % inactive
% for k=1:numTFs
% %
%   % posterior of the TF-link being *inactive* using the 32 models
%   indSingle = find(combConstr(:,k)==0);
%   linkNegativeMargPosteriors{1}(:,k) = logsumexp(results_b.marlls(:, indSingle), 2) - ...
%                           logsumexp(results_b.marlls(:, setdiff(ind32, indSingle)), 2);
        
%   % posterior of the TF-link being *inactive* using the 32 models plus prior  
%   linkNegativeMargPosteriors{2}(:,k) = logsumexp(withprior(:, indSingle), 2) - ...
%                         logsumexp(withprior(:, setdiff(ind32, indSingle)), 2); 
                      
%   % posterior of the TF-link being *inactive* using the 16 models   
%   linkNegativeMargPosteriors{3}(:,k) = logsumexp(results_b.marlls(:, comb16(:, k)==0), 2) - ...
%                           logsumexp(results_b.marlls(:, comb16(:, k)==1), 2);
               
%   % posterior of the TF-link being *inactive* using the 16 models plus prior  
%   linkNegativeMargPosteriors{4}(:,k) = logsumexp(withprior(:, comb16(:, k)==0), 2) - ...
%                         logsumexp(withprior(:, comb16(:, k)==1), 2);
 
%   % find the *single* index inside comb where only "k" is *active* 
%   indIndiv = find(combConstr(:,k)==1 & sum(combConstr,2)==1 );     
  
%   % posterior of the TF-link being *inactive* using the 2 models              
%   linkNegativeMargPosteriors{5}(:,k) = results_b.marlls(:, ind0) - results_b.marlls(:, indIndiv);
  
%   % posterior of the TF-link being *inactive* using the 2 models plus prior             
%   linkNegativeMargPosteriors{6}(:,k) = withprior(:, ind0)  - withprior(:, indIndiv);
  
%   % posterior of the TF-link being *inactive* using maximum likelihood model             
%   linkNegativeMargPosteriors{7}(:,k) = baseline_a.marlls(:, 1) - baseline_a.marlls(:, k+1);
% %  
% end


% % 4 PLOT ---------------------------------------------------------------
% % GLOBAL RANKING BASED ON THE ABSENCE OF A SINGLE LINK
% % STARTS    ------------------------------------------------------------
% %
% % separate prior for eacn link being active 
% %priorSingleAbsentTF = 1 - priorSingleTF;
% %prioraccsSingleAbsentTF(1) = mean(priorSingleAbsentTF);
% %prioraccsSingleAbsentTF(2) = sum(priorSingleAbsentTF .* priorSingleAbsentTF) / sum(priorSingleAbsentTF).^2;
% for k=1:numTFs
%     indSingle = find(combConstr(:,k)==0);
%     priorSingleAbsentTF(k) = sum(prior32(indSingle));
% end

% % Bayes correct classification rate based on a classifier that  uses
% % a unifrom prior (as classification rule)
% prioraccsSingleAbsentTF(1) = mean(priorSingleAbsentTF);
% % Bayes correct classification rate based on a classifier that uses 
% % the empirical prior
% prioraccsSingleAbsentTF(2) = sum(priorSingleAbsentTF.*priorSingleAbsentTF) / sum(priorSingleAbsentTF);


% r4 = zeros(length(T2), length(linkNegativeMargPosteriors));
% pvals4 = r4;
% for k=1:length(linkNegativeMargPosteriors), 
%     BestLink = [];
%     BestPost = []; 
%     % for loop over the number of genes 
%     % that computes the best single-TF-absent-link for each gene
%     for n=1:size(linkNegativeMargPosteriors{k},1)
%        [foo, BL] = sort(linkNegativeMargPosteriors{k}(n,:), 'descend');
%        BestLink(n) = BL(1); 
%        BestPost(n) = foo(1);
%     end
%     [foo, I] = sort(BestPost, 'descend');
%     for l=1:length(T2),
%        r4(l,k) = 0;  
%        nanCnt =0;
%        for j=1:T2(l)
%           gene = I(j);   
%           TF = BestLink(I(j));
%           MM = M(gene, TF);
%           if ~isnan(MM)
%              r4(l,k) = r4(l,k) + 1-MM;
%           else
%              nanCnt = nanCnt + 1;
%           end
%        end
%        pvals4(l, k) = 1 - binocdf(r4(l,k) - 1, ...
% 			          T2(l)-nanCnt, ...
% 			          priorSingleAbsentTF (baselines2(k)));
%        r4(l,k) = r4(l,k)/(T2(l)-nanCnt); 
%     end
%     %[foo, I] = sort(linkNegativeMargPosteriors{k}(:), 'descend');
%     %for l=1:length(T2),
%     %   r4(l,k) = 0;  
%     %   nanCnt =0;
%     %   for j=1:T2(l)
%     %      gene = mod(I(j), numGenes);
%     %      gene(gene==0)=numGenes;          
%     %      TF = floor(I(j)/numGenes) + 1;
%     %      TF(TF==6)=5; 
%     %      MM = M(gene, TF);
%     %      if ~isnan(MM)
%     %         r4(l,k) = r4(l,k) + 1-MM;
%     %      else
%     %         nanCnt = nanCnt + 1;
%     %      end
%     %   end
%     %   %pvals3(l, k) = 1 - binocdf(nansum(M(I(1:T2(l)))) - 1, ...
% 	%   %		          sum(~isnan(M(I(1:T2(l))))), ...
% 	%   %		          prioraccs2(baselines2(k)));
%     %   pvals4(l, k) = 1 - binocdf(r4(l,k) - 1, ...
% 	%		          T2(l)-nanCnt, ...
% 	%		          priorSingleAbsentTF (baselines2(k)));
%     %   r4(l,k) = r4(l,k)/(T2(l)-nanCnt); 
%     %end
% %
% end
% % plots bars 
% h4 = figure;
% set(gca, 'FontSize', fontSize);
% h = bar(100*r4(:, plotRest==1));
% set(gca, 'XTickLabel', T2);
% hold on
% v = axis;
% v(3:4) = [0 100];
% axis(v)
% plot(v(1:2), 100*prioraccsSingleAbsentTF(1)*[1 1], 'b');
% if incPrior,
%   plot(v(1:2), 100*prioraccsSingleAbsentTF(2)*[1 1], 'g');
% end
% hold off
% legends = {'Posterior-32', 'Posterior-32 + prior', 'Posterior-16', 'Posterior-16 + prior', 'Posterior-2', 'Posterior-2 + prior',...
%        'Baseline', 'Uniform prior', 'Empirical prior', ...
%        'Location', 'EastOutside'};
% legend(legends([plotRest, 1, incPrior]==1));  
% axis(v)
% xlabel('# of top predictions')
% ylabel('Enrichment (%)')
% drosStarBars(h, pvals4(:, plotRest==1));
% set(gca, 'FontSize', fontSize);
% set(gcf, 'PaperUnits', 'centimeters')
% set(gcf, 'PaperPosition', [0, 0, figSize])
% % 4 PLOT ---------------------------------------------------------------
% % GLOBAL RANKING BASED ON THE ABSENCE OF A SINGLE LINK
% % END    ---------------------------------------------------------------


% % Computation of marginal posterior probability over pairs of links 
% % being both inactive (non-regulators)
% linkNegativePairPosteriors = {};
% cnt = 0;
% for k=1:numTFs
%   for g=(k+1):numTFs
%   %   
%   cnt = cnt + 1;
  
%   pairs(cnt,:) = [k g];
%   % find all the indices inside comb where the TFs "k" and "g" are *inactive* 
%   indPair = find(combConstr(:,k)==0 & combConstr(:,g)==0);
   
%   % posterior probability of the TF-pair being *inactive* under the 32 hypotheses case
%   linkNegativePairPosteriors{1}(:, cnt) = logsumexp(results_b.marlls(:,  indPair), 2) - ...
% 			                      logsumexp(results_b.marlls(:, setdiff(ind32, indPair)), 2);
    
%   % posterior probability of the TF-pair being *inactive*  under the 32 hypotheses case
%   linkNegativePairPosteriors{2}(:, cnt) = logsumexp(withprior(:,  indPair), 2) - ...
% 			                      logsumexp(withprior(:, setdiff(ind32, indPair)), 2);
                                  
%   % find the indices inside comb (restricted to 16 hypotheses) where 
%   % the TFs "k" and "g" are *inactive* 
%   indPairSingle = find(combConstr(:,k)==0 & combConstr(:,g)==0 & sum(combConstr,2)<=2 );
 
%   % posterior probability of the TF-pair  being *inactive*  under the 16 hypotheses case
%   linkNegativePairPosteriors{3}(:, cnt) = logsumexp(results_b.marlls(:,  indPairSingle), 2) - ...
% 			              logsumexp(results_b.marlls(:, setdiff(ind16, indPairSingle)), 2);                   
                      
%   % posterior probability of the TF-pair  being *inactive*  under the 16 hypotheses case
%   linkNegativePairPosteriors{4}(:, cnt) = logsumexp(withprior(:,  indPairSingle), 2) - ...
% 			              logsumexp(withprior(:, setdiff(ind16, indPairSingle)), 2);
  
%   % find the indices inside comb restricted to 4 models 
%   %(0 0; k 0; 0 g; k g )
%   ind4 = find(  sum(combConstr(:, setdiff(1:numTFs, [k g]) ), 2)==0  );
  
%   % posterior probability of the TF-pair being *inactive* under the 4 hypotheses case
%   linkNegativePairPosteriors{5}(:, cnt) = results_b.marlls(:, ind0) - ...
%   			             logsumexp(results_b.marlls(:, setdiff(ind4, ind0)), 2);
          
%   % posterior probability of the TF-pair being *inactive* under the 4 hypotheses case
%   linkNegativePairPosteriors{6}(:, cnt) = withprior(:, ind0) - ...
%   			             logsumexp(withprior(:, setdiff(ind4, ind0)), 2);
                     
%   % baseline maximum likelihood model [[k=0.g=0, rest=1] minus the zero model)
%   indPSinglebase = find(baselinecomb(:,k)==0 & baselinecomb(:,g)==0 & sum(baselinecomb,2)==1 ); 
%   linkNegativePairPosteriors{7}(:, cnt) = logsumexp(baseline_a.marlls(:, indPSinglebase), 2) - baseline_a.marlls(:, 1);  
%   %             
%   end
% end



% % 5 PLOT ---------------------------------------------------------------
% % GLOBAL RANKING BASED ON THE ABSENCE OF A PAIR OF LINKS
% % STARTS ---------------------------------------------------------------
% cnt = 0;
% for k=1:numTFs
%   for g=(k+1):numTFs
%       cnt = cnt + 1;
%       indAbsentPair = find(combConstr(:,k)==0 & combConstr(:,g)==0);
%       priorPairAbsentTF(cnt) = sum(prior32(indAbsentPair));
%   end
% end

% % Bayes correct classification rate based on a classifier that  uses
% % a unifrom prior (as classification rule)
% prioraccsPairAbsentTF(1) = mean(priorPairAbsentTF);
% % Bayes correct classification rate based on a classifier that uses 
% % the empirical prior
% prioraccsPairAbsentTF(2) = sum(priorPairAbsentTF.*priorPairAbsentTF) / sum(priorPairAbsentTF);

% % Antti's way (not clear waht Bayes correct rate is that)
% %prioraccsPairAbsentTF(1) = mean(priorPairAbsentTF);
% %prioraccsPairAbsentTF(2) = sum(priorPairAbsentTF .* priorPairAbsentTF) / sum(priorPairAbsentTF).^2;

% r5 = zeros(length(T2), length(linkNegativePairPosteriors));
% pvals5 = r5;
% for k=1:length(linkNegativePairPosteriors),
%     BestLink = [];
%     BestPost = []; 
%     % for loop over the number of genes 
%     % that computes the besrt pair-TF-absent-link for each gene
%     for n=1:size(linkNegativePairPosteriors{k},1)
%        [foo, BL] = sort(linkNegativePairPosteriors{k}(n,:), 'descend');
%        BestLink(n) = BL(1); 
%        BestPost(n) = foo(1);
%     end
%     [foo, I] = sort(BestPost, 'descend');
%     for l=1:length(T2),
%        r5(l,k) = 0;  
%        nanCnt =0;
%        for j=1:T2(l)
%           gene = I(j);   
%           TFpair = BestLink(I(j)); 
%           MM =  prod(1 - M(gene, pairs(TFpair,:)));
%           if ~isnan(MM)
%              r5(l,k) = r5(l,k) + MM;
%           else
%              nanCnt = nanCnt + 1;
%           end          
%        end
%        pvals5(l, k) = 1 - binocdf(r5(l,k) - 1, ...
% 			          T2(l)-nanCnt, ...
% 			          prioraccsPairAbsentTF(baselines2(1)));                  
%        r5(l,k) = r5(l,k)/(T2(l)-nanCnt);
%     end
%     %[foo, I] = sort(linkNegativePairPosteriors{k}(:), 'descend');
%     %for l=1:length(T2),
%     %   r5(l,k) = 0;  
%     %   nanCnt = 0;
%     %   for j=1:T2(l)
%     %      gene = mod(I(j), numGenes);
%     %      gene(gene==0)=numGenes;          
%     %      TFpair = floor(I(j)/numGenes) + 1;
%     %      TFpair(TFpair==11)=10;
%     % 
%     %      MM =  prod(1 - M(gene, pairs(TFpair,:)));
%     %      if ~isnan(MM)
%     %         r5(l,k) = r5(l,k) + MM;
%     %      else
%     %         nanCnt = nanCnt + 1;
%     %      end
%     %   end
%     %   pvals5(l, k) = 1 - binocdf(r5(l,k) - 1, ...
% 	%		          T2(l)-nanCnt, ...
%     %		        prioraccsPairAbsentTF(baselines2(1)));                  
%     %  r5(l,k) = r5(l,k)/(T2(l)-nanCnt);
%     %end
% %
% end
% % plots bars 
% h5 = figure;
% set(gca, 'FontSize', fontSize);
% h = bar(100*r5(:, plotRest==1));
% set(gca, 'XTickLabel', T2);
% hold on
% v = axis;
% v(3:4) = [0 100];
% axis(v)
% plot(v(1:2), 100*prioraccsPairAbsentTF(1)*[1 1], 'b');
% if incPrior,
%   plot(v(1:2), 100*prioraccsPairAbsentTF(2)*[1 1], 'g');
% end
% hold off
% legends = {'Posterior-32', 'Posterior-32 + prior', 'Posterior-16', 'Posterior-16 + prior','Posterior-4', 'Posterior-4 + prior',...
%        'Baseline', 'Uniform prior', 'Empirical prior', ...
%        'Location', 'EastOutside'};
% legend(legends([plotRest, 1, incPrior]==1));  
% axis(v)
% xlabel('# of top predictions')
% ylabel('Enrichment (%)')
% drosStarBars(h, pvals5(:, plotRest==1));
% set(gca, 'FontSize', fontSize);
% set(gcf, 'PaperUnits', 'centimeters')
% set(gcf, 'PaperPosition', [0, 0, figSize])
% % 5 PLOT ---------------------------------------------------------------
% % GLOBAL RANKING BASED ON THE ABSENCE OF A PAIR OF LINKS
% % END    ---------------------------------------------------------------



% print Plots
property = 'Constrained';
if flag ~= 1
    property = 'Unconstrained'; 
end
if printPlot 
   print(h1, '-depsc2', [ddir 'drosophilaBars_' 'EnrichmentGlobalMAP' num2str(2) '_', property '.eps']);
   print(h1b, '-depsc2', [ddir 'drosophilaBars_' 'EnrichmentGlobalMAPIN_SITU' num2str(2) '_', property '.eps']);
   print(h2, '-depsc2', [ddir 'drosophilaBars_' 'EnrichmentSingleLinks' num2str(sum(plotRest)) '_', property '.eps']); 
   print(h3, '-depsc2', [ddir 'drosophilaBars_' 'EnrichmentPairLinks' num2str(sum(plotRest)) '_', property '.eps']);
   %print(h4, '-depsc2', [ddir 'drosophilaBars_' 'EnrichmentNegativeSingleLinks' num2str(sum(plotRest)) '_', property '.eps']);
   %print(h5, '-depsc2', [ddir 'drosophilaBars_' 'EnrichmentNegativePairLinks' num2str(sum(plotRest)) '_', property '.eps']);
end
