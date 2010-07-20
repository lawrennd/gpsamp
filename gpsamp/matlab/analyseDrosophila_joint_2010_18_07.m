% USER-specified: Do you want constrained (positive interactions weights)
% or unconstrained models? 
% !! Warning: Baseline models was run only for the constrained case !!
flag = 1; % "1" for  constrained; "anything else" for unconstrained 

% USER-specified: Sizes of the ranking sets for the first plot
T1 = [200, 400, 800, 1600, 3200, 6000];
% Sizes of the ranking sets for the second plot
T2 = [100, 200, 400, 800, 1600, 3200];

% USER-specified:: Directory whwre you store figures
ddir = 'figures/';
printPlot = 0; % 0 means not printing

% models (as stored  in the results variables; see below) 
% correspidng to 5 TFs being active/inactive 
combConstr = [0 0 0 0 0; 0 0 1 0 0; 0 0 0 0 1; 0 0 1 0 1;
	1 0 0 0 0; 0 1 0 0 0; 0 0 0 1 0;
	1 1 0 0 0; 1 0 1 0 0; 1 0 0 1 0; 1 0 0 0 1;
	0 1 1 0 0; 0 1 0 1 0; 0 1 0 0 1;
	0 0 1 1 0;
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
    % load resutls for 4 models using at most pairs of twist and mef2
    % (the decay-zeroth (first) model had an error... )
    results_b = sortResults(load('results/multitf5a_2010-04-22_summary.mat'));
    % load the correct result for the zeroth model
    results_b2 = sortResults(load('results/multitf5a_null_2010-06-07_summary.mat'));
    % load the 12 models that correspond to a most of 2 TFs
    results_b3 = sortResults(load('results/multitf6a_2010-06-24_summary.mat'));
    % correct the zeroth model
    results_b.marlls(:, 1) = results_b2.marlls;
    % append the 4-models with the remaining 12 models
    results_b.marlls = [results_b.marlls, results_b3.marlls];
    % append with the remaining 16 models to get the total 32 models
    results_b4 = sortResults(load('results/multitf6b_2010-07-13_summary.mat'));
    results_b.marlls = [results_b.marlls, results_b4.marlls];

    % load the 15 baseline models (zerhoth model is excluded)
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

% number of TFs
numTFs = size(combConstr,2);
M = drosMakeValidationMatrix(chipdistances, results_b.genes, 2000);


% compuation of the prior over the network
prior = zeros(1, size(combConstr, 1));
for k=1:size(combConstr, 1),
  prior(k) = mean(all(M(:, combConstr(k, :)==1), 2));
end
prioracc = prior;
prior = prior / sum(prior);
% end of compuating the prior 


% the WEIRD prioracc computation (I have no clue what Antti is doing here)
prioraccs(1) = mean(prioracc(2:end));
prioraccs(2) = sum(prioracc(2:end) .* prior(2:end)) / sum(prior(2:end));


% non-uniform prior over models 
withprior = results_b.marlls + repmat(log(prior), [size(results_b.marlls, 1), 1]);

% indices of all 32 models 
ind32 = 1:size(results_b.marlls, 2); 
% indices that correspond to 16 models (at most 2 TFs) 
ind16 = find(sum(combConstr,2)<=2); 
% % indices that correspond to 6 models (at most 1 TF) -> single TFs models  
%ind6 = find(sum(combConstr,2)<=1);
% index of the zeroth model 
ind0 = find(sum(combConstr,2)==0);
linkMargPosteriors = {};
comb16 = combConstr(ind16,:);
%comb6 = combConstr(ind6,:);

% the next two computations are numerically unstable  
%posteriors{1} = results_b.marlls - repmat(logsumexp(results_b.marlls, 2), [1, size(results_b.marlls, 2)]);
%posteriors{2} = withprior - repmat(logsumexp(withprior, 2), [1, size(withprior, 2)]);

posteriors = {};
for k=1:length(ind32)
    % posterior of each model in the 32-model 
    posteriors{1}(:,k) = results_b.marlls(:,k) - logsumexp(results_b.marlls(:, setdiff(ind32, k)), 2);
    % posterior of each model in the 32-model case plus prior 
    posteriors{2}(:,k) = withprior(:,k) - logsumexp(withprior(:, setdiff(ind32, k)), 2);
end

prior16 = zeros(1, size(comb16, 1));
for k=1:size(comb16, 1),
  prior16(k) = mean(all(M(:, comb16(k, :)==1), 2));
end
prior16 = prior16 / sum(prior16);
% non-uniform prior over models 
withprior16 = results_b.marlls(:,ind16) + repmat(log(prior16), [size(results_b.marlls(:,ind16), 1), 1]);
for k=1:length(ind16)
    % posterior of each model in the 32-model 
    posteriors{3}(:,k) = results_b.marlls(:,ind16(k)) - logsumexp(results_b.marlls(:, setdiff(ind16, k)), 2);
    % posterior of each model in the 32-model case plus prior 
    posteriors{4}(:,k) = withprior16(:,k) - logsumexp(withprior16(:, [1:k-1,k+1:length(ind16)]), 2);
end


%prior6 = zeros(1, size(comb6, 1));
%for k=1:size(comb6, 1),
%  prior6(k) = mean(all(M(:, comb6(k, :)==1), 2));
%end
%prior6 = prior6 / sum(prior6);
% % non-uniform prior over models 
%withprior6 = results_b.marlls(:,ind6) + repmat(log(prior6), [size(results_b.marlls(:,ind6), 1), 1]);
%for k=1:length(ind6)
%    % posterior of each model in the 32-model 
%    posteriors{5}(:,k) = results_b.marlls(:,ind6(k)) - logsumexp(results_b.marlls(:, setdiff(ind6, k)), 2);
%    % posterior of each model in the 32-model case plus prior 
%    posteriors{6}(:,k) = withprior6(:,k) - logsumexp(withprior6(:, [1:k-1,k+1:length(ind6)]), 2);
%end


% baseline models likelihood ratios 
prioraccinds = [1, 2, 1, 2, 1, 2, 1];
posteriors{end+1} = baseline_a.marlls - repmat(baseline_a.marlls(:, 1), ...
					   [1, size(baseline_a.marlls, 2)]);

scores = {}; I={}; sscores={}; J={};
r = zeros(length(posteriors), length(T1));
r2 = zeros(length(posteriors), length(T1));
pvals = zeros(length(posteriors), length(T1));
for k=1:length(posteriors)   
  if k <=2, % 32 models 
    mycomb = combConstr;
  elseif  (k>2 & k<=4), % 16 models 
    mycomb = comb16;
  %elseif  (k>4 & k<=6), % 6 models 
  %  mycomb = comb6;
  else
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
          acc2 = acc2 + all(M(J{k}(m), :) == mycomb(I{k}(J{k}(m)), :), 2);
	      count = count + 1;
      end
    end
    r(k, l) = acc / count;
    r2(k, l) = acc2 / count;
    pvals(k, l) = 1 - binocdf(acc - 1, count, prioraccs(prioraccinds(k)));
  end
end

h1 = figure;
% plot bars
h = bar(100*r');

set(gca, 'XTickLabel', T1);
hold on
v = axis;
v(3:4) = [0 50];
axis(v)
plot(v(1:2), 100*prioraccs(1)*[1 1], 'b');
plot(v(1:2), 100*prioraccs(2)*[1 1], 'g');
hold off
legend('MAP-32 + uniform prior', 'MAP-32 + prior', 'MAP-16 + uniform prior', 'MAP-16 + prior', 'Baseline', 'Uniform random', 'Random from prior', 'Location', 'EastOutside')
axis(v)
xlabel('# of top genes')
ylabel('Enrichment (%)')
drosStarBars(h, pvals');

for k=1:numTFs
% 
  indSingle = find(combConstr(:,k)==1);
  % posterior of the TF-link being active using the 32 models   
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
%  
end

prior2 = nanmean(M);
prioraccs2(1) = mean(prior2);
prioraccs2(2) = sum(prior2 .* prior2) / sum(prior2).^2;

baselines2 = [1, 2, 1, 2, 1. 2, 1];

r3 = zeros(length(T2), length(linkMargPosteriors));
pvals3 = r3;
for k=1:length(linkMargPosteriors),
  [foo, I] = sort(linkMargPosteriors{k}(:), 'descend');
  for l=1:length(T2),
    r3(l, k) = nanmean(M(I(1:T2(l))));
    pvals3(l, k) = 1 - binocdf(nansum(M(I(1:T2(l)))) - 1, ...
			       sum(~isnan(M(I(1:T2(l))))), ...
			       prioraccs2(baselines2(k)));
  end
end


% plots bars 
h2 = figure;
h = bar(100*r3);
set(gca, 'XTickLabel', T2);
hold on
v = axis;
v(3:4) = [0 60];
axis(v)
plot(v(1:2), 100*prioraccs2(1)*[1 1], 'b');
plot(v(1:2), 100*prioraccs2(2)*[1 1], 'g');
hold off
legend('Posterior-32', 'Posterior-32 + prior', 'Posterior-16', 'Posterior-16 + prior', ...
       'Posterior-2', 'Posterior-2 + prior',...
       'Baseline', 'Uniform random', 'Random from prior', ...
       'Location', 'EastOutside');
axis(v)
xlabel('# of top predictions')
ylabel('Enrichment (%)')

drosStarBars(h, pvals3);

% print Plots
property = 'Constrained';
if flag ~= 1
    property = 'Unconstrained'; 
end
dd = date;
if printPlot 
   print(h1, '-depsc2', [ddir 'drosophilaBars_' 'EnrichmentGlobalMAP_', property dd '.eps']);
   print(h2, '-depsc2', [ddir 'drosophilaBars_' 'EnrichmentMargLinks_', property dd '.eps']); 
end