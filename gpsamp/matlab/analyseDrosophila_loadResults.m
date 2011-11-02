function [posteriors, linkMargPosteriors, linkPairPosteriors, priors, M, post_genes] = analyseDrosophila_loadResults()
% USER-specified: Do you want constrained (positive interactions weights)
% or unconstrained models? 
% !! Warning: Baseline models was run only for the constrained case !!
flag = 1; % "1" for  constrained; "anything else" for unconstrained 

load datasets/drosophila_data;
load datasets/testset;

analyseDrosophila_constants;

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
regression_a.genes = drosexp.genes(testset.indices); 

% load inferelator results
load inferelator_models.mat; 

inf_genes = cell(size(all_rows));
for k=1:length(all_rows),
  inf_genes{k} = drosexp.genes{find(strcmp(all_rows(k), drosexp.probes))};
end

% Permute the rows to match other results
[dummy, I] = sort(inf_genes);
inf_genes = inf_genes(I);
singles = singles(I, :);
doubles = doubles(I, :);
all_rows = all_rows(I);
regression_a.marlls = regression_a.marlls(I, :);
regression_a.genes = regression_a.genes(I);

[pair_a, pair_b] = strtok(doubles_cols, '.');
inf_pairs = zeros(length(pair_a), 2);
for k=1:length(pair_b),
  pair_b{k} = pair_b{k}(2:end);
  inf_pairs(k,:) = [find(strcmp(pair_a{k}, singles_cols)),
		    find(strcmp(pair_b{k}, singles_cols))];
end

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
singles = singles(mask==1,:);
doubles = doubles(mask==1,:);
all_rows = all_rows(mask==1);
inf_genes = inf_genes(mask==1);

% number of TFs
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
all_rows = all_rows(mask==1);
inf_genes = inf_genes(mask==1);
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
priors.prior32 = prior32;
priors.accinds = [1, 2, 1, 2, 1, 1, 1];
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

    
linkMargPosteriors = {};

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


% baseline models likelihood ratios 
posteriors{end+1} = baseline_a.marlls - repmat(baseline_a.marlls(:, 1), ...
					   [1, size(baseline_a.marlls, 2)]);

posteriors{end+1} = regression_a.marlls;

          
% compute global ranking perforamcne using random prediction. Compute

accRand = 0;
countRand = 0;
accRand2 = 0;
countRand2 = 0;
accPrior = 0;
countPrior = 0;
accPrior2 = 0;
countPrior2 = 0;
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
    accRand2 = accRand2 + all( M(perm(1),:) == Mrand );
    countRand2 = countRand2 + 1;
    if sum(Mprior(1,:),2) > 0 
	   accPrior = accPrior + all( M(perm(1), Mprior==1), 2);
	   countPrior = countPrior + 1;
    end
    accPrior2 = accPrior2 + all( M(perm(1),:) == Mprior );
    countPrior2 = countPrior2 + 1;
end
priors.accs31(1) =  accRand/countRand;
priors.accs31(2) =  accPrior/countPrior;
priors.accs31_neg(1) =  accRand2/countRand2;
priors.accs31_neg(2) =  accPrior2/countPrior2;

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
priors.focused_accs31(1) =  accRand/countRand;
priors.focused_accs31(2) =  accPrior/countPrior;


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
  
  %inf_pairs =  [3 5; 2 4; 2 3; 4 5; 3 4; 1 4; 1 5];
  % posterior of the TF-link being active using inferelator           
  %linkMargPosteriors{8}(:,k) = abs(singles(:,k));
  linkMargPosteriors{8}(:,k) = max(abs([singles(:,k), doubles(:, any(inf_pairs == k, 2))]), [], 2);
%  
end



% Computation of marginal posterior probability over pairs of links 
linkPairPosteriors = {};
cnt = 0;
for k=1:numTFs
  for g=(k+1):numTFs
  %   
  cnt = cnt + 1;
  
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
  
  % posterior probability of the TF-pair by combining *Heuristically*
  % the single TFs models 
  indk = find(combConstr(:,k)==1  &  sum(combConstr,2)==1);
  indg = find(combConstr(:,g)==1  &  sum(combConstr,2)==1);
  
  linkPairPosteriors{9}(:, cnt) = ...
      results_b.marlls(:, indk) + results_b.marlls(:, indg)  - ...
      results_b.marlls(:, ind0) - logsumexp(results_b.marlls(:, [ind0  indk indg]),2);
  %             
  end
end

% Inferelator predictions
linkPairPosteriors{8} = zeros(size(linkPairPosteriors{7}));

for k=1:size(doubles, 2),
  linkPairPosteriors{8}(:, all(bsxfun(@eq, inf_pairs(k,:), pairs), 2)) = abs(doubles(:,k));
end

% compute also Antti's doubles 
for i=1:numGenes
   for k=1:10
       linkPairPosteriors{8}(i,k) = max(linkPairPosteriors{8}(i,k), min(abs(singles(i,pairs(k,1))),  abs(singles(i,pairs(k,2)))));
   end
end




post_genes = results_b.genes;
