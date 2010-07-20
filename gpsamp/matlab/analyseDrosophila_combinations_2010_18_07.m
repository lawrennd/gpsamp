
% USER-specified: Do you want constrained (positive interactions weights)
% or unconstrained models? 
% !! Warning: Baseline models was run only for the constrained case !!
flag = 1; % "1" for  constrained; "anything else" for unconstrained 

% USER-specified: Sizes of the ranking sets 
T = [20, 50, 100, 150, 200];

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
    % load results for 4 models using at most pairs of twist and mef2
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


numTFs = size(combConstr,2);
M = drosMakeValidationMatrix(chipdistances, results_b.genes, 2000);

% indices of all 32 models 
ind32 = 1:size(results_b.marlls, 2); 
% indices that correspond to 16 models (at most 2 TFs) 
ind16 = find(sum(combConstr,2)<=2);
% index of the zeroth model 
ind0 = find(sum(combConstr,2)==0);

% computation of pair-probabilities using set of 32 hypotheses (models)
J_joint32 = {};
% computation of pair-probabilities using set of 16 hypotheses (models)
J_joint16 = {};
% computation of pair-probabilities using set of 4 hypotheses (models)
J_joint4 = {};
% computation of pair-probabilities using set of 2 hypotheses (zeroth
% model and the 2-TF model)
J_joint2 = {};
cnt = 0;
%
for k=1:numTFs
  for g=(k+1):numTFs
  %   
  cnt = cnt + 1;
  
  % find all the indices inside comb where the TFs "k" and "g" are active 
  indPair = find(combConstr(:,k)==1 & combConstr(:,g)==1);
   
  % posterior probability of the TF-pair under the 32 hypotheses case
  [foo, J_joint32{cnt}] = sort(logsumexp(results_b.marlls(:,  indPair), 2) - ...
			              logsumexp(results_b.marlls(:, setdiff(ind32, indPair)), 2), 'descend');
    
  % find the *single* index inside comb (restricted to 16 hypotheses) where 
  % the TFs "k" and "g" are active 
  indPairSingle = find(combConstr(:,k)==1 & combConstr(:,g)==1 & sum(combConstr,2)==2 );
 
  % posterior probability of the TF-pair under the 16 hypotheses case
  [foo, J_joint16{cnt}] = sort(results_b.marlls(:,  indPairSingle) - ...
			              logsumexp(results_b.marlls(:, setdiff(ind16, indPairSingle)), 2), 'descend');
  
  % find the indices inside comb restricted to 4 models 
  %(0 0; k 0; 0 g; k g )
  ind4 = find(  sum(combConstr(:, setdiff(1:numTFs, [k g]) ), 2)==0  );
  
  % posterior probability of the TF-pair under the 4 hypotheses case
  [foo, J_joint4{cnt}] = sort(results_b.marlls(:, indPairSingle) - ...
  			             logsumexp(results_b.marlls(:, setdiff(ind4, indPairSingle)), 2), 'descend');
          
  % posterior probability of the TF-pair under the 2 hypotheses case (NOT PLOTTED)
  [foo, J_joint2{cnt}] = sort(results_b.marlls(:, indPairSingle) - ...
  			             results_b.marlls(:, ind0), 'descend');        

  % baseline maximum likelihood model
  indPSinglebase = find(baselinecomb(:,k)==1 & baselinecomb(:,g)==1 & sum(baselinecomb,2)==2 ); 
  [foo, J_jointbase{cnt}] = sort(baseline_a.marlls(:, indPSinglebase) - baseline_a.marlls(:, 1),'descend');  
  %             
  end
end

J_indiv32 = {};
J_indiv16 = {};
J_indiv2 = {};
J_indbase = {};
comb16 = combConstr(ind16,:);
for k=1:numTFs,
%
  % single TF probability computed from 32 models 
  
  indSingle = find(combConstr(:,k)==1);
  
  [foo, J_indiv32{k}] = sort(logsumexp(results_b.marlls(:,indSingle),2) - ...
			            logsumexp(results_b.marlls(:,setdiff(ind32, indSingle)),2),'descend');
    
  % single TF probability computed from 16 models                
  [foo, J_indiv16{k}] = sort(logsumexp(results_b.marlls(:,comb16(:, k)==1),2) - ...
			            logsumexp(results_b.marlls(:,comb16(:, k)==0),2),'descend');
      
  % find the *single* index inside comb where only "k" is active 
  indIndiv = find(combConstr(:,k)==1 & sum(combConstr,2)==1 );     
  
  % single TF probability computed from 2 models             
  [foo, J_indiv2{k}] = sort(results_b.marlls(:, indIndiv) - results_b.marlls(:, ind0),'descend');
  
  % baseline maximum likelihood evaluation         
  [foo, J_indbase{k}] = sort(baseline_a.marlls(:, k+1) - baseline_a.marlls(:, 1), 'descend');
 
%  
end


h1 = figure; 
cnt = 0;
for k=1:numTFs
  for g=(k+1):numTFs
  cnt = cnt + 1;    
  figure(h1);
  subplot(2, 5, cnt);
  %drosPlotAccuracyBars({J_joint32{cnt}, J_joint16{cnt}, J_joint4{cnt},  J_joint2{cnt}, J_jointbase{cnt}}, prod(M(:, [k g]), 2), T);
  drosPlotAccuracyBars({J_joint32{cnt}, J_joint16{cnt}, J_joint4{cnt},  J_jointbase{cnt}}, prod(M(:, [k g]), 2), T);
  title(sprintf('%s + %s', drosTF.names{k}, drosTF.names{g}));
  end
end
%legend('Posterior-32', 'Posterior-16', 'Posterior-4', 'Baseline');
%legend('Posterior-32', 'Posterior-16', 'Posterior-4', 'Posterior-2', 'Baseline');

h2 = figure; 
for k=1:numTFs,
  figure(h2);
  subplot(2, 5, k);
  drosPlotAccuracyBars({J_indiv32{k}, J_indiv16{k}, J_indiv2{k}, J_indbase{k}}, M(:, k), T);
  title(sprintf('%s', drosTF.names{k}));
end
subplot(2, 5, 8);
bar(rand(4));
axis([-10 -9 -10 -9]);
axis off;
legend('Posterior-32', 'Posterior-16', 'Posterior-4/2(pair/single)', 'Baseline');
%legend('Posterior-32', 'Posterior-16',  'Posterior-2', 'Baseline');

% print Plots
property = 'Constrained';
if flag ~= 1
    property = 'Unconstrained'; 
end
dd = date;
if printPlot 
   print(h1, '-depsc2', [ddir 'drosophilaBars_' 'PairsOfTFs_', property dd '.eps']);
   print(h2, '-depsc2', [ddir 'drosophilaBars_' 'SingleTFs_', property dd '.eps']); 
end

