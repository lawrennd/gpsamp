%results = load('results/multitf3_2010-03-05_summary.mat');
results_a = sortResults(load('results/multitf5a_2010-04-22_summary2.mat'));
results_b = sortResults(load('results/multitf5a_2010-04-22_summary.mat'));
results_a2 = sortResults(load('results/multitf5a_null_2010-06-07_summary2.mat'));
results_b2 = sortResults(load('results/multitf5a_null_2010-06-07_summary.mat'));
results_b3 = sortResults(load('results/multitf6a_2010-06-24_summary.mat'));
results_a.marlls(:, 1) = results_a2.marlls;
results_b.marlls(:, 1) = results_b2.marlls;
results_b.marlls = [results_b.marlls, results_b3.marlls];

comb = [0 0 0 0 0; 0 0 1 0 0; 0 0 0 0 1; 0 0 1 0 1;
	1 0 0 0 0; 0 1 0 0 0; 0 0 0 1 0;
	1 1 0 0 0; 1 0 1 0 0; 1 0 0 1 0; 1 0 0 0 1;
	0 1 1 0 0; 0 1 0 1 0; 0 1 0 0 1;
	0 0 1 1 0;
	0 0 0 1 1];

combinds = find(sum(comb, 2) == 2);

results_ca = sortResults(load('results/multitf5c_2010-04-22_summary2.mat'));
results_cb = sortResults(load('results/multitf5c_2010-04-22_summary.mat'));
results_ca2 = sortResults(load('results/multitf5c_null_2010-06-07_summary2.mat'));
results_cb2 = sortResults(load('results/multitf5c_null_2010-06-07_summary.mat'));
results_ca.marlls(:, 1, :) = results_ca2.marlls;
results_cb.marlls(:, 1, :) = results_cb2.marlls;

cvals_a = results_a.marlls - mean(results_ca.marlls, 3);
cvals_b = results_b.marlls(:, 1:4) - mean(results_cb.marlls, 3);

M = drosMakeValidationMatrix(chipdistances, results_a.genes, 2000);

I = find(sum(comb, 2) == 1);
tfI = comb(I, :) * (1:5)';

J_joint = {};
J_joint0 = {};
J_indep = {};
for k=1:length(combinds),
  [foo, J_joint{k}] = sort(results_b.marlls(:,combinds(k)) - ...
			   logsumexp(results_b.marlls(:,setdiff(1:size(results_b.marlls, 2), combinds(k))),2), 'descend');
  [foo, J_joint0{k}] = sort(results_b.marlls(:,combinds(k)) - ...
			   results_b.marlls(:,1), 'descend');
  tfs = find(comb(combinds(k), :));
  
  [foo, J_indep{k}] = sort(results_b.marlls(:,I(tfI==tfs(1))) + ...
			   results_b.marlls(:,I(tfI==tfs(2))) - ...
			   logsumexp(results_b.marlls(:, ...
						  setdiff([1; I], ...
						  [I(tfI==tfs(1)), I(tfI==tfs(2))])),2), 'descend');
end

J_indiv = {};
J_indiv0 = {};
J_indiv00 = {};
J_old = {};
for k=1:length(I),
  [foo, J_indiv{k}] = sort(logsumexp(results_b.marlls(:,comb(:, tfI(k))==1),2) - ...
			   logsumexp(results_b.marlls(:,comb(:, tfI(k))==0),2), ...
			   'descend');
  [foo, J_indiv0{k}] = sort(logsumexp(results_b.marlls(:,I(tfI==k)),2) - ...
			   logsumexp(results_b.marlls(:,1),2), ...
			   'descend');
  [foo, J_indiv00{k}] = sort(logsumexp(results_b.marlls(:,I(tfI==k)),2) - ...
			   0*logsumexp(results_b.marlls(:,1),2), ...
			   'descend');
end
posterior_indices = {[2, 4], [3, 4]};
for k=1:2,
  [foo, J_old{k}] = sort(logsumexp(results_b.marlls(:,posterior_indices{k}),2) - ...
			 logsumexp(results_b.marlls(:, setdiff(1:4, posterior_indices{k})),2), ...
			 'descend');
end
for k=3:5,
  J_old{k} = [];
end


[foo, J_basic] = sort(results_b.marlls(:,4) - ...
		      logsumexp(results_b.marlls(:,[3, 2, 1]),2), 'descend');

[foo, J_cv] = sort(cvals_b(:, 4) - max(cvals_b(:, 1:3), [], 2), 'descend');
%[foo, J_cv2] = sort(cvals_b(:, 4) - logsumexp(cvals_b(:, 1:3), 2), 'descend');
[foo, J_cv2] = sort(cvals_b(:, 4) ...
		    - logsumexp(results_b.marlls(:, 1:3), 2) ...
		    + mean(logsumexp(results_cb.marlls(:, 1:3, :), 2), 3), 'descend');

pnasrankings = load('~/mlprojects/disimrank/matlab/results/rankings.mat');

% val_basic = prod(M(J_basic, [3,5]), 2); val_basic = val_basic(~isnan(val_basic));
% val_cv = prod(M(J_cv, [3,5]), 2); val_cv = val_cv(~isnan(val_cv));

% drosPlotROC(val_basic', sum(val_basic), length(val_basic), 'b')
% hold on
% drosPlotROC(val_cv', sum(val_cv), length(val_cv), 'm')
% hold off

% acc_cv = cumsum(val_cv) ./ (1:length(val_cv))';
% acc_basic = cumsum(val_basic) ./ (1:length(val_basic))';

T = [20, 50, 100, 150, 200];

%clf;
%plot(T, acc_cv(T), 'm-*');
%hold on
%plot(T, acc_basic(T), 'b-s');

for k=1:length(combinds),
  figure(1);
  subplot(2, 5, k);
  drosPlotAccuracyBars({J_joint{k}, J_joint0{k}, J_indep{k}}, prod(M(:, comb(combinds(k), :)==1), 2), T);
  tfs = find(comb(combinds(k), :));
  title(sprintf('%s + %s', drosTF.names{tfs(1)}, drosTF.names{tfs(2)}));
  %drosPlotAccuracyBars({J_basic, J_cv, J_cv2}, prod(M(:, [3, 5]), 2), T)
end

for k=1:length(I),
  figure(2);
  subplot(1, 5, tfI(k));
  drosPlotAccuracyBars({J_indiv{k}, J_indiv0{k}, J_indiv00{k}}, M(:, tfI(k)), T);
  title(sprintf('%s', drosTF.names{tfI(k)}));
  %drosPlotAccuracyBars({J_basic, J_cv, J_cv2}, prod(M(:, [3, 5]), 2), T)
end
