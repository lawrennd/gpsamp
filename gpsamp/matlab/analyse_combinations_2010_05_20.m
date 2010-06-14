%results = load('results/multitf3_2010-03-05_summary.mat');
results_a = sortResults(load('results/multitf5a_2010-04-22_summary2.mat'));
results_b = sortResults(load('results/multitf5a_2010-04-22_summary.mat'));
results_a2 = sortResults(load('results/multitf5a_null_2010-06-07_summary2.mat'));
results_b2 = sortResults(load('results/multitf5a_null_2010-06-07_summary.mat'));
results_a.marlls(:, 1) = results_a2.marlls;
results_b.marlls(:, 1) = results_b2.marlls;

results_ca = sortResults(load('results/multitf5c_2010-04-22_summary2.mat'));
results_cb = sortResults(load('results/multitf5c_2010-04-22_summary.mat'));
results_ca2 = sortResults(load('results/multitf5c_null_2010-06-07_summary2.mat'));
results_cb2 = sortResults(load('results/multitf5c_null_2010-06-07_summary.mat'));
results_ca.marlls(:, 1, :) = results_ca2.marlls;
results_cb.marlls(:, 1, :) = results_cb2.marlls;

cvals_a = results_a.marlls - mean(results_ca.marlls, 3);
cvals_b = results_b.marlls - mean(results_cb.marlls, 3);

M = drosMakeValidationMatrix(chipdistances, results_a.genes, 2000);

[foo, J_basic] = sort(results_a.marlls(:,4) - ...
		      logsumexp(results_a.marlls(:,[3, 2, 1]),2), 'descend');

[foo, J_cv] = sort(cvals_a(:, 4) - max(cvals_a(:, 1:3), [], 2), 'descend');
%[foo, J_cv2] = sort(cvals_a(:, 4) - logsumexp(cvals_a(:, 1:3), 2), 'descend');
[foo, J_cv2] = sort(cvals_a(:, 4) ...
		    - logsumexp(results_a.marlls(:, 1:3), 2) ...
		    + mean(logsumexp(results_ca.marlls(:, 1:3, :), 2), 3), 'descend');

% pnasrankings = load('~/mlprojects/disimrank/matlab/results/rankings.mat');

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

drosPlotAccuracyBars({J_basic, J_cv, J_cv2}, prod(M(:, [3, 5]), 2), T)
