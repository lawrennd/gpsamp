%results = load('results/multitf3_2010-03-05_summary.mat');
results_a = sortResults(load('results/multitf5a_2010-04-22_summary2.mat'));
results_b = sortResults(load('results/multitf5a_2010-04-22_summary.mat'));
results_ca = sortResults(load('results/multitf5c_2010-04-22_summary2.mat'));
results_cb = sortResults(load('results/multitf5c_2010-04-22_summary.mat'));

cvals_a = results_a.marlls - mean(results_ca.marlls, 3);
cvals_b = results_b.marlls - mean(results_cb.marlls, 3);

M = drosMakeValidationMatrix(chipdistances, results_a.genes, 2000);

[foo, J_basic] = sort(results_a.marlls(:,4) - ...
		      logsumexp(results_a.marlls(:,[3, 2, 1]),2), 'descend');

[foo, J_cv] = sort(cvals_a(:, 4) - max(cvals_a(:, 1:3), [], 2), 'descend');

pnasrankings = load('~/mlprojects/disimrank/matlab/results/rankings.mat');

val_basic = prod(M(J_basic, [3,5]), 2); val_basic = val_basic(~isnan(val_basic));
val_cv = prod(M(J_cv, [3,5]), 2); val_cv = val_cv(~isnan(val_cv));

drosPlotROC(val_basic', sum(val_basic), length(val_basic), 'b')
hold on
drosPlotROC(val_cv', sum(val_cv), length(val_cv), 'm')
hold off
