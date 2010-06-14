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

resultinds = zeros(size(results_a.genes));
for k=1:length(results_a.genes),
  J = drosFindGeneinds(drosexp, results_a.genes(k));
  resultinds(k) = J(1);
end

labels = {'twi', 'mef2', 'both'};
marll_indices = {2, 3, 4};
posterior_indices = {[2, 4], [3, 4], [4]};
chipindices = {5, 3, [3, 5]};

I_chip = drosFindGeneinds(chipdistances, results_a.genes, 1)';
I_insitu = drosFindGeneinds(drosinsitu, results_a.genes, 1)';

hasinsitu = (I_insitu ~= 0);
isinsitu = zeros(size(hasinsitu));
isinsitu(hasinsitu) = any(drosinsitu.data(I_insitu(I_insitu ~= 0), :), 2);

DATAHEADINGS = {'FBgn', 'ischip', 'hasinsitu', 'isinsitu'};
SCOREHEADINGS = {'marll1', 'marll2', 'posterior1', 'posterior2', 'crossval1', 'crossval2', 'crosspost1', 'crosspost2', 'baseline1', 'baseline2'};

for k=1:length(labels),
  ischip = zeros(size(hasinsitu));
  ischip(I_chip ~= 0) = all(chipdistances.data(I_chip(I_chip ~= 0), chipindices{k}) < 2000, 2);
  scores = zeros(length(ischip), 10);
  scores(:, 1) = results_b.marlls(:, marll_indices{k});
  scores(:, 2) = results_a.marlls(:, marll_indices{k});
  scores(:, 3) = logsumexp(results_b.marlls(:,posterior_indices{k}),2) - ...
      logsumexp(results_b.marlls(:, setdiff(1:4, posterior_indices{k})),2);
      %logsumexp(results_b.marlls,2);
  scores(:, 4) = logsumexp(results_a.marlls(:,posterior_indices{k}),2) - ...
      logsumexp(results_a.marlls(:, setdiff(1:4, posterior_indices{k})),2);
      %logsumexp(results_a.marlls,2);
  scores(:, 5) = cvals_b(:, marll_indices{k});
  scores(:, 6) = cvals_a(:, marll_indices{k});
  scores(:, 7) = scores(:, 3) - ...
      mean(logsumexp(results_cb.marlls(:, posterior_indices{k}, :), 2), 3) + ...
      mean(logsumexp(results_cb.marlls(:, setdiff(1:4, posterior_indices{k}), :), 2), 3);
      %logsumexp(cvals_b,2);
  scores(:, 8) = scores(:, 4) - ...
      mean(logsumexp(results_ca.marlls(:, posterior_indices{k}, :), 2), 3) + ...
      mean(logsumexp(results_ca.marlls(:, setdiff(1:4, posterior_indices{k}), :), 2), 3);
      %logsumexp(cvals_a,2);
  scores(:, 9) = results_b.marlls(:, 1);
  scores(:, 10) = results_a.marlls(:, 1);

  fid = fopen(sprintf('dros_%s_results.txt', labels{k}), 'w');
  fprintf(fid, '%s ', DATAHEADINGS{:});
  fprintf(fid, '%s ', SCOREHEADINGS{:});
  fprintf(fid, '\n');
  for l=1:length(hasinsitu),
    fprintf(fid, '%s %d %d %d', results_a.genes{l}, ischip(l), hasinsitu(l), isinsitu(l));
    fprintf(fid, ' %f', scores(l, :));
    fprintf(fid, '\n');
  end
  fclose(fid);
end
