% version 2010-07-01

results_b = sortResults(load('results/multitf5a_2010-04-22_summary.mat'));
results_b2 = sortResults(load('results/multitf5a_null_2010-06-07_summary.mat'));
results_b3 = sortResults(load('results/multitf6a_2010-06-24_summary.mat'));
results_b.marlls(:, 1) = results_b2.marlls;
results_b.marlls = [results_b.marlls, results_b3.marlls];

comb = [0 0 0 0 0; 0 0 1 0 0; 0 0 0 0 1; 0 0 1 0 1;
	1 0 0 0 0; 0 1 0 0 0; 0 0 0 1 0;
	1 1 0 0 0; 1 0 1 0 0; 1 0 0 1 0; 1 0 0 0 1;
	0 1 1 0 0; 0 1 0 1 0; 0 1 0 0 1;
	0 0 1 1 0;
	0 0 0 1 1];

combinds = find(sum(comb, 2) == 2);

I = find(sum(comb, 2) == 1);
tfI = comb(I, :) * (1:5)';

labels = {};
for k=2:size(comb, 1),
  switch (sum(comb(k, :)))
   case 1,
    labels{k-1} = drosTF.names{find(comb(k, :))};
   case 2,
    J = find(comb(k, :));
    labels{k-1} = sprintf('%s+%s', drosTF.names{J(1)}, drosTF.names{J(2)});
  end
end

I_chip = drosFindGeneinds(chipdistances, results_a.genes, 1)';
I_insitu = drosFindGeneinds(drosinsitu, results_a.genes, 1)';

hasinsitu = (I_insitu ~= 0);
isinsitu = zeros(size(hasinsitu));
isinsitu(hasinsitu) = any(drosinsitu.data(I_insitu(I_insitu ~= 0), :), 2);

DATAHEADINGS = {'FBgn', 'ischip', 'hasinsitu', 'isinsitu'};
SCOREHEADINGS = {'marll', 'posterior_single', 'posterior_many', 'baseline'};

for k=1:length(labels),
  ischip = zeros(size(hasinsitu));
  ischip(I_chip ~= 0) = all(chipdistances.data(I_chip(I_chip ~= 0), find(comb(k+1, :))) < 2000, 2);
  scores = zeros(length(ischip), 4);
  scores(:, 1) = results_b.marlls(:, k+1);
  scores(:, 2) = results_b.marlls(:, k+1) - ...
      results_b.marlls(:, 1);
  switch (sum(comb(k+1, :)))
   case 1,
    tfI = find(comb(k+1, :));
    scores(:, 3) = logsumexp(results_b.marlls(:,comb(:, tfI)==1),2) - ...
	logsumexp(results_b.marlls(:,comb(:, tfI)==0),2);
   case 2,
    scores(:, 3) = results_b.marlls(:, k+1) - ...
	logsumexp(results_b.marlls(:, setdiff(1:16, k+1)),2);
  end
  scores(:, 4) = results_b.marlls(:, 1);

  fid = fopen(sprintf('results_2010-07-01/dros_%s_results.txt', labels{k}), 'w');
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
