results{1} = sortResults(load('results/multitfToyOneCond6a_2010-07-09_summary.mat'));
results{2} = sortResults(load('results/multitfToyTwoConds6a_2010-07-09_summary.mat'));
load toy4TFs28_June_10.mat

tcomb = [0 0 0; 1 0 0; 0 1 0; 0 0 1;
	 1 1 0; 1 0 1; 0 1 1;
	 1 1 1];

auc = {};
M = Net(results{1}.genes, 1:3);
linkprobs = {};
I = {};
val = {};
val1 = {};
figure(1);
for k=1:2,
  subplot(1, 2, k);
  for l=1:3,
    linkprobs{k}(:, l) = logsumexp(results{k}.marlls(:, tcomb(:, l)==1), 2) - ...
	logsumexp(results{k}.marlls(:, tcomb(:, l)==0), 2);
  end
  p = linkprobs{k}(1:500, :);
  [foo, I{k}] = sort(p(:), 'descend');
  N = M(1:500, :);
  val1{k} = N(I{k});
  auc{1}(k) = drosPlotROC(val1{k}', sum(val1{k}), length(val1{k}), 'r');
  hold on
  [foo, I{k}] = sort(linkprobs{k}(:), 'descend');
  val{k} = M(I{k});
  auc{2}(k) = drosPlotROC(val{k}', sum(val{k}), length(val{k}));
  hold off
end

modelprobs = {};
modelinds = {};
for k=1:2,
  subplot(1, 2, k);
  [foo, modelinds{k}] = max(results{k}.marlls, [], 2);
  J = sub2ind(size(results{k}.marlls), 1:length(results{k}.genes), modelinds{k}');
  modelprobs{k} = results{k}.marlls(J') - ...
	logsumexp(results{k}.marlls, 2);

  val{k} = all(M == tcomb(modelinds{k}, :), 2);
  [foo, I1{k}] = sort(modelprobs{k}(1:500), 'descend');
  mean(val{k}(I1{k}))
  [foo, I2{k}] = sort(modelprobs{k}, 'descend');
  mean(val{k}(I2{k}))
end
