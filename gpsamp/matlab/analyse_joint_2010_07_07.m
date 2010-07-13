%results = load('results/multitf3_2010-03-05_summary.mat');
results_b = sortResults(load('results/multitf5a_2010-04-22_summary.mat'));
results_b2 = sortResults(load('results/multitf5a_null_2010-06-07_summary.mat'));
results_b3 = sortResults(load('results/multitf6a_2010-06-24_summary.mat'));
results_b.marlls(:, 1) = results_b2.marlls;
results_b.marlls = [results_b.marlls, results_b3.marlls];

baseline_a = sortResults(load('results/multitf5a_baseline_2010-07-02_summary.mat'));
baseline_a2 = sortResults(load('results/multitf5a_baseline2_2010-07-02_summary.mat'));

baseline_a.marlls = [baseline_a2.marlls, baseline_a.marlls];

comb = [0 0 0 0 0; 0 0 1 0 0; 0 0 0 0 1; 0 0 1 0 1;
	1 0 0 0 0; 0 1 0 0 0; 0 0 0 1 0;
	1 1 0 0 0; 1 0 1 0 0; 1 0 0 1 0; 1 0 0 0 1;
	0 1 1 0 0; 0 1 0 1 0; 0 1 0 0 1;
	0 0 1 1 0;
	0 0 0 1 1];

baselinecomb = [0 0 0 0 0;
		1 0 0 0 0; 0 1 0 0 0; 0 0 1 0 0 ; 0 0 0 1 0; 0 0 0 0 1;
		1 1 0 0 0; 1 0 1 0 0; 1 0 0 1 0; 1 0 0 0 1;
		0 1 1 0 0; 0 1 0 1 0; 0 1 0 0 1;
		0 0 1 1 0; 0 0 1 0 1;
		0 0 0 1 1];

combinds = find(sum(comb, 2) == 2);

M = drosMakeValidationMatrix(chipdistances, results_b.genes, 2000);

I = find(sum(comb, 2) == 1);

% Order the baseline results the same as the others
[foo, IA, IB] = intersect(comb, baselinecomb, 'rows');
[foo, J] = sort(IB);
combbaseinds = IA(J);

prior = zeros(1, size(comb, 1));
for k=1:size(comb, 1),
  prior(k) = mean(all(M(:, comb(k, :)==1), 2));
end
prioracc = prior;
prior = prior / sum(prior);

prioraccs(1) = mean(prioracc(2:end));
prioraccs(2) = sum(prioracc(2:end) .* prior(2:end)) / sum(prior(2:end));

posteriors = {};
posteriors{1} = results_b.marlls;
posteriors{2} = results_b.marlls + repmat(log(prior), [size(results_b.marlls, 1), 1]);

T = [200, 400, 800, 1600, 3200, 6000];

prioraccinds = [1, 2, 1];
posteriors{3} = baseline_a.marlls - repmat(baseline_a.marlls(:, 1), ...
					   [1, size(baseline_a.marlls, 2)]);

scores = {}; I={}; sscores={}; J={};
r = zeros(3, length(T));
r2 = zeros(3, length(T));
pvals = zeros(3, length(T));
for k=1:3,
  if k < 3,
    posteriors{k} = posteriors{k} - repmat(logsumexp(posteriors{k}, 2), [1, 16]);
    mycomb = comb;
  else
    mycomb = baselinecomb;
  end

  [scores{k}, I{k}] = max(posteriors{k}, [], 2);
  [sscores{k}, J{k}] = sort(scores{k}, 'descend');

  for l=1:length(T),
    acc = 0; acc2 = 0;
    count = 0;
    for m=1:T(l),
      if I{k}(J{k}(m)) ~= 1,
	acc = acc + all(M(J{k}(m), mycomb(I{k}(J{k}(m)), :)==1), 2);
	acc2 = acc2 + all(M(J{k}(m), :) == mycomb(I{k}(J{k}(m)), :), 2);
	count = count + 1;
      end
    end
    r(k, l) = acc / count;
    r2(k, l) = acc2 / count;
    pvals(k, l) = 1 - binocdf(acc - 1, count, prioraccs(prioraccinds(k)));
  end
end

figure(1);
h = bar(100*r');
set(gca, 'XTickLabel', T);
hold on
v = axis;
v(3:4) = [0 50];
axis(v)
plot(v(1:2), 100*prioraccs(1)*[1 1], 'b');
plot(v(1:2), 100*prioraccs(2)*[1 1], 'g');
hold off
legend('No prior', 'Empirical prior', 'Baseline', 'Uniform random', 'Random from prior', 'Location', 'EastOutside')
axis(v)
xlabel('# of top genes')
ylabel('Enrichment (%)')

drosStarBars(h, pvals');


scores = {}; I={}; sscores={}; J={};
r2 = zeros(3, length(T));
pvals2 = zeros(3, length(T));
for k=1:3,
  if k < 3,
    posteriors{k} = posteriors{k} - repmat(logsumexp(posteriors{k}, 2), [1, 16]);
    mycomb = comb;
  else
    mycomb = baselinecomb;
  end

  [scores{k}, I{k}] = max(posteriors{k}(:, 1:6), [], 2);
  [sscores{k}, J{k}] = sort(scores{k}, 'descend');

  for l=1:length(T),
    acc = 0;
    count = 0;
    for m=1:T(l),
      if I{k}(J{k}(m)) ~= 1,
	acc = acc + all(M(J{k}(m), mycomb(I{k}(J{k}(m)), :)==1), 2);
	count = count + 1;
      end
    end
    r2(k, l) = acc / count;
    pvals2(k, l) = 1 - binocdf(acc - 1, count, prioraccs(prioraccinds(k)));
  end
end



T = [100, 200, 400, 800, 1600, 3200];

% Single TF probabilities
% 1 = posterior many
% 2 = posterior single

withprior = bsxfun(@plus, results_b.marlls, log(prior));

sposteriors = {};
for k=1:5,
  sposteriors{1}(:, k) = logsumexp(results_b.marlls(:, comb(:, k)==1), 2) - ...
      logsumexp(results_b.marlls(:, comb(:, k)==0), 2);
  sposteriors{2}(:, k) = logsumexp(withprior(:, comb(:, k)==1), 2) - ...
      logsumexp(withprior(:, comb(:, k)==0), 2);
  E = zeros(1, 5);
  E(k) = 1;
  sposteriors{3}(:, k) = results_b.marlls(:, all(bsxfun(@eq, comb, E), 2)) - ...
      results_b.marlls(:, 1);
  sposteriors{4}(:, k) = withprior(:, all(bsxfun(@eq, comb, E), 2)) - ...
      withprior(:, 1);
  sposteriors{5}(:, k) = baseline_a.marlls(:, k+1) - ...
      baseline_a.marlls(:, 1);
end

prior2 = nanmean(M);
prioraccs2(1) = mean(prior2);
prioraccs2(2) = sum(prior2 .* prior2) / sum(prior2).^2;

baselines2 = [1, 2, 1, 2, 1];

r3 = zeros(length(T), length(sposteriors));
pvals3 = r3;
for k=1:length(sposteriors),
  [foo, I] = sort(sposteriors{k}(:), 'descend');
  for l=1:length(T),
    r3(l, k) = nanmean(M(I(1:T(l))));
    pvals3(l, k) = 1 - binocdf(nansum(M(I(1:T(l)))) - 1, ...
			       sum(~isnan(M(I(1:T(l))))), ...
			       prioraccs2(baselines2(k)));
  end
end

figure(2);
h2 = bar(100*r3);
set(gca, 'XTickLabel', T);
hold on
v = axis;
v(3:4) = [0 50];
axis(v)
plot(v(1:2), 100*prioraccs2(1)*[1 1], 'b');
plot(v(1:2), 100*prioraccs2(2)*[1 1], 'g');
hold off
legend('Posterior many', '+prior', 'Posterior single', '+prior', ...
       'Baseline', 'Uniform random', 'Random from prior', ...
       'Location', 'EastOutside')
axis(v)
xlabel('# of top predictions')
ylabel('Enrichment (%)')

drosStarBars(h2, pvals3);
