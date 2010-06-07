results = load('results/multitf3_2010-03-05_summary.mat');
results_a = load('results/multitf5a_2010-04-22_summary.mat');
results_b = load('results/multitf5b_2010-04-22_summary.mat');
post = drosComputePosterior(results_a.marlls);
twiprob = sum(post(:, [2, 4]), 2);
mefprob = sum(post(:, [3, 4]), 2);
%M = drosMakeValidationMatrix(chipdistances, results_a.genes, 2000);
M = drosMakeValidationMatrix(chipdistances, results_a.genes, 2000);
M_drosexp = drosMakeValidationMatrix(chipdistances, drosexp.genes, 2000);
%nanmean(M(twiprob > .999, :))
%nanmean(M(mefprob > .999, :))

%nanmean(M(twiprob > .99, :))
%nanmean(M(mefprob > .99, :))

%nanmean(M(twiprob > .95, :))
%nanmean(M(mefprob > .95, :))

%[foo, J_twi] = sort(twiprob, 'descend');
%[foo, J_mef] = sort(mefprob, 'descend');

[foo, J_twi] = sort(logsumexp(results_a.marlls(:,[4, 2]),2) - ...
		    logsumexp(results_a.marlls(:,[3, 1]),2), 'descend');
[foo, J_mef] = sort(logsumexp(results_a.marlls(:,[4, 3]),2) - ...
		    logsumexp(results_a.marlls(:,[2, 1]),2), 'descend');
%[foo, J_twi] = sort(logsumexp(results_a.marlls(:,[2]),2) - ...
%		     0*logsumexp(results_a.marlls(:,[3, 1]),2), 'descend');
%[foo, J_mef] = sort(logsumexp(results_a.marlls(:,[3]),2) - ...
%		     0*logsumexp(results_a.marlls(:,[2, 1]),2), 'descend');

val_twi = M(J_twi, 3); val_twi = val_twi(~isnan(val_twi));
val_mef = M(J_mef, 5); val_mef = val_mef(~isnan(val_mef));

subplot(1, 2, 1)
drosPlotROC(val_twi', sum(val_twi), length(val_twi))
hold on;
subplot(1, 2, 2)
drosPlotROC(val_mef', sum(val_mef), length(val_mef))
hold on;

twipost = drosComputePosterior(results_a.marlls(:, [1,2]));
mefpost = drosComputePosterior(results_a.marlls(:, [1,3]));
%[foo, J2_twi] = sort(twipost(:, 2), 'descend');
%[foo, J2_mef] = sort(mefpost(:, 2), 'descend');
[foo, J2_twi] = sort(logsumexp(results_b.marlls(:,[4, 2]),2) - ...
		     logsumexp(results_b.marlls(:,[3, 1]),2), 'descend');
[foo, J2_mef] = sort(logsumexp(results_b.marlls(:,[4, 3]),2) - ...
		     logsumexp(results_b.marlls(:,[2, 1]),2), 'descend');
%[foo, J2_twi] = sort(logsumexp(results_b.marlls(:,[2]),2) - ...
%		     0*logsumexp(results_b.marlls(:,[3, 1]),2), 'descend');
%[foo, J2_mef] = sort(logsumexp(results_b.marlls(:,[3]),2) - ...
%		     0*logsumexp(results_b.marlls(:,[2, 1]),2), 'descend');
%[foo, J2_twi] = sort(results_a.marlls(:, 4) - results_a.marlls(:, 3), 'descend');
%[foo, J2_mef] = sort(results_a.marlls(:, 4) - results_a.marlls(:, 2), 'descend');
% [foo, J2_twi] = sort(logsumexp(results_a.marlls(:,[4, 2]),2) - ...
%     logsumexp(results_a.marlls(:,[3, 1]),2), 'descend');
% [foo, J2_mef] = sort(logsumexp(results_a.marlls(:,[4, 3]),2) - ...
%     logsumexp(results_a.marlls(:,[2, 1]),2), 'descend');

jointpost = results_a.marlls(:,4) - ...
    logsumexp(results_a.marlls(:,[3, 2, 1]),2);

val2_twi = M(J2_twi, 3); val2_twi = val2_twi(~isnan(val2_twi));
val2_mef = M(J2_mef, 5); val2_mef = val2_mef(~isnan(val2_mef));

subplot(1, 2, 1)
drosPlotROC(val2_twi', sum(val2_twi), length(val2_twi), 'g')
subplot(1, 2, 2)
drosPlotROC(val2_mef', sum(val2_mef), length(val2_mef), 'g')

[foo, J3_twi] = sort(logsumexp(results_a.marlls(:,[2]),2) - ...
		     logsumexp(results_a.marlls(:,[1]),2), 'descend');
[foo, J3_mef] = sort(logsumexp(results_a.marlls(:,[3]),2) - ...
		     logsumexp(results_a.marlls(:,[1]),2), 'descend');
%[foo, J3_twi] = sort(results_a.marlls(:, 2), 'descend');
%[foo, J3_mef] = sort(results_a.marlls(:, 3), 'descend');

val3_twi = M(J3_twi, 3); val3_twi = val3_twi(~isnan(val3_twi));
val3_mef = M(J3_mef, 5); val3_mef = val3_mef(~isnan(val3_mef));

subplot(1, 2, 1)
drosPlotROC(val3_twi', sum(val3_twi), length(val3_twi), 'r')
subplot(1, 2, 2)
drosPlotROC(val3_mef', sum(val3_mef), length(val3_mef), 'r')

rankings = load('~/mlprojects/disimrank/matlab/results/rankings.mat');
[A, B, C] = intersect(drosexp.genes(rankings.disimrank.twi), results_a.genes);
twirank = rankings.disimrank.twi(sort(B));
pnasval_twi = M_drosexp(twirank, 3); pnasval_twi = pnasval_twi(~isnan(pnasval_twi));

[A, B, C] = intersect(drosexp.genes(rankings.disimrank.mef2), results_a.genes);
mefrank = rankings.disimrank.mef2(sort(B));
pnasval_mef = M_drosexp(mefrank, 5); pnasval_mef = pnasval_mef(~isnan(pnasval_mef));

[A, B, C] = intersect(drosexp.genes(rankings.indrank.twi), results_a.genes);
twirank = rankings.indrank.twi(sort(B));
pnasval2_twi = M_drosexp(twirank, 3); pnasval2_twi = pnasval2_twi(~isnan(pnasval2_twi));

[A, B, C] = intersect(drosexp.genes(rankings.indrank.mef2), results_a.genes);
mefrank = rankings.indrank.mef2(sort(B));
pnasval2_mef = M_drosexp(mefrank, 5); pnasval2_mef = pnasval2_mef(~isnan(pnasval2_mef));

subplot(1, 2, 1)
drosPlotROC(pnasval_twi', sum(pnasval_twi), length(pnasval_twi), 'm')
plot([0 1], [0 1], 'k--')
axis([0 .2 0 .4]);
hold off;
subplot(1, 2, 2)
drosPlotROC(pnasval_mef', sum(pnasval_mef), length(pnasval_mef), 'm')
plot([0 1], [0 1], 'k--')
axis([0 .2 0 .4]);
hold off;
