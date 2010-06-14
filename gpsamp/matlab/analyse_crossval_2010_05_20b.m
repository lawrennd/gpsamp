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

resultinds = zeros(size(results_a.genes));
for k=1:length(results_a.genes),
  J = drosFindGeneinds(drosexp, results_a.genes(k));
  resultinds(k) = J(1);
end

%post = drosComputePosterior(results_a.marlls);
%twiprob = sum(post(:, [2, 4]), 2);
%mefprob = sum(post(:, [3, 4]), 2);
%M = drosMakeValidationMatrix(chipdistances, results_a.genes, 2000);
%M = drosMakeValidationMatrix(chipdistances, results_a.genes, 2000);
%M_drosexp = drosMakeValidationMatrix(chipdistances, drosexp.genes, 2000);
%nanmean(M(twiprob > .999, :))
%nanmean(M(mefprob > .999, :))

%nanmean(M(twiprob > .99, :))
%nanmean(M(mefprob > .99, :))

%nanmean(M(twiprob > .95, :))
%nanmean(M(mefprob > .95, :))

%[foo, J_twi] = sort(twiprob, 'descend');
%[foo, J_mef] = sort(mefprob, 'descend');

%[foo, multitfrank1.twi] = sort(logsumexp(cvals_a(:,2),2) - ...
%			       0*logsumexp(cvals_a(:,1),2), 'descend');
%[foo, multitfrank1.mef2] = sort(logsumexp(cvals_a(:,3),2) - ...
%				0*logsumexp(cvals_a(:,1),2), 'descend');
[foo, multitfrank1.twi] = sort(logsumexp(results_a.marlls(:,[4, 2]),2) - ...
		     logsumexp(results_a.marlls(:,[3, 1]),2), 'descend');
[foo, multitfrank1.mef2] = sort(logsumexp(results_a.marlls(:,[4, 3]),2) - ...
		     logsumexp(results_a.marlls(:,[2, 1]),2), 'descend');
%[foo, J_twi] = sort(logsumexp(results_a.marlls(:,[2]),2) - ...
%		     0*logsumexp(results_a.marlls(:,[3, 1]),2), 'descend');
%[foo, J_mef] = sort(logsumexp(results_a.marlls(:,[3]),2) - ...
%		     0*logsumexp(results_a.marlls(:,[2, 1]),2), 'descend');

% val_twi = M(J_twi, 3); val_twi = val_twi(~isnan(val_twi));
% val_mef = M(J_mef, 5); val_mef = val_mef(~isnan(val_mef));

% subplot(1, 2, 1)
% drosPlotROC(val_twi', sum(val_twi), length(val_twi))
% hold on;
% subplot(1, 2, 2)
% drosPlotROC(val_mef', sum(val_mef), length(val_mef))
% hold on;

% twipost = drosComputePosterior(results_a.marlls(:, [1,2]));
% mefpost = drosComputePosterior(results_a.marlls(:, [1,3]));
%[foo, J2_twi] = sort(twipost(:, 2), 'descend');
%[foo, J2_mef] = sort(mefpost(:, 2), 'descend');
% [foo, multitfrank2.twi] = sort(logsumexp(cvals_a(:,[2,4]),2) - ...
% 			       logsumexp(cvals_a(:,[1,3]),2), 'descend');
% [foo, multitfrank2.mef2] = sort(logsumexp(cvals_a(:,[3,4]),2) - ...
% 				logsumexp(cvals_a(:,[1,2]),2), 'descend');

posterior_indices = [2,4];
scores = logsumexp(results_a.marlls(:,posterior_indices),2) - ...
	 logsumexp(results_a.marlls(:, setdiff(1:4, posterior_indices)),2) - ...
	 mean(logsumexp(results_ca.marlls(:, posterior_indices, :), 2), 3) + ...
	 mean(logsumexp(results_ca.marlls(:, setdiff(1:4, posterior_indices), :), 2), 3);

[foo, multitfrank2.twi] = sort(scores, 'descend');

posterior_indices = [3,4];
scores = logsumexp(results_a.marlls(:,posterior_indices),2) - ...
	 logsumexp(results_a.marlls(:, setdiff(1:4, posterior_indices)),2) - ...
	 mean(logsumexp(results_ca.marlls(:, posterior_indices, :), 2), 3) + ...
	 mean(logsumexp(results_ca.marlls(:, setdiff(1:4, posterior_indices), :), 2), 3);

[foo, multitfrank2.mef2] = sort(scores, 'descend');
% [foo, multitfrank2.twi] = sort(max(cvals_a(:,[2,4]),[], 2) - ...
% 			       max(cvals_a(:,[1,3]),[], 2), 'descend');
% [foo, multitfrank2.mef2] = sort(max(cvals_a(:,[3,4]),[], 2) - ...
% 				max(cvals_a(:,[1,2]),[], 2), 'descend');


% [foo, multitfrank2.twi] = sort(logsumexp(results_b.marlls(:,[4, 2]),2) - ...
% 		     logsumexp(results_b.marlls(:,[3, 1]),2), 'descend');
% [foo, multitfrank2.mef2] = sort(logsumexp(results_b.marlls(:,[4, 3]),2) - ...
% 		     logsumexp(results_b.marlls(:,[2, 1]),2), 'descend');
%[foo, J2_twi] = sort(results_a.marlls(:, 4) - results_a.marlls(:, 3), 'descend');
%[foo, J2_mef] = sort(results_a.marlls(:, 4) - results_a.marlls(:, 2), 'descend');
% [foo, J2_twi] = sort(logsumexp(results_a.marlls(:,[4, 2]),2) - ...
%     logsumexp(results_a.marlls(:,[3, 1]),2), 'descend');
% [foo, J2_mef] = sort(logsumexp(results_a.marlls(:,[4, 3]),2) - ...
%     logsumexp(results_a.marlls(:,[2, 1]),2), 'descend');

% jointpost = results_a.marlls(:,4) - ...
%     logsumexp(results_a.marlls(:,[3, 2, 1]),2);

% val2_twi = M(J2_twi, 3); val2_twi = val2_twi(~isnan(val2_twi));
% val2_mef = M(J2_mef, 5); val2_mef = val2_mef(~isnan(val2_mef));

% subplot(1, 2, 1)
% drosPlotROC(val2_twi', sum(val2_twi), length(val2_twi), 'g')
% subplot(1, 2, 2)
% drosPlotROC(val2_mef', sum(val2_mef), length(val2_mef), 'g')

%[foo, multitfrank1.twi] = sort(logsumexp(cvals_a(:,[4,2]),2) - ...
%			       logsumexp(cvals_a(:,[1,3]),2), 'descend');
%[foo, multitfrank1.mef2] = sort(logsumexp(cvals_a(:,[4,3]),2) - ...
%				logsumexp(cvals_a(:,[1,2]),2), 'descend');
% [foo, multitfrank3.twi] = sort(logsumexp(results_a.marlls(:,[2]),2) - ...
% 		     logsumexp(results_a.marlls(:,[1]),2), 'descend');
% [foo, multitfrank3.mef2] = sort(logsumexp(results_a.marlls(:,[3]),2) - ...
% 		     logsumexp(results_a.marlls(:,[1]),2), 'descend');
%[foo, J3_twi] = sort(results_a.marlls(:, 2), 'descend');
%[foo, J3_mef] = sort(results_a.marlls(:, 3), 'descend');

% val3_twi = M(J3_twi, 3); val3_twi = val3_twi(~isnan(val3_twi));
% val3_mef = M(J3_mef, 5); val3_mef = val3_mef(~isnan(val3_mef));

% subplot(1, 2, 1)
% drosPlotROC(val3_twi', sum(val3_twi), length(val3_twi), 'r')
% subplot(1, 2, 2)
% drosPlotROC(val3_mef', sum(val3_mef), length(val3_mef), 'r')

pnasrankings = load('~/mlprojects/disimrank/matlab/results/rankings.mat');
% [A, B, C] = intersect(drosexp.genes(rankings.disimrank.twi), results_a.genes);
% twirank = rankings.disimrank.twi(sort(B));
% pnasval_twi = M_drosexp(twirank, 3); pnasval_twi = pnasval_twi(~isnan(pnasval_twi));

% [A, B, C] = intersect(drosexp.genes(rankings.disimrank.mef2), results_a.genes);
% mefrank = rankings.disimrank.mef2(sort(B));
% pnasval_mef = M_drosexp(mefrank, 5); pnasval_mef = pnasval_mef(~isnan(pnasval_mef));

% [A, B, C] = intersect(drosexp.genes(rankings.indrank.twi), results_a.genes);
% twirank = rankings.indrank.twi(sort(B));
% pnasval2_twi = M_drosexp(twirank, 3); pnasval2_twi = pnasval2_twi(~isnan(pnasval2_twi));

% [A, B, C] = intersect(drosexp.genes(rankings.indrank.mef2), results_a.genes);
% mefrank = rankings.indrank.mef2(sort(B));
% pnasval2_mef = M_drosexp(mefrank, 5); pnasval2_mef = pnasval2_mef(~isnan(pnasval2_mef));

% subplot(1, 2, 1)
% drosPlotROC(pnasval_twi', sum(pnasval_twi), length(pnasval_twi), 'm')
% plot([0 1], [0 1], 'k--')
% axis([0 .2 0 .4]);
% hold off;
% subplot(1, 2, 2)
% drosPlotROC(pnasval_mef', sum(pnasval_mef), length(pnasval_mef), 'm')
% plot([0 1], [0 1], 'k--')
% axis([0 .2 0 .4]);
% hold off;

multitfrank1.twi = resultinds(multitfrank1.twi);
multitfrank1.mef2 = resultinds(multitfrank1.mef2);
multitfrank2.twi = resultinds(multitfrank2.twi);
multitfrank2.mef2 = resultinds(multitfrank2.mef2);
%multitfrank3.twi = resultinds(multitfrank3.twi);
%multitfrank3.mef2 = resultinds(multitfrank3.mef2);

drosPlotEvaluation3
