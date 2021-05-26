%zthresholds = 0:.1:9;
%for t=1:length(zthresholds),
%  zacc(t) = nanmean(M(results.zscores' > zthresholds(t)));
%end

pplrthresholds = 0:.1:1;
for t=1:length(pplrthresholds),
  val = M(results.pplrs' > pplrthresholds(t));
  pplracc(t) = nanmean(val);
  ppval(t) = 1 - hygecdf(nansum(val)-1, sum(~isnan(M(:))), nansum(M(:)), sum(~isnan(val)));
end

scorethresholds = 10.^[-5:.5:0];
for t=1:length(scorethresholds),
  val = M(scores > scorethresholds(t));
  scoreacc(t) = nanmean(val);
  spval(t) = 1 - hygecdf(nansum(val)-1, sum(~isnan(M(:))), nansum(M(:)), sum(~isnan(val)));
end

subplot(4, 2, 1);

loglog(scorethresholds, spval)
xlabel('score threshold')
title('p-values')
%plot(pplrthresholds, pplracc)
%xlabel('pplr score threshold')
%ylabel('frequency of positive ChIP validation')

subplot(4, 2, 2);
semilogx(scorethresholds, scoreacc)
xlabel('score threshold')
title('accuracies')

% subplot(4, 3, 3);
% regthresholds = 0.9:.001:1;
% for t=1:length(regthresholds),
%   val = M(regression_probs > regthresholds(t));
%   regacc(t) = nanmean(val);
%   rpval(t) = 1 - hygecdf(nansum(val)-1, sum(~isnan(M(:))), nansum(M(:)), sum(~isnan(val)));
% end

% plot(regthresholds, regacc)
% xlabel('regression probability threshold')
% axis tight

TFs = {'bin', 'twi', 'mef2'};

for k=1:3,
  %for t=1:length(zthresholds),
  %  zaccs{k}(t) = nanmean(M(results.zscores(k, :) > zthresholds(t), k));
  %end

  for t=1:length(scorethresholds),
    val = M(scores(:, k) > scorethresholds(t), k);
    scoreaccs{k}(t) = nanmean(val);
    pvals{k}(t) = 1 - hygecdf(nansum(val)-1, sum(~isnan(M(:, k))), nansum(M(:, k)), sum(~isnan(val)));
  end
%  for t=1:length(pplrthresholds),
%    val = M(results.pplrs(k, :) > pplrthresholds(t), k);
%    pplraccs{k}(t) = nanmean(val);
%    pvals{k}(t) = 1 - hygecdf(nansum(val)-1, sum(~isnan(M(:, k))), nansum(M(:, k)), sum(~isnan(val)));
%  end

  subplot(4, 2, 1 + 2*k);
  
  loglog(scorethresholds, pvals{k})
  xlabel('score threshold')
  ylabel(TFs{k})
  %ylabel('frequency of positive ChIP validation')

  subplot(4, 2, 2 + 2*k);

  semilogx(scorethresholds, scoreaccs{k})
  xlabel('score threshold')
  axis tight
  %ylabel('frequency of positive ChIP validation')

%   subplot(4, 3, 1 + 3*k);
  
%   semilogy(pplrthresholds, pvals{k})
%   xlabel('pplr score threshold')
%   ylabel(TFs{k})
%   %ylabel('frequency of positive ChIP validation')

%   subplot(4, 3, 2 + 3*k);

%   plot(pplrthresholds, pplraccs{k})
%   xlabel('pplr score threshold')
%   %ylabel('frequency of positive ChIP validation')

end
