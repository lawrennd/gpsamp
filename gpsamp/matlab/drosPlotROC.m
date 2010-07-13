function auc = drosPlotROC(r, totpositives, totgenes, varargin),

tp = cumsum(r) / totpositives;
fp = cumsum(~r) / (totgenes - totpositives);

if totpositives == sum(r),
  tp = [0, tp, 1];
  fp = [0, fp, 1];
else
  tp = [0, tp];
  fp = [0, fp];
end

plot(fp, tp, varargin{:});

if max(fp, tp) < 1,
  tp = [tp, 1];
  fp = [fp, 1];
end

auc = trapz(fp, tp);
