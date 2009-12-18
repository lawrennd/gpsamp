N = 16;

genes = {};
pplrs = {};
lls = {};

for k=1:N,
  fname = sprintf('multitf_2009-12-17_m%d_r%d.mat', N, k);
  r = load(fname);
  genes{k} = r.mygenes;
  pplrs{k} = gpmtfSummariseResults(r.testGene, 'pplr2');
  zscores{k} = gpmtfSummariseResults(r.testGene, 'zscore');
  lls{k} = gpmtfSummariseResults(r.testGene, 'loglike');
end

results.genes = cat(1, genes{:});
results.pplrs = cat(2, pplrs{:});
results.zscores = cat(2, zscores{:});
results.lls = cat(2, lls{:});

% save multitf_2009-12-17_summary1 -struct results
