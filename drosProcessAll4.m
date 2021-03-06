function drosProcessAll4(filename, N, demo),

genes = {};
marlls = {};

for k=1:N,
  fname = sprintf('%s_m%d_r%d.mat', filename, N, k);
  try,
    r = load(fname);
  catch
    warning('error reading file %s, skipping', fname)
    continue
  end
  fprintf('Processing file %s, %d/%d...\n', fname, k, N);
  [Genes, GenesVar, TFs, models] = feval(demo, N, k, '', 0);
  genes{k} = r.mygenes(1:size(r.testGene, 1));
  marlls{k} = gpmtfSummariseResults(r.testGene, 'margLogLik1', Genes, GenesVar, TFs, models);
end

if iscell(genes{1}),
  results.genes = cat(1, genes{:});
else
  results.genes = cat(2, genes{:});
end
results.marlls = cat(1, marlls{:});

save([filename '_summary.mat'], '-struct', 'results');
