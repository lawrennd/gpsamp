function drosProcessAll5(filename, N, demo),

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
  marlls{k} = gpmtfSummariseResults(r.testGene, 'margLogLik2', Genes, GenesVar, TFs, models);
end

results.genes = cat(1, genes{:});
results.marlls = cat(1, marlls{:});

save([filename '_summary2.mat'], '-struct', 'results');
