function drosProcessAll4(filename, N),

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
  genes{k} = r.mygenes(1:size(r.testGene, 1));
  marlls{k} = gpmtfSummariseResults(r.testGene, 'margLogLik1');
end

results.genes = cat(1, genes{:});
results.marlls = cat(1, marlls{:});

save([filename '_summary.mat'], '-struct', 'results');
