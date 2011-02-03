function drosProcess32Figures(filename1, N, demo, filestem, my_k),

genes = {};
marlls = {};

if nargin < 4,
  filestem = 'dros';
end

if nargin < 5,
  my_k = 1:N;
end

for k=my_k,
  fname1 = sprintf('%s_m%d_r%d.mat', filename1, N, k);
  fprintf('Loading results from %s\n', fname1);
  try,
    r = load(fname1);
  catch
    warning('error reading file %s, skipping', fname1)
    continue
  end
  fprintf('Processing file %s, %d/%d...\n', fname1, k, N);
  [Genes, GenesVar, TFs, models, mygenes] = feval(demo, N, k, '', 0);

  for l=1:length(mygenes),
    close all;
    gpmtfTestPlot16Models(r.testGene(l, :), ...
			  Genes(l, :, :), GenesVar(l, :, :), ...
			  TFs, models, mygenes{l}, filestem, 1);
  end
end
