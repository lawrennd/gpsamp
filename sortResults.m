function r = sortResults(r),

[r.genes, I] = sort(r.genes);
r.marlls = r.marlls(I, :, :);
