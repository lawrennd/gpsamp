function pvals = p_from_rand(r, rrand)

[d1, d2] = size(r);
d3 = size(rrand, 3);

pvals = zeros(size(r));

for i=1:d1,
    for j=1:d2,
        v = sort(squeeze(rrand(i, j, :)), 'descend');
        pvals(i, j) = min(find(r(i,j) > v)) / d3;
    end
end
