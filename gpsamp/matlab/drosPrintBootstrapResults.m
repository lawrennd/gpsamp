function drosPrintBootstrapResults(m, values, labels, sizes, file, ignore),

if nargin < 5,
    file = [];
end

if nargin < 6,
    ignore = [];
end

N = size(m, 3);

if isempty(file),
    fid = 1;
else
    fid = fopen(file, 'w');
end

print_it(fid, m(values, values, :), sizes, labels, ignore);

if fid ~=1,
    fclose(fid);
end


function print_it(fid, m, mysize, labels, ignore),

[M, N, O] = size(m);

LABELS = labels;

for k=1:length(LABELS)-1,
    fprintf(fid, ' &');
end
fprintf(fid, ' \\\\\n');

for k=1:O-1,
    fprintf(fid, '\\multicolumn{%d}{c|}{Top %d} &\n', length(labels)+1, mysize(k));
end
fprintf(fid, '\\multicolumn{%d}{c}{Top %d}\\\\\n', length(labels)+1, mysize(O));

for k=1:O-1,
    fprintf(fid, '    &');
    fprintf(fid, ' %s &', LABELS{:});
    fprintf(fid, '\n');
end
fprintf(fid, '    &');
fprintf(fid, ' %s &', LABELS{1:end-1});
fprintf(fid, ' %s \\\\\n', LABELS{end});

for k=1:M,
  for j=1:O,
    fprintf(fid, '%s', LABELS{k});
    for l=1:N,
      if any(k == ignore) || any(l == ignore),
	fprintf(fid, ' & -  ');
      elseif m(k,l,j) > .99,
	fprintf(fid, ' & ***');
      elseif m(k,l,j) > .95,
	fprintf(fid, ' & ** ');
      elseif m(k,l,j) > .9,
	fprintf(fid, ' & *  ');
      elseif m(k,l,j) > .8,
	fprintf(fid, ' & +  ');
      elseif m(k,l,j) > .7,
	fprintf(fid, ' & .  ');
      else
	fprintf(fid, ' &    ');
      end
    end
    if j<O,
      fprintf(fid, ' & ');
    end
    fprintf(fid, '\n');
  end
  fprintf(fid, '\\\\\n');
end




% drosPrintBootstrapMatrices(matrices{1}, 1:4)
% drosPrintBootstrapMatrices(matrices{2}, 1:4, 4)
% drosPrintBootstrapMatrices(matrices{3}, 1:4)
% drosPrintBootstrapMatrices(matrices{4}, 1:4)
% drosPrintBootstrapMatrices(matrices{5}, 1:4, 4)
% drosPrintBootstrapMatrices(matrices{6}, 1:4)
