function M = drosMakeValidationMatrix(chipdistances, genes, threshold),

tforder = {'tin', 'bin', 'twi', 'bap', 'mef2'};

if nargin < 3,
  threshold = 2000;
end

I = zeros(size(tforder));
for k=1:length(tforder),
  I(k) = find(~cellfun(@isempty, strfind(chipdistances.labels, tforder{k})));
end

J = drosFindGeneinds(chipdistances, genes, 1);

M = zeros(length(J), length(I));
M(J==0, :) = NaN;
M(J~=0, :) = chipdistances.data(J(J~=0), I) < threshold;
