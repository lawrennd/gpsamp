function [res, accs] = analyseDrosophila_runBootstrap(rankings, validation, repeats, T)

% DROSBOOTSTRAPEVALUATION Performs bootstrap sampling of rankings
% FORMAT
% DESC Performs bootstrap sampling of rankings to assess significance of differences.
% ARG rankings : a cell array of R rankings, each of which is an array
% of indices of genes in drosexp in order of descending preference
% ARG validation : a binary vector of size(drosexp.genes) indicating
% the validation results of each gene, or NaN if the gene should be ignored
% ARG repeats : the number of repeats N
% ARG threshold : a vector of T n's in top-n to consider
% RETURN res : And RxRxT array where element (i,j,t) is the frequency that
% method i was better than method j for t'th threshold
% RETURN accs : An NxRxT array where element (n,i,t) is the accuracy of
% method i on repeat n for t'th threshold
%
% SEEALSO : drosLoadData, demRunRankings, drosDoBootstrap
%
% COPYRIGHT : Antti Honkela, 2009-2011

res = zeros([length(rankings), length(rankings), length(T)]);
accs = zeros([repeats, length(rankings), length(T)]);

N = size(validation, 1);
val = validation;
r = rankings;
r2 = cell(size(r));
for k=1:length(r),
    r2{k} = sub2ind(size(val), r{k}(:, 1), r{k}(:, 2));
end

if any(any(isnan(val))),
    for j=1:length(T),
        fprintf('Running iteration %d/%d\n', j, length(T));
        for n=1:repeats,
            I = randi(N, N, 1);
            inds = accumarray(I, 1, [N, 1]);

            myacc = zeros(1, length(rankings));
            for k=1:length(rankings),
                J = cumsum(inds(r{k}(:, 1)));
                t = find(J >= T(j), 1);
                myinds = zeros(size(r{k}(:, 1)));
                myinds(1:t-1) = inds(r{k}(1:t-1, 1));
                myinds(t) = inds(r{k}(t, 1)) - J(t) + T(j);
                myacc(k) = nansum(myinds(myinds ~= 0) .* val(r2{k}(myinds ~= 0))) ./ sum(myinds .* ~isnan(val(r2{k})));
            end
            accs(n, :, j) = myacc;
        end

        for k=1:length(rankings),
            for l=k+1:length(rankings),
                res(k, l, j) = mean(accs(:, k, j) >= accs(:, l, j));
                res(l, k, j) = mean(accs(:, k, j) < accs(:, l, j));
            end
        end
    end
else
    for j=1:length(T),
        fprintf('Running iteration %d/%d\n', j, length(T));
        for n=1:repeats,
            I = randi(N, N, 1);
            inds = accumarray(I, 1, [N, 1]);

            myacc = zeros(1, length(rankings));
            for k=1:length(rankings),
                J = cumsum(inds(r{k}(:, 1)));
                t = find(J >= T(j), 1);
                myinds = zeros(size(r{k}(:, 1)));
                myinds(1:t-1) = inds(r{k}(1:t-1, 1));
                myinds(t) = inds(r{k}(t, 1)) - J(t) + T(j);
                myacc(k) = sum(myinds(myinds ~= 0) .* val(r2{k}(myinds ~= 0))) ./ sum(myinds);
            end
            accs(n, :, j) = myacc;
        end

        for k=1:length(rankings),
            for l=k+1:length(rankings),
                res(k, l, j) = mean(accs(:, k, j) >= accs(:, l, j));
                res(l, k, j) = mean(accs(:, k, j) < accs(:, l, j));
            end
        end
    end
end
