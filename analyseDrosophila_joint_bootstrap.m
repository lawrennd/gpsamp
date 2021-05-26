N_REPEATS = 100000;
validation = 'droid';

plotMAP = [1 0 0 0 1 1 1];
% for the rest plots 
plotRest = [1 0 0 0 1 0 1 1]; 

analyseDrosophila_constants;
if (~exist('posteriors')),
  [posteriors, linkMargPosteriors, linkPairPosteriors, priors, M, post_genes] = analyseDrosophila_loadResults(validation);
else
  fprintf('Re-using old posteriors, clear ''posteriors'' to force reload/recompute\n');
end

T1(end) = size(M, 1);

if (~exist('drosinsitu')),
  load datasets/drosophila_data;
end

s = RandStream.create(RandStream.getDefaultStream.Type, 'seed', 5489);
RandStream.setDefaultStream(s);

% 1 PLOT ---------------------------------------------------------------
% GLOBAL RANKING BASED ON THE MAP MODEL 
% START  ---------------------------------------------------------------

fprintf('Starting task 1/4...\n');

mycomb = combConstr;
myM = zeros(size(M, 1), size(mycomb, 1));
for j=2:size(mycomb, 1),
    myM(:, j) = all(M(:, mycomb(j,:)==1),2);
end
myM(:, 1) = NaN;

rr1 = cell(size(posteriors));
for k=1:length(posteriors)   
  [scores, I] = max(posteriors{k}, [], 2);
  [sscores, J] = sort(scores, 'descend');
  rr1{k} = [J, I(J)];
end
V = (linkMargPosteriors{8} > 0);
[C, I, J] = unique([mycomb; V], 'rows', 'first');
% First elements of J give the new sorting of mycomb
assert(all(all(mycomb == C(J(1:32), :))));
% Remaining elements of I(J) give the rows of mycomb
% corresponding to V
assert(all(all(mycomb(I(J(33:end)), :) == V)));
rr1{k+1} = [(1:size(V, 1))', I(J(33:end))];

[res1, accs1] = analyseDrosophila_runBootstrap(rr1(plotMAP~=0), myM, N_REPEATS, T1);

% 1 PLOT ---------------------------------------------------------------
% GLOBAL RANKING BASED ON THE MAP MODEL 
% END    ---------------------------------------------------------------

% 1C PLOT --------------------------------------------------------------
% GLOBAL RANKING BASED ON THE MAP MODEL WITH NEGATIVE VALIDATION
% START  ---------------------------------------------------------------

fprintf('Starting task 2/4...\n');

mycomb = combConstr;
myM2 = zeros(size(M, 1), size(mycomb, 1));
for j=1:size(mycomb, 1),
    myM2(:, j) = all(bsxfun(@eq, M, mycomb(j,:)),2);
end

[res1c, accs1c] = analyseDrosophila_runBootstrap(rr1(plotMAP~=0), myM2, N_REPEATS, T1);

% 1C PLOT --------------------------------------------------------------
% GLOBAL RANKING BASED ON THE MAP MODEL WITH NEGATIVE VALIDATION
% END    ---------------------------------------------------------------


% 2 PLOT ---------------------------------------------------------------
% GLOBAL RANKING BASED ON PRESENCE OF SINGLE LINKS
% START  ---------------------------------------------------------------
%

fprintf('Starting task 3/4...\n');

rr2 = cell(size(linkMargPosteriors));
for k=1:length(linkMargPosteriors), 
    % compute the best single-TF-link for each gene
    [foo, v] = sort(linkMargPosteriors{k}(:), 'descend');
    [v1, v2] = ind2sub(size(M), v);
    rr2{k} = [v1, v2];
end

[res2, accs2] = analyseDrosophila_runBootstrap(rr2(plotRest~=0), M, N_REPEATS, T2);

% 2 PLOT ---------------------------------------------------------------
% GLOBAL RANKING BASED ON PRESENCE OF SINGLE LINKS
% END    ---------------------------------------------------------------


% 3 PLOT ---------------------------------------------------------------
% GLOBAL RANKING BASED ON THE PRESENCE OF A PAIR OF  LINKS
% START  ---------------------------------------------------------------

fprintf('Starting task 4/4...\n');

pairM = zeros(size(M, 1), size(pairs, 1));
for j=1:size(pairs, 1),
  pairM(:, j) = all(M(:, pairs(j,:)),2);
end

rr3 = cell(size(linkPairPosteriors));
for k=1:length(linkPairPosteriors), 
    % compute the best single-TF-link for each gene
    [foo, v] = sort(linkPairPosteriors{k}(:), 'descend');
    [v1, v2] = ind2sub(size(M), v);
    rr3{k} = [v1, v2];
end

[res3, accs3] = analyseDrosophila_runBootstrap(rr3(plotRest~=0), pairM, N_REPEATS, T2);
% 3 PLOT ---------------------------------------------------------------
% GLOBAL RANKING BASED ON THE PRESENCE OF A PAIR OF  LINKS
% END    ---------------------------------------------------------------

switch validation,
  case 'droid',
    save results/bootstrap_droid res1 res1c res2 res3 accs1 accs1c accs2 accs3
  case 'chipchip',
    save results/bootstrap_chipchip res1 res1c res2 res3 accs1 accs1c accs2 accs3
  otherwise,
    save results/bootstrap_unknown res1 res1c res2 res3 accs1 accs1c accs2 accs3
end
