% USER-specified: Sizes of the ranking sets for the first plot
T1 = [200, 400, 800, 1600, 3200, 6000];
% Sizes of the ranking sets for the second plot
T2 = [100, 200, 400, 800, 1600, 3200];

N_REPEATS = 100000;

analyseDrosophila_constants;
if (~exist('posteriors')),
  [posteriors, linkMargPosteriors, linkPairPosteriors, priors, M, post_genes] = analyseDrosophila_loadResults;
else
  fprintf('Re-using old posteriors, clear ''posteriors'' to force reload/recompute\n');
end

if (~exist('drosinsitu')),
  load datasets/drosophila_data;
end

% 1 PLOT ---------------------------------------------------------------
% GLOBAL RANKING BASED ON THE MAP MODEL 
% START  ---------------------------------------------------------------             
r1 = zeros(length(posteriors), length(T1));
r1rand = zeros(length(posteriors), length(T1), N_REPEATS);
for k=1:length(posteriors)   
  if ((k<=2) || (k==6)), % 32 models 
    mycomb = combConstr;
  elseif  (k>2 && k<=4), % 16 models 
    mycomb = comb16;
  elseif (k==5)
    mycomb = baselinecomb;
  end
  myM = zeros(size(M, 1), size(mycomb, 1));
  for j=1:size(mycomb, 1),
    myM(:, j) = all(M(:, mycomb(j,:)==1),2);
  end
  [scores, I] = max(posteriors{k}, [], 2);
  [sscores, J] = sort(scores, 'descend');
  NMAX = max(T1);
  hits = myM(sub2ind(size(myM), J(1:NMAX), I(J(1:NMAX))));
  hits(I(J(1:NMAX))==1) = NaN;
  for l=1:length(T1),
    r1(k, l) = nanmean(hits(1:T1(l)));
  end
  for n=1:N_REPEATS,
    if (rem(n, 1000)==0),
      fprintf('Task 1, k=%d/%d, repeat #%d/%d\n', k, length(posteriors), ...
              n, N_REPEATS);
    end
    R = randperm(size(M, 1));
    hits = myM(sub2ind(size(myM), J(R(1:NMAX)), I(J(1:NMAX))));
    hits(I(J(1:NMAX))==1) = NaN;
    for l=1:length(T1),
      r1rand(k, l, n) = nanmean(hits(1:T1(l)));
    end
  end
end
% 1 PLOT ---------------------------------------------------------------
% GLOBAL RANKING BASED ON THE MAP MODEL 
% END    ---------------------------------------------------------------

% 1B PLOT ---------------------------------------------------------------
% GLOBAL RANKING BASED ON THE MAP MODEL, IN-SITU FILTERING
% START  ---------------------------------------------------------------

% Find the genes with positive in situ annotations
[C, IA, IB] = intersect(drosinsitu.genes(any(drosinsitu.data, 2)), post_genes);

isIB = zeros(length(post_genes), 1);
isIB(IB) = 1;
isIB = logical(isIB);
IBind = cumsum(isIB);
IBind(IBind==0) = 1;

r1b = zeros(length(posteriors), length(T1));
r1brand = zeros(length(posteriors), length(T1), N_REPEATS);
for k=1:length(posteriors)   
  if ((k <=2) || (k==6)), % 32 models 
    mycomb = combConstr;
  elseif  (k>2 && k<=4), % 16 models 
    mycomb = comb16;
  elseif (k==5) % 32 models
    mycomb = baselinecomb;
  end
  myM = zeros(size(M, 1), size(mycomb, 1));
  for j=1:size(mycomb, 1),
    myM(:, j) = all(M(:, mycomb(j,:)==1),2);
  end
  [scores, I] = max(posteriors{k}, [], 2);
  [sscores, J] = sort(scores, 'descend');
  NMAX = max(T1);
  hits = myM(sub2ind(size(myM), J(1:NMAX), I(J(1:NMAX))));
  hits(I(J(1:NMAX))==1) = NaN;
  hits(~isIB(J(1:NMAX))) = NaN;
  for l=1:length(T1),
    r1b(k, l) = nanmean(hits(1:T1(l)));
  end
  for n=1:N_REPEATS,
    if (rem(n, 1000)==0),
      fprintf('Task 2, k=%d/%d, repeat #%d/%d\n', k, length(posteriors), ...
              n, N_REPEATS);
    end
    R = randperm(length(IB));
    hits = myM(sub2ind(size(myM), IB(R(IBind(J(1:NMAX)))), I(J(1:NMAX))));
    hits(I(J(1:NMAX))==1) = NaN;
    hits(~isIB(J(1:NMAX))) = NaN;
    % for m=1:NMAX,
    %   if I(J(m)) ~= 1 && isIB(J(m)),
    %     t = IBind(J(m));
    %     hits(m) = myM(IB(R(t)), I(J(m)));
    %   else
    %     hits(m) = NaN;
    %   end
    % end
    for l=1:length(T1),
      r1brand(k, l, n) = nanmean(hits(1:T1(l)));
    end
  end
end
% 1B PLOT ---------------------------------------------------------------
% GLOBAL RANKING BASED ON THE MAP MODEL, IN-SITU FILTERING
% END    ---------------------------------------------------------------

% 2 PLOT ---------------------------------------------------------------
% GLOBAL RANKING BASED ON PRESENCE OF SINGLE LINKS
% START  ---------------------------------------------------------------
%

numGenes = size(linkMargPosteriors{1},1);
r2 = zeros(length(T2), length(linkMargPosteriors));
r2rand = zeros(length(T2), length(linkMargPosteriors), N_REPEATS);
for k=1:length(linkMargPosteriors), 
    % compute the best single-TF-link for each gene
    [BestPost, BestLink] = max(linkMargPosteriors{k}, [], 2);
    [foo, I] = sort(BestPost, 'descend');
    NMAX = max(T2);
    hits = M(sub2ind(size(M), I(1:NMAX), BestLink(I(1:NMAX))));
    for l=1:length(T2),
      r2(l, k) = mean(hits(1:T2(l)));
    end
    for n=1:N_REPEATS,
      if (rem(n, 1000)==0),
        fprintf('Task 2, k=%d/%d, repeat #%d/%d\n', k, ...
                length(linkMargPosteriors), n, N_REPEATS);
      end
      J = randperm(size(M, 1));
      hits = M(sub2ind(size(M), I(J(1:NMAX)), BestLink(I(1:NMAX))));
      for l=1:length(T2),
	r2rand(l, k, n) = mean(hits(1:T2(l)));
      end
    end
end

% 2 PLOT ---------------------------------------------------------------
% GLOBAL RANKING BASED ON PRESENCE OF SINGLE LINKS
% END    ---------------------------------------------------------------


% 3 PLOT ---------------------------------------------------------------
% GLOBAL RANKING BASED ON THE PRESENCE OF A PAIR OF  LINKS
% START  ---------------------------------------------------------------
%
% separate prior for each pairs of links being active 
cnt = 0;
for k=1:numTFs
  for g=(k+1):numTFs
      cnt = cnt + 1;
      indPair = find(combConstr(:,k)==1 & combConstr(:,g)==1);
      priorPairTF(cnt) = sum(priors.prior32(indPair));
  end
end

% Bayes correct classification rate based on a classifier that  uses
% a unifrom prior (as classification rule)
prioraccsPairTF(1) = mean(priorPairTF);
% Bayes correct classification rate based on a classifier that uses 
% the empirical prior
prioraccsPairTF(2) = sum(priorPairTF.*priorPairTF) / sum(priorPairTF);


r3 = zeros(length(T2), length(linkPairPosteriors));
r3rand = zeros(length(T2), length(linkPairPosteriors), N_REPEATS);
for k=1:length(linkPairPosteriors),
    % compute the best pair-TF-link for each gene
    [BestPost, BestLink] = max(linkPairPosteriors{k}, [], 2);
    [foo, I] = sort(BestPost, 'descend');
    NMAX = max(T2);
    hits = prod([M(sub2ind(size(M), I(1:NMAX), ...
                           pairs(BestLink(I(1:NMAX))', 1))), ...
		 M(sub2ind(size(M), I(1:NMAX), ...
                           pairs(BestLink(I(1:NMAX))', 2)))], 2);
    for l=1:length(T2),
      r3(l, k) = mean(hits(1:T2(l)));
    end
    for n=1:N_REPEATS,
      if (rem(n, 1000)==0),
        fprintf('Task 3, k=%d/%d, repeat #%d/%d\n', k, ...
                length(linkPairPosteriors), n, N_REPEATS);
      end
      J = randperm(size(M, 1));
      hits = prod([M(sub2ind(size(M), I(J(1:NMAX)), ...
                             pairs(BestLink(I(1:NMAX))', 1))), ...
                   M(sub2ind(size(M), I(J(1:NMAX)), ...
                             pairs(BestLink(I(1:NMAX))', 2)))], 2);
      for l=1:length(T2),
	r3rand(l, k, n) = mean(hits(1:T2(l)));
      end
    end
%
end
% 3 PLOT ---------------------------------------------------------------
% GLOBAL RANKING BASED ON THE PRESENCE OF A PAIR OF  LINKS
% END    ---------------------------------------------------------------

% save results/randomisation r1 r1rand r1b r1brand r2 r2rand r3 r3rand
