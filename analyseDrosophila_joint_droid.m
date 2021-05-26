% USER-specified:: Directory where you store figures
ddir = 'figures/';
validation = 'droid';
printPlot = 1; % 0 means not printing

% plot MAP models
plotMAP = [1 0 0 0 1 1 1];
%plotMAP = [1 1 0 0 0]; 
% for the rest plots 
plotRest = [1 0 0 0 1 0 1 1]; 
%plotRest = [1 1 0 0 1 1 0]; 

% USER-specified: whether to include empirical prior line
incPrior = 0;
figSize = [7 5];
fontSize = 7;

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


% 1 PLOT ---------------------------------------------------------------
% GLOBAL RANKING BASED ON THE MAP MODEL 
% START  ---------------------------------------------------------------
scores = {}; I={}; sscores={}; J={};
r1 = zeros(length(posteriors)+1, length(T1));
pvals1 = r1;
for k=1:length(posteriors)   
  if ((k<=2) | (k==6)), % 32 models 
    mycomb = combConstr;
  elseif  (k>2 & k<=4), % 16 models 
    mycomb = comb16;
  elseif (k==5)
    mycomb = baselinecomb;
  end
  [scores{k}, I{k}] = max(posteriors{k}, [], 2);
  [sscores{k}, J{k}] = sort(scores{k}, 'descend');
  for l=1:length(T1),
    acc = 0; acc2 = 0;
    count = 0;
    for m=1:T1(l),
      if I{k}(J{k}(m)) ~= 1,
	      acc = acc + all( M(J{k}(m), mycomb( I{k}(J{k}(m)), :)==1), 2); 
          %acc2 = acc2 + all(M(J{k}(m), :) == mycomb(I{k}(J{k}(m)), :), 2);
	      count = count + 1;
      end
    end
    r1(k, l) = acc / count;
    %r2(k, l) = acc2 / count;
    pvals1(k, l) = 1 - binocdf(acc - 1, count, priors.accs31(priors.accinds(k)));
    %pvals1(k, l) = 1 - hygecdf(acc-1, numel(M), sum(M(:)), count);
  end
end
Jinf = ~any(linkMargPosteriors{8} > 0, 2);
V = all(~((linkMargPosteriors{8} > 0) & ~M), 2);
r1(end, end) = mean(V(~Jinf));
pvals1(end, :) = 1;
pvals1(end, end) = 1 - binocdf(sum(V(~Jinf)) - 1, sum(~Jinf), priors.accs31(priors.accinds(k)));
h1 = figure;
set(gca, 'FontSize', fontSize);
% plot bars
h = bar(100*r1(plotMAP==1,:)');
set(gca, 'XTickLabel', T1);
hold on
v = axis;
v(3:4) = [0 80];
axis(v);
plot(v(1:2), 100*priors.accs31(1)*[1 1], 'b');
if incPrior,
  plot(v(1:2), 100*priors.accs31(2)*[1 1], 'g');
end
hold off
legends = {'MAP-32', 'MAP-32 + prior', 'MAP-16', 'MAP-16 + prior', ...
           'ML-Baseline', 'Regression', ...
           sprintf('Inferelator (only for %d genes)', T1(end)), ...
           'Uniform prior', 'Empirical prior', 'Location', 'EastOutside'};
legend(legends([plotMAP, 1, incPrior]==1));
axis(v)
xlabel('# of global top genes')
ylabel('Enrichment (%)')
drosStarBars(h, pvals1(plotMAP==1,:)');
set(gca, 'FontSize', fontSize);
set(gcf, 'PaperUnits', 'centimeters')
set(gcf, 'PaperPosition', [0, 0, figSize])
% 1 PLOT ---------------------------------------------------------------
% GLOBAL RANKING BASED ON THE MAP MODEL 
% END    ---------------------------------------------------------------

% 1C PLOT ---------------------------------------------------------------
% GLOBAL RANKING BASED ON THE MAP MODEL, INCLUDING NEGATIVE PREDICTIONS
% START  ---------------------------------------------------------------             
r1c = zeros(length(posteriors)+1, length(T1));
pvals1c = r1c;
for k=1:length(posteriors)   
  if ((k<=2) | (k==6)), % 32 models 
    mycomb = combConstr;
  elseif  (k>2 & k<=4), % 16 models 
    mycomb = comb16;
  elseif (k==5)
    mycomb = baselinecomb;
  end
  for l=1:length(T1),
    acc = 0; acc2 = 0;
    count = 0;
    for m=1:T1(l),
        %acc = acc + all( M(J{k}(m), mycomb( I{k}(J{k}(m)), :)==1), 2); 
        acc = acc + all(M(J{k}(m), :) == mycomb(I{k}(J{k}(m)), :), 2);
        count = count + 1;
    end
    r1c(k, l) = acc / count;
    %r2(k, l) = acc2 / count;
    pvals1c(k, l) = 1 - binocdf(acc - 1, count, priors.accs31_neg(priors.accinds(k)));
    %pvals1(k, l) = 1 - hygecdf(acc-1, numel(M), sum(M(:)), count);
  end
end
V = all((linkMargPosteriors{8} > 0) == M, 2);
r1c(end, end) = mean(V);
pvals1c(end, :) = 1;
pvals1c(end, end) = 1 - binocdf(sum(V) - 1, length(V), priors.accs31_neg(priors.accinds(k)));
h1c = figure;
set(gca, 'FontSize', fontSize);
% plot bars
h = bar(100*r1c(plotMAP==1,:)');
set(gca, 'XTickLabel', T1);
hold on
v = axis;
v(3:4) = [0 80];
axis(v);
plot(v(1:2), 100*priors.accs31_neg(1)*[1 1], 'b');
if incPrior,
  plot(v(1:2), 100*priors.accs31_neg(2)*[1 1], 'g');
end
hold off
legends = {'MAP-32', 'MAP-32 + prior', 'MAP-16', 'MAP-16 + prior', 'ML-Baseline', 'Regression', sprintf('Inferelator (only for %d genes)', T1(end)), 'Uniform prior', 'Empirical prior', 'Location', 'EastOutside'};
legend(legends([plotMAP, 1, incPrior]==1));
axis(v)
xlabel('# of top genes')
ylabel('Enrichment (%)')
drosStarBars(h, pvals1c(plotMAP==1,:)');
set(gca, 'FontSize', fontSize);
set(gcf, 'PaperUnits', 'centimeters')
set(gcf, 'PaperPosition', [0, 0, figSize])
% 1C PLOT ---------------------------------------------------------------
% GLOBAL RANKING BASED ON THE MAP MODEL, INCLUDING NEGATIVE PREDICTIONS
% END    ---------------------------------------------------------------

% 1B PLOT ---------------------------------------------------------------
% GLOBAL RANKING BASED ON THE MAP MODEL, IN-SITU FILTERING
% START  ---------------------------------------------------------------

% % Find the genes with positive in situ annotations
% [C, IA, IB] = intersect(drosinsitu.genes(any(drosinsitu.data, 2)), post_genes);

% scores = {}; I={}; sscores={}; J={};
% r1 = zeros(length(posteriors), length(T1));
% pvals1 = r1;
% for k=1:length(posteriors)   
%   if ((k <=2) | (k==6)), % 32 models 
%     mycomb = combConstr;
%   elseif  (k>2 & k<=4), % 16 models 
%     mycomb = comb16;
%   elseif (k==5) % 32 models
%     mycomb = baselinecomb;
%   end
%   [scores{k}, I{k}] = max(posteriors{k}, [], 2);
%   [sscores{k}, J{k}] = sort(scores{k}, 'descend');
%   for l=1:length(T1),
%     acc = 0; acc2 = 0;
%     count = 0;
%     for m=1:T1(l),
%       if I{k}(J{k}(m)) ~= 1 && any(IB == J{k}(m)),
% 	      acc = acc + all( M(J{k}(m), mycomb( I{k}(J{k}(m)), :)==1), 2); 
%           %acc2 = acc2 + all(M(J{k}(m), :) == mycomb(I{k}(J{k}(m)), :), 2);
% 	      count = count + 1;
%       end
%     end
%     r1(k, l) = acc / count;
%     %r2(k, l) = acc2 / count;
%     pvals1(k, l) = 1 - binocdf(acc - 1, count, priors.focused_accs31(priors.accinds(k)));
%   end
% end
% h1b = figure;
% set(gca, 'FontSize', fontSize);
% % plot bars
% h = bar(100*r1(plotMAP==1,:)');
% set(gca, 'XTickLabel', T1);
% hold on
% v = axis;
% v(3:4) = [0 100];
% axis(v);
% plot(v(1:2), 100*priors.focused_accs31(1)*[1 1], 'b');
% %plot(v(1:2), 100*prioraccs31(1)*[1 1], 'b--');
% if incPrior,
%   plot(v(1:2), 100*priors.focused_accs31(2)*[1 1], 'g');
%   %plot(v(1:2), 100*prioraccs31(2)*[1 1], 'g--');
% end
% hold off
% %legends = {'MAP-32', 'MAP-32 + prior', 'MAP-16', 'MAP-16 + prior', 'Baseline', 'Focused prior', 'Global prior', 'Focused Empirical prior', 'Global Empirical prior', 'Location', 'EastOutside'};
% legends = {'MAP-32', 'MAP-32 + prior', 'MAP-16', 'MAP-16 + prior', 'Baseline', 'Regression', 'Uniform prior', 'Empirical prior', 'Location', 'EastOutside'};
% %legend(legends([plotMAP, 1, 1, incPrior, incPrior]==1));
% legend(legends([plotMAP, 1, incPrior]==1));
% axis(v)
% xlabel('# of global top genes')
% ylabel('Enrichment (%)')
% drosStarBars(h, pvals1(plotMAP==1,:)');
% set(gca, 'FontSize', fontSize);
% set(gcf, 'PaperUnits', 'centimeters')
% set(gcf, 'PaperPosition', [0, 0, figSize])
% 1B PLOT ---------------------------------------------------------------
% GLOBAL RANKING BASED ON THE MAP MODEL, IN-SITU FILTERING
% END    ---------------------------------------------------------------


% 2 PLOT ---------------------------------------------------------------
% GLOBAL RANKING BASED ON PRESENCE OF SINGLE LINKS
% START  ---------------------------------------------------------------
%
% separate prior for eacn link being active 
%priorSingleTF = nanmean(M);
%prioraccsSingleTF(1) = mean(priorSingleTF);
%prioraccsSingleTF(2) = sum(priorSingleTF .* priorSingleTF) / sum(priorSingleTF).^2;
for k=1:numTFs
    indSingle = find(combConstr(:,k)==1);
    priorSingleTF(k) = sum(priors.prior32(indSingle));
end

% Bayes correct classification rate based on a classifier that  uses
% a unifrom prior (as classification rule)
prioraccsSingleTF(1) = mean(priorSingleTF);
% Bayes correct classification rate based on a classifier that uses 
% the empirical prior
prioraccsSingleTF(2) = sum(priorSingleTF .* priorSingleTF) / sum(priorSingleTF);

%prioraccsSingleTF(1) = mean(priorSingleTF);
%prioraccsSingleTF(2) = sum(priorSingleTF .* priorSingleTF) / sum(priorSingleTF).^2;
baselines2 = [1, 2, 1, 2, 1. 2, 1, 1];
numGenes = size(linkMargPosteriors{1},1);
r2 = zeros(length(T2), length(linkMargPosteriors));
pvals2 = r2;
for k=1:length(linkMargPosteriors), 
    % compute the best single-TF-link for each gene
    %[BestPost, BestLink] = max(, [], 2);
    [foo, I] = sort(linkMargPosteriors{k}(:), 'descend');
    for l=1:length(T2),
       r2(l,k) = 0;  
       nanCnt =0;
       for j=1:T2(l)
          MM = M(I(j));
          if ~isnan(MM)
             r2(l,k) = r2(l,k) + MM;
          else
             nanCnt = nanCnt + 1;
          end
       end
       % pvals2(l, k) = 1 - binocdf(r2(l,k) - 1, ...
       %  		          T2(l)-nanCnt, ...
       %    		          prioraccsSingleTF(baselines2(k)));
       pvals2(l, k) = 1 - hygecdf(r2(l,k)-1, numel(M), sum(M(:)), T2(l)-nanCnt);
       r2(l,k) = r2(l,k)/(T2(l)-nanCnt); 
    end
end
% plots bars 
h2 = figure;
set(gca, 'FontSize', fontSize);
h = bar(100*r2(:, plotRest==1));
set(gca, 'XTickLabel', T2);
hold on
v = axis;
v(3:4) = [0 80];
axis(v)
plot(v(1:2), 100*prioraccsSingleTF(1)*[1 1], 'b');
if incPrior,
  plot(v(1:2), 100*prioraccsSingleTF(2)*[1 1], 'g');
end
hold off
legends = {'Posterior-32', 'Posterior-32 + prior', 'Posterior-16', 'Posterior-16 + prior', ...
       'Posterior-2', 'Posterior-2 + prior',...
       'ML-Baseline', 'Inferelator', 'Uniform prior', 'Empirical prior'};
legend(legends([plotRest, 1, incPrior]==1));  
%legend('Posterior-32', 'Posterior-32 + prior', 'Posterior-16', 'Posterior-16 + prior', ...
%       'Posterior-2', 'Posterior-2 + prior',...
%       'Baseline', 'Uniform random', 'Random from prior', ...
%       'Location', 'EastOutside');
axis(v)
xlabel('# of top predictions')
ylabel('Enrichment (%)')
drosStarBars(h, pvals2(:, plotRest==1));
set(gca, 'FontSize', fontSize);
set(gcf, 'PaperUnits', 'centimeters')
set(gcf, 'PaperPosition', [0, 0, figSize])
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


pairM = zeros(size(M, 1), size(pairs, 1));
for j=1:size(pairs, 1),
  pairM(:, j) = all(M(:, pairs(j,:)),2);
end

baselines2 = [1, 2, 1, 2, 1. 2, 1, 1];
r3 = zeros(length(T2), length(linkPairPosteriors));
pvals3 = r3;
for k=1:length(linkPairPosteriors),
    % compute the best pair-TF-link for each gene
    [foo, I] = sort(linkPairPosteriors{k}(:), 'descend');
    for l=1:length(T2),
       r3(l,k) = 0;  
       nanCnt =0;
       for j=1:T2(l)
          MM =  pairM(I(j));
          if ~isnan(MM)
             r3(l,k) = r3(l,k) + MM;
          else
             nanCnt = nanCnt + 1;
          end
       end
       % pvals3(l, k) = 1 - binocdf(r3(l,k) - 1, ...
       %  		          T2(l)-nanCnt, ...
       %  		          prioraccsPairTF(baselines2(1)));
       pvals3(l, k) = 1 - hygecdf(r3(l,k)-1, numel(pairM), ...
                                  sum(pairM(:)), T2(l)-nanCnt);
       r3(l,k) = r3(l,k)/(T2(l)-nanCnt);
    end
    %[foo, I] = sort(linkPairPosteriors{k}(:), 'descend');
    %for l=1:length(T2),
    %   r3(l,k) = 0;  
    %   nanCnt = 0;
    %   for j=1:T2(l)
    %      gene = mod(I(j), numGenes);
    %      gene(gene==0)=numGenes;          
    %      TFpair = floor(I(j)/numGenes) + 1;
    %      TFpair(TFpair==11)=10;
    % 
    %      MM =  prod(M(gene, pairs(TFpair,:)));
    %      if ~isnan(MM)
    %         r3(l,k) = r3(l,k) + MM;
    %      else
    %         nanCnt = nanCnt + 1;
    %      end
    %   end
    %   pvals3(l, k) = 1 - binocdf(r3(l,k) - 1, ...
	%		          T2(l)-nanCnt, ...
	%		          prioraccsPairTF(baselines2(1)));                  
    %   r3(l,k) = r3(l,k)/(T2(l)-nanCnt);
    %end
%
end
% plots bars 
h3 = figure;
set(gca, 'FontSize', fontSize);
h = bar(100*r3(:, plotRest==1));
set(gca, 'XTickLabel', T2);
hold on
v = axis;
v(3:4) = [0 20];
axis(v)
plot(v(1:2), 100*prioraccsPairTF(1)*[1 1], 'b');
if incPrior,
  plot(v(1:2), 100*prioraccsPairTF(2)*[1 1], 'g');
end
hold off
legends = {'Posterior-32', 'Posterior-32 + prior', 'Posterior-16', 'Posterior-16 + prior', 'Posterior-4','Posterior-4 + prior',...
       'ML-Baseline', 'Inferelator', 'Uniform prior', 'Empirical prior'};
legend(legends([plotRest, 1, incPrior]==1));  
%legend('Posterior-32', 'Posterior-32 + prior', 'Posterior-16', 'Posterior-16 + prior', 'Posterior-4','Posterior-4 + prior',...
%       'Baseline', 'Uniform random', 'Random from prior', ...
%       'Location', 'EastOutside');
axis(v)
xlabel('# of top predictions')
ylabel('Enrichment (%)')
drosStarBars(h, pvals3(:, plotRest==1));
set(gca, 'FontSize', fontSize);
set(gcf, 'PaperUnits', 'centimeters')
set(gcf, 'PaperPosition', [0, 0, figSize])
% 3 PLOT ---------------------------------------------------------------
% GLOBAL RANKING BASED ON THE PRESENCE OF A PAIR OF  LINKS
% END    ---------------------------------------------------------------


% print Plots
property = 'Constrained';
if flag ~= 1
    property = 'Unconstrained'; 
end

switch validation,
  case 'droid',
    namestem = 'drosophilaDroidBars_';
  case 'chipchip',
    namestem = 'drosophilaBars_';
  otherwise,
    error('unknown validation method');
end

if printPlot 
   print(h1, '-depsc2', [ddir namestem 'EnrichmentGlobalMAP' num2str(2) '_', property '.eps']);
   %print(h1b, '-depsc2', [ddir namestem 'EnrichmentGlobalMAPIN_SITU' num2str(2) '_', property '.eps']);
   print(h1c, '-depsc2', [ddir namestem 'EnrichmentGlobalMAPNeg' num2str(2) '_', property '.eps']);
   print(h2, '-depsc2', [ddir namestem 'EnrichmentSingleLinks' num2str(sum(plotRest)) '_', property '.eps']); 
   print(h3, '-depsc2', [ddir namestem 'EnrichmentPairLinks' num2str(sum(plotRest)) '_', property '.eps']);
end
