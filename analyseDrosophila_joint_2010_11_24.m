% USER-specified: Sizes of the ranking sets for the first plot
T1 = [200, 400, 800, 1600, 3200, 6000];
% Sizes of the ranking sets for the second plot
T2 = [100, 200, 400, 800, 1600, 3200];

% USER-specified:: Directory where you store figures
ddir = 'figures/';
printPlot = 1; % 0 means not printing

% plot MAP models
plotMAP = [1 0 0 0 1 1]; 
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
scores = {}; I={}; sscores={}; J={};
r1 = zeros(length(posteriors), length(T1));
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
  end
end
h1 = figure;
set(gca, 'FontSize', fontSize);
% plot bars
h = bar(100*r1(plotMAP==1,:)');
set(gca, 'XTickLabel', T1);
hold on
v = axis;
v(3:4) = [0 100];
axis(v);
plot(v(1:2), 100*priors.accs31(1)*[1 1], 'b');
if incPrior,
  plot(v(1:2), 100*priors.accs31(2)*[1 1], 'g');
end
hold off
legends = {'MAP-32', 'MAP-32 + prior', 'MAP-16', 'MAP-16 + prior', 'Baseline', 'Regression', 'Uniform prior', 'Empirical prior', 'Location', 'EastOutside'};
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

% 1B PLOT ---------------------------------------------------------------
% GLOBAL RANKING BASED ON THE MAP MODEL, IN-SITU FILTERING
% START  ---------------------------------------------------------------

% Find the genes with positive in situ annotations
[C, IA, IB] = intersect(drosinsitu.genes(any(drosinsitu.data, 2)), post_genes);

scores = {}; I={}; sscores={}; J={};
r1 = zeros(length(posteriors), length(T1));
pvals1 = r1;
for k=1:length(posteriors)   
  if ((k <=2) | (k==6)), % 32 models 
    mycomb = combConstr;
  elseif  (k>2 & k<=4), % 16 models 
    mycomb = comb16;
  elseif (k==5) % 32 models
    mycomb = baselinecomb;
  end
  [scores{k}, I{k}] = max(posteriors{k}, [], 2);
  [sscores{k}, J{k}] = sort(scores{k}, 'descend');
  for l=1:length(T1),
    acc = 0; acc2 = 0;
    count = 0;
    for m=1:T1(l),
      if I{k}(J{k}(m)) ~= 1 && any(IB == J{k}(m)),
	      acc = acc + all( M(J{k}(m), mycomb( I{k}(J{k}(m)), :)==1), 2); 
          %acc2 = acc2 + all(M(J{k}(m), :) == mycomb(I{k}(J{k}(m)), :), 2);
	      count = count + 1;
      end
    end
    r1(k, l) = acc / count;
    %r2(k, l) = acc2 / count;
    pvals1(k, l) = 1 - binocdf(acc - 1, count, priors.focused_accs31(priors.accinds(k)));
  end
end
h1b = figure;
set(gca, 'FontSize', fontSize);
% plot bars
h = bar(100*r1(plotMAP==1,:)');
set(gca, 'XTickLabel', T1);
hold on
v = axis;
v(3:4) = [0 100];
axis(v);
plot(v(1:2), 100*priors.focused_accs31(1)*[1 1], 'b');
%plot(v(1:2), 100*prioraccs31(1)*[1 1], 'b--');
if incPrior,
  plot(v(1:2), 100*priors.focused_accs31(2)*[1 1], 'g');
  %plot(v(1:2), 100*prioraccs31(2)*[1 1], 'g--');
end
hold off
%legends = {'MAP-32', 'MAP-32 + prior', 'MAP-16', 'MAP-16 + prior', 'Baseline', 'Focused prior', 'Global prior', 'Focused Empirical prior', 'Global Empirical prior', 'Location', 'EastOutside'};
legends = {'MAP-32', 'MAP-32 + prior', 'MAP-16', 'MAP-16 + prior', 'Baseline', 'Regression', 'Uniform prior', 'Empirical prior', 'Location', 'EastOutside'};
%legend(legends([plotMAP, 1, 1, incPrior, incPrior]==1));
legend(legends([plotMAP, 1, incPrior]==1));
axis(v)
xlabel('# of global top genes')
ylabel('Enrichment (%)')
drosStarBars(h, pvals1(plotMAP==1,:)');
set(gca, 'FontSize', fontSize);
set(gcf, 'PaperUnits', 'centimeters')
set(gcf, 'PaperPosition', [0, 0, figSize])
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
    [BestPost, BestLink] = max(linkMargPosteriors{k}, [], 2);
    [foo, I] = sort(BestPost, 'descend');
    for l=1:length(T2),
       r2(l,k) = 0;  
       nanCnt =0;
       for j=1:T2(l)
          gene = I(j);   
          TF = BestLink(I(j));
          MM = M(gene, TF);
          if ~isnan(MM)
             r2(l,k) = r2(l,k) + MM;
          else
             nanCnt = nanCnt + 1;
          end
       end
       pvals2(l, k) = 1 - binocdf(r2(l,k) - 1, ...
			          T2(l)-nanCnt, ...
	  		          prioraccsSingleTF(baselines2(k)));
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
v(3:4) = [0 100];
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


baselines2 = [1, 2, 1, 2, 1. 2, 1, 1];
r3 = zeros(length(T2), length(linkPairPosteriors));
pvals3 = r3;
for k=1:length(linkPairPosteriors),
    % compute the best pair-TF-link for each gene
    [BestPost, BestLink] = max(linkPairPosteriors{k}, [], 2);
    [foo, I] = sort(BestPost, 'descend');
    for l=1:length(T2),
       r3(l,k) = 0;  
       nanCnt =0;
       for j=1:T2(l)
          gene = I(j);   
          TFpair = BestLink(I(j));
          MM =  prod(M(gene, pairs(TFpair,:)));
          if ~isnan(MM)
             r3(l,k) = r3(l,k) + MM;
          else
             nanCnt = nanCnt + 1;
          end
       end
       pvals3(l, k) = 1 - binocdf(r3(l,k) - 1, ...
			          T2(l)-nanCnt, ...
			          prioraccsPairTF(baselines2(1)));                  
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
v(3:4) = [0 100];
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
if printPlot 
   print(h1, '-depsc2', [ddir 'drosophilaBars_' 'EnrichmentGlobalMAP' num2str(2) '_', property '.eps']);
   print(h1b, '-depsc2', [ddir 'drosophilaBars_' 'EnrichmentGlobalMAPIN_SITU' num2str(2) '_', property '.eps']);
   print(h2, '-depsc2', [ddir 'drosophilaBars_' 'EnrichmentSingleLinks' num2str(sum(plotRest)) '_', property '.eps']); 
   print(h3, '-depsc2', [ddir 'drosophilaBars_' 'EnrichmentPairLinks' num2str(sum(plotRest)) '_', property '.eps']);
end
