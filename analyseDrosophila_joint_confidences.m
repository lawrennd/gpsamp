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

if (~exist('r1rand')),
  load results/randomisation;
end


% 1 PLOT ---------------------------------------------------------------
% GLOBAL RANKING BASED ON THE MAP MODEL 
% START  ---------------------------------------------------------------

pvals1 = p_from_rand(r1, r1rand);
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

pvals1b = p_from_rand(r1b, r1brand);
h1b = figure;
set(gca, 'FontSize', fontSize);
% plot bars
h = bar(100*r1b(plotMAP==1,:)');
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
drosStarBars(h, pvals1b(plotMAP==1,:)');
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

pvals2 = p_from_rand(r2, r2rand);
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

pvals3 = p_from_rand(r3, r3rand);
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
   print(h1, '-depsc2', [ddir 'drosophilaBars_' 'EnrichmentGlobalMAP' num2str(2) '_', property '_newconf.eps']);
   print(h1b, '-depsc2', [ddir 'drosophilaBars_' 'EnrichmentGlobalMAPIN_SITU' num2str(2) '_', property '_newconf.eps']);
   print(h2, '-depsc2', [ddir 'drosophilaBars_' 'EnrichmentSingleLinks' num2str(sum(plotRest)) '_', property '_newconf.eps']); 
   print(h3, '-depsc2', [ddir 'drosophilaBars_' 'EnrichmentPairLinks' num2str(sum(plotRest)) '_', property '_newconf.eps']);
end
