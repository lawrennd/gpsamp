results{1} = sortResults(load('results/multitfToyOneCond6a_2010-07-09_summary.mat'));
results{2} = sortResults(load('results/multitfToyTwoConds6a_2010-07-09_summary.mat'));
load datasets/toy4TFs28_June_10.mat

tcomb = [0 0 0; 1 0 0; 0 1 0; 0 0 1;
	 1 1 0; 1 0 1; 0 1 1;
	 1 1 1];
gstyles = {'m--', 'g--'};
styles = {'r-.', 'b'};
styles2 = {'r*-.', 'bx-'};
FONTSIZE=6;


if ~exist('resultsGroundTr'),
  % compute optimal likelihoods/probalilities using the ground truth parameters 
  Grtr_kineticsTF = modelTest.Likelihood.kineticsTF; 
  testsetIndices = 31:1:1030;
  Genes = Genes(testsetIndices, :, :); 
  GenesVar = GenesVar(testsetIndices, :, :);     
  % loop for the different conditions 
  for k=1:2
    if k == 1 
      TestGenes = Genes(:,:,1);
      TestGenesVar = GenesVar(:,:,1);
    else
      TestGenes = Genes;
      TestGenesVar = GenesVar;
    end
    % loop for the different models 
    for m=1:size(tcomb)
    %
        TFset = find(tcomb(m,:));
        numTFs = sum(tcomb(m,:));
        % model options
        options = gpmtfOptions(ones(1,size(TimesG,2),k), numTFs);
        options.jointAct = 'sigmoid';
        noiseM = {'pumaWhite'};
        %options.spikePriorW = 'yes';
        options.noiseModel = noiseM;
        %options.constraints.spaceW = 'positive';
        options.tauMax = 0; % no delays
        % define the dense discretized grid in the time axis for the TF latent functions
        [options, TimesF] = gpmtfDiscretize(TimesG, options);
     
        % loop over genes 
        for j=1:size(TestGenes,1)
            % compute the log likelihood            
            tmp = [];
            if numTFs > 0
               modelTest = gpmtfCreate(TestGenes(j,:,:), TestGenesVar(j,:,:), [], [], TimesG, TimesF, options);         
               modelTest.Likelihood.kinetics = GrTruth{j+30}.kinetics; 
               modelTest.Likelihood.kineticsTF = Grtr_kineticsTF(TFset,:);
               for r=1:modelTest.Likelihood.numReplicas
                  modelTest.Likelihood.TF(:,:,r) = TFs(TFset,:,r);
               end
               modelTest.Likelihood.W = GrTruth{j+30}.W(TFset);
               modelTest.Likelihood.W0 = GrTruth{j+30}.W0;          
               %modelTest.Likelihood
               %[j m]
               %pause
               for r=1:modelTest.Likelihood.numReplicas
                  tmp(r) = gpmtfLogLikelihoodGene(modelTest.Likelihood, zeros(size(modelTest.Likelihood.TF)), r, 1);
               end
            %   
            else 
               % only decay model
               B = GrTruth{j+30}.kinetics(1);
               D = GrTruth{j+30}.kinetics(2);
               A = GrTruth{j+30}.kinetics(4);
               PredGenes = B/D  + (A - B/D)*exp(-TimesG*D);
               for r=1:k
                   tmp(r) = - 0.5*sum(log(2*pi*TestGenesVar(1,:,r) ),2)....
                            - 0.5*sum(((TestGenes(j,:,r) - PredGenes).^2)./TestGenesVar(1,:,r),2);
               end  
            end
            resultsGroundTr{k}.marlls(j,m) = sum(tmp);  
        end
        %    
    end
  end
end




auc = {};
M = Net(results{1}.genes, 1:3);
linkprobs = {};
I = {};
val = {};
val1 = {};
figure(1);
for k=1:2,
  %subplot(1, 2, k);
  for l=1:3,
    linkprobs{k}(:, l) = logsumexp(results{k}.marlls(:, tcomb(:, l)==1), 2) - ...
	logsumexp(results{k}.marlls(:, tcomb(:, l)==0), 2);
  end
  for l=1:3,
    trlinkprobs{k}(:, l) = logsumexp(resultsGroundTr{k}.marlls(:, tcomb(:, l)==1), 2) - ...
	logsumexp(resultsGroundTr{k}.marlls(:, tcomb(:, l)==0), 2);
  end
  [foo, I{k}] = sort(linkprobs{k}(:), 'descend');
  val{k} = M(I{k});
  auc{1}(k) = drosPlotROC(val{k}', sum(val{k}), length(val{k}), styles{k});
  hold on
  [foo, I{k}] = sort(trlinkprobs{k}(:), 'descend');
  val{k} = M(I{k});
  auc{2}(k) = drosPlotROC(val{k}', sum(val{k}), length(val{k}), gstyles{k});
  %p = linkprobs{k}(1:500, :);
  %[foo, I{k}] = sort(p(:), 'descend');
  %N = M(1:500, :);
  %val1{k} = N(I{k});
  %auc{1}(k) = drosPlotROC(val1{k}', sum(val1{k}), length(val1{k}), styles{k});
  set(gca, 'FontSize', FONTSIZE);
end
plot([0,1], [0,1], 'k:')
set(gca, 'FontSize', FONTSIZE);
legend(sprintf('Single condition\n(AUC=%.2f)', auc{1}(1)), ...
       sprintf('Single cond. ground truth\n(AUC=%.2f)', auc{2}(1)), ...
       sprintf('Two conditions\n(AUC=%.2f)', auc{1}(2)), ...
       sprintf('Two cond. ground truth\n(AUC=%.2f)', auc{2}(2)), ...
       'Random', 'Location', 'EastOutside')
set(gca, 'FontSize', FONTSIZE);
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [20 20])
set(gcf, 'PaperPosition', [0 0 10 7])
hold off

modelprobs = {};
modelinds = {};
T = 50:50:1000;
figure(2);
for k=1:2,
  %subplot(1, 2, k);
  [foo, modelinds{k}] = max(results{k}.marlls, [], 2);
  J = sub2ind(size(results{k}.marlls), 1:length(results{k}.genes), modelinds{k}');
  modelprobs{k} = results{k}.marlls(J') - ...
	logsumexp(results{k}.marlls, 2);

  val{k} = all(M == tcomb(modelinds{k}, :), 2);
  %[foo, I1{k}] = sort(modelprobs{k}(1:500), 'descend');
  %mean(val{k}(I1{k}))
  [foo, I2{k}] = sort(modelprobs{k}, 'descend');
  mean(val{k}(I2{k}));
  averates = cumsum(val{k}(I2{k})) ./ (1:1000)';
  plot(T, 100*averates(T), styles2{k});
  set(gca, 'FontSize', FONTSIZE);
  hold on
end
plot([50,1000], (100/8)*[1,1], 'k:')
xlabel('# of most confident genes')
ylabel('Accuracy of MAP predictions (%)')
set(gca, 'FontSize', FONTSIZE);
legend('Single condition', ...
       'Two conditions', ...
       'Random', 'Location', 'EastOutside')
set(gca, 'FontSize', FONTSIZE);
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [20 20])
set(gcf, 'PaperPosition', [0 0 10 7])
axis([50 1000 0 100])
hold off


% figure(1); print -depsc2 figures/toy_link_roc_2010-07-14.eps
% figure(2); print -depsc2 figures/toy_map_accuracy_2010-07-14.eps
