
dataName = 'regressionToy';
expNo = 1;
storeRes = 0;

% load outputs and inputs
% create some data or load data from somewhere
%[Y, X] = loadDataset(dataName)
X = (0:0.01:1)';
sigma2 = 0.3^2; 
Y = sin(2*pi*X) + cos((5*pi*X)) + sqrt(sigma2)*randn(101,1);
plot(X, Y, '+k', 'lineWidth', 3); 
[n D] = size(X); 
T = size(Y,2);


% model options
options = gpsampOptions('regression'); 

% create the model
model = gpsampCreate(Y, X, options);

mcmcoptions = mcmcOptions('controlPnts');
mcmcoptions.adapt.incrNumContrBy = 1;
% adaption phase
[model PropDist samples accRates] = gpsampControlAdapt(model, mcmcoptions.adapt);
% training/sampling phase
[model PropDist samples accRates] = gpsampControlTrain(model, PropDist, mcmcoptions.train);

if storeRes == 1
    d = date; 
    save(['dem' dataName num2str(experimentNo) d '.mat'], 'model', 'PropDist','samples','accRates');
end

% plot the samples
gpsampPlot(model, samples);
