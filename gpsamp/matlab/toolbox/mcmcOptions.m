function mcmcoptions = mcmcOptions(algorithm) 
%


switch algorithm 
    case 'controlPnts'
        
        % MCMC OPTIONS (you can change these options if you like)   
        mcmcoptions.train.StoreEvery = 100; % store samples after burn in  every StoreEvery iterations
        mcmcoptions.train.Burnin = 5000;  % burn in time
        mcmcoptions.train.T = 50000; % sampling time
        mcmcoptions.train.Store = 0;  % store the results regularly in a file 
        
        mcmcoptions.train.disp = 0; 
        
        % options for the adaptive phase in MCMC (you could change this, although not recommended) 
        mcmcoptions.adapt.T = 200;          
        mcmcoptions.adapt.Burnin = 100;
        mcmcoptions.adapt.StoreEvery = 10; 
        mcmcoptions.adapt.disp = 1;
        mcmcoptions.adapt.initialNumContrPnts = 3; 
        mcmcoptions.adapt.incrNumContrBy = 1;
end