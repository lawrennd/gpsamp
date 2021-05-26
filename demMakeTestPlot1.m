
load datasets/drosophila_data.mat;
load datasets/trainScaleDros;

printPlots = 0; 

G = [];

disp('mef2 Top-ranked genes');
for n=1:size(mef2_r,2)
   % 
   %
   I = find(strcmp(mef2_r{n}.genes, drosexp.fbgns));
   prScore = mean(drosexp.fitmean(I, :) ./ sqrt(drosexp.fitvar(I, :)), 2);   
   % select one probe  
   [val j] = max(prScore);
   I = I(j); 
   
   
   Gene = drosexp.fitmean(I, :)/sc;
   Gene = reshape(Gene,1,12,3);
   GeneVar = drosexp.fitvar(I, :)/(sc.^2); 
   GeneVar = reshape(GeneVar,1,12,3);
   
   % make the plot 
   gpmtfTestPlot(modelTest, Gene, GeneVar,  mef2_r{n}.genes, mef2_r{n}.testGene, TFs, 'dros', printPlots);
   
   fprintf(1,'Gene #%2d Average log likelihood %5f\n',n, mean(mef2_r{n}.testGene.LogL));     
   disp('Press any key to continue');
   
   %pause;
   
   G = [G, I]; 
   close all; 
   
   %
   %
end


for n=1:size(twi_r, 2)
   %
   %
   
   I = find(strcmp(twi_r {n}.genes, drosexp.fbgns));
   
   prScore = mean(drosexp.fitmean(I, :) ./ sqrt(drosexp.fitvar(I, :)), 2);   
   % select one probe  
   [val j] = max(prScore);
   I = I(j); 
   
   Gene = drosexp.fitmean(I, :)/sc;
   Gene = reshape(Gene,1,12,3);
   GeneVar = drosexp.fitvar(I, :)/(sc.^2); 
   GeneVar = reshape(GeneVar,1,12,3);
      
   % make the plot 
   gpmtfTestPlot(modelTest, Gene, GeneVar,  twi_r {n}.genes, twi_r{n}.testGene, TFs, 'dros', printPlots);
   
   
   fprintf(1,'Gene #%2d Average log likelihood %5f\n',n, mean(twi_r{n}.testGene.LogL));     
   disp('Press any key to continue');
   
   %pause;
   
   G = [G, I]; 
    
   close all; 
   %
   %
end


