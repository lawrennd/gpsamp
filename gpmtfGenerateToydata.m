function [LikParams lengthscale] = gpmtfGenerateToydata(NumGs, TimesG, NumTFs, TimesF)
%
%

% generate NumTFs different GP functions using the exponential kernel
n = size(TimesF(:),1);
SizG = size(TimesG(:),1)
sigma2 = 0.01; 
sigmaf = 1;
lengthscale = [4 4 3.5];
contr = [n+1 n+2 n+3 n+4 n+5 n+6 n+7 n+8]; 
sc = max(TimesG)/11;
FF1 = [-5 -2 0.5 -1 -2 -3 -3 -4;
      -5 -3 0 -0.5 -1 -2 -3 -3;
      -5 -4 -4 -1 0 -4 -4 -5];
  
FF2 = [-5 -2 0.5 -1 -2 -3 -3 -4;
      -5 -4 -4 -3 -3 -1 0 -3;
      -5 -4 -4 -1 0 -4 -4 -5];

FF3 = [-5 -1 -0.5 0 -0.5 -2 -3 -4;
      -5 -4 -4 -3 -3 -1 0 -3;
      -5 -4 -4 -1 0 -4 -4 -5];
% put some noise 
FF1 = FF1 + 0.05*randn(size(FF1));
FF2 = FF2 + 0.05*randn(size(FF2));
  
step = (TimesF(end) - TimesF(1))/(size(contr,2)-1);  
Xu = TimesF(1):step:TimesF(end);


FF = FF1(1:NumTFs,:);
un = size(Xu,2);
m = zeros(n+un,1);
for i =1:NumTFs
    covfunParams.logtheta = [0.5*log(lengthscale(i)) 0.5*log(sigmaf) 0.5*log(sigma2)];
    K = kernCompute(covfunParams, [TimesF(:); Xu(:)], [TimesF(:); Xu(:)]);
   
    [cmuMinus, cSigma, KInvKu] = gaussianFastConditional(m', K, 1:n,contr);
    [L,er]=jitterChol(cSigma);
    if er>0, L = real(sqrtm(cSigma)); end
    cmu = cmuMinus + FF(i,:)*KInvKu;
    F(i,:) = gaussianFastSample(1, cmu, L);
    if i ==1
    plot(TimesF,log(1+exp(F(i,:))));
    hold on;
    elseif i==2
    plot(TimesF,log(1 + exp(F(i,:))),'r');
    elseif i==3
    plot(TimesF,log(1 + exp(F(i,:))),'g');
    else
    plot(TimesF,log(1 + exp(F(i,:))),'c');    
    end
end

F1 = F;
pause

FF = FF2(1:NumTFs,:);
un = size(Xu,2);
m = zeros(n+un,1);
for i =1:NumTFs
    covfunParams.logtheta = [0.5*log(lengthscale(i)) 0.5*log(sigmaf) 0.5*log(sigma2)];
    K = kernCompute(covfunParams, [TimesF(:); Xu(:)], [TimesF(:); Xu(:)]);
   
    [cmuMinus, cSigma, KInvKu] = gaussianFastConditional(m', K, 1:n,contr);
    [L,er]=jitterChol(cSigma);
    if er>0, L = real(sqrtm(cSigma)); end
    cmu = cmuMinus + FF(i,:)*KInvKu;
    F(i,:) = gaussianFastSample(1, cmu, L);
    if i ==1
    plot(TimesF,log(1+exp(F(i,:))));
    hold on;
    elseif i==2
    plot(TimesF,log(1 + exp(F(i,:))),'r');
    elseif i==3
    plot(TimesF,log(1 + exp(F(i,:))),'g');
    else
    plot(TimesF,log(1 + exp(F(i,:))),'c');    
    end
end

F2 = F;


FF = FF3(1:NumTFs,:);
un = size(Xu,2);
m = zeros(n+un,1);
for i =1:1%NumTFs
    covfunParams.logtheta = [0.5*log(lengthscale(i)) 0.5*log(sigmaf) 0.5*log(sigma2)];
    K = kernCompute(covfunParams, [TimesF(:); Xu(:)], [TimesF(:); Xu(:)]);
   
    [cmuMinus, cSigma, KInvKu] = gaussianFastConditional(m', K, 1:n,contr);
    [L,er]=jitterChol(cSigma);
    if er>0, L = real(sqrtm(cSigma)); end
    cmu = cmuMinus + FF(i,:)*KInvKu;
    F(i,:) = gaussianFastSample(1, cmu, L);
    if i ==1
    plot(TimesF,log(1+exp(F(i,:))));
    hold on;
    elseif i==2
    plot(TimesF,log(1 + exp(F(i,:))),'r');
    elseif i==3
    plot(TimesF,log(1 + exp(F(i,:))),'g');
    else
    plot(TimesF,log(1 + exp(F(i,:))),'c');    
    end
end


F3 = F(1,:);

clear F;
F(:,:,1) = F1;
F(:,:,2) = F2;
save FF2 F F3;
pause
return; 






% produce Random interaction Weights (the last column are the bias terms) 
W = randn(NumGs,NumTFs);
ok = zeros(NumGs,NumTFs);
for j=1:NumGs  
  ch = randperm(NumTFs);
  ok(j,ch(1)) = 1;
  if (NumTFs > 1) & (rand > 0.75)
      ok(j,ch(2)) = 1;
  end
end
W = W.*ok;

Net_X = (-1)*(W<0) + 1*(W>0); 
Net_X = ok;    
  
W0 = randn(NumGs,1);

% produce random kinetic parameters for each gene 
for j=1:NumGs
    %B D S A 
    %ok = randperm(4);% randn;
    %Kinetics(j,:) = rand(1,4).*ok; 
    Kinetics(j,:) = 2*rand(1,4)-1;
    Kinetics(j,:) = exp(Kinetics(j,:));%exp(rand(1,4));
    % decays fixed to one
    %Kinetics(j,2) = 1; 
    Kinetics(j,4) = Kinetics(j,1)/Kinetics(j,2);  
end


% produce random kinetic parameters for each TF-gene 
for j=1:NumTFs
    %D S A 
    KineticsTF(j,:) = 2*rand(1,2)-1;
    KineticsTF(j,:) = exp(KineticsTF(j,:));
    %
end
KineticsTF(:,2) = ones(NumTFs,1);

% solve numerically the differential equation
LikParams.numTFs = NumTFs;
LikParams.F = F;
LikParams.Genes = zeros(NumGs,SizG);
LikParams.GenesTF = zeros(NumTFs,SizG);
LikParams.TimesG = TimesG; 
LikParams.TimesF = TimesF;
LikParams.kinetics = Kinetics;
LikParams.kineticsTF = KineticsTF;
LikParams.Net_X = Net_X;
LikParams.W = W;
LikParams.W0 = W0;
LikParams.sigmas = sigma2*ones(NumGs,SizG);
LikParams.sigmasTF = sigma2*ones(NumTFs,SizG);
LikParams.singleAct = 'logOnePlusExp';
LikParams.jointAct = 'genHill';
LikParams.jointActBin = 0;
LikParams.startTime = 1;

step = TimesF(2)-TimesF(1);
NofTaus = floor(abs(TimesF(1))/step);
NofTaus
ok = 0:step:TimesG(end);
if ok(end) < TimesG(end)
    ok = [ok, TimesG(end)];
end

LikParams.sizTime = size(ok,2);
sizTime = LikParams.sizTime; 
LikParams.startTime = size(TimesF,2) - sizTime + 1;
LikParams.startTime
TimesF(LikParams.startTime)
LikParams.Tausreadme = 'taus are the delays in the ODEs'; 

for j=1:NumGs
  LikParams.Taus(j) = 0; 
 % ok = randperm(NofTaus)-1;
  ok = 0;%ok(1); 
 % size(TimesF,2) - sizTime + 1 - ok;
   LikParams.Tausindex(j) = size(TimesF,2) - sizTime + 1 - ok;
 % LikParams.Taus(j) = ok*step; 
end 
LikParams.Taus
LikParams.Tausindex

uu = LikParams.TimesF(LikParams.startTime:end);
[commonSlots, comInds] = intersect(uu,LikParams.TimesG);
LikParams.comInds = comInds; 

LikParams.TF = [];
[loglikval, Genes] = gpmtfLogLikelihoodgene(LikParams, F, 1, 1:NumGs);
[commonSlots, comInds] = intersect(TimesF,TimesG);
uu = TimesF(LikParams.startTime:end);
[commonSlots, comInds] = intersect(uu,TimesG);
Genes = Genes(:,comInds);
% add random noise 
Genes = Genes + sqrt(sigma2)*randn(size(Genes));

uu = LikParams.TimesF;
[commonSlots, comIndsTF] = intersect(uu,LikParams.TimesG);
LikParams.comIndsTF = comIndsTF;
[loglikval, GenesTF] = gpmtfLogLikelihoodgeneTF(LikParams, F, 1, 1:NumTFs);

% add random noise 
GenesTF = GenesTF + sqrt(sigma2)*randn(size(GenesTF));
uu = TimesF;
[commonSlots, comInds] = intersect(uu,TimesG);
GenesTF = GenesTF(:,comInds);

LikParams.TF = gpmtfComputeTF(LikParams, F,  1:NumTFs);
LikParams.Genes = Genes; 
LikParams.GenesTF = GenesTF;