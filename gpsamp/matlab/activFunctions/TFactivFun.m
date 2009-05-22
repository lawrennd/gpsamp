function fx = TFactivFun(LikParams,ff,J)
%
%
%
%


SizF = size(ff,2);
NumGenes = size(J(:),1);

% perform the individual transforms first  
ff = feval(LikParams.TFsingleAct,ff);
W = LikParams.W(J,:); 
W0 = LikParams.W0(J); 
if strcmp(LikParams.TFjointAct,'michMenten')
Net_X = LikParams.Net_X(J,:);
end
switch LikParams.TFjointAct
    case 'lin'
       fx = W*ff + repmat(W0,[1 SizF]); 
    case 'sigmoid'
       xp = W*ff + repmat(W0,[1 SizF]);
       % pass it through the sigmoid function 
       fx = sigmoid(xp);
    case 'michMenten'
       fx = michMenten(ff, W, Net_X);
    %case 'michMentenAct'
    %   xp = repmat(ff,[NumGenes 1]); 
    %   fx = xp./(xp + repmat(exp(-W0),[1 SizF]));
    %case 'michMentenRepres'
    %   xp = repmat(ff,[NumGenes 1]); 
    %   fx = 1./(xp + repmat(exp(-W0),[1 SizF])); 
end

% binarize the outputs if necessary 
if LikParams.TFjointActBin == 1
    fx = round(fx);
end

