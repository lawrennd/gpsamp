function fx = jointactFunc(LikParams,ff,J)
%
%
%
%

SizF = size(ff,2);

%NumGenes = size(J(:),1);
W = LikParams.W(J,:);
W0 = LikParams.W0(J); 

%
switch LikParams.jointAct
    case 'genHill'
       % generalized hill function function 
       %ff(ff==0)=eps;
       %fx = W*log(ff) + W0(:, ones(1, SizF)); %repmat(W0,[1 SizF]);
       fx = W*log(ff+1e-100) + W0(:, ones(1, SizF)); %repmat(W0,[1 SizF]);
       %W0 = sum(W,2).*log(Par.Gammas(J)); 
       %
       fx = 1 ./ (1 + exp(-fx));
       %fx = sigmoid(fx); 
    case 'lin'
       % 
       fx = W*ff + W0(:, ones(1, SizF)); %repmat(W0,[1 SizF]); 
       %
    case 'sigmoid'
       %
       fx = W*ff + W0(:, ones(1, SizF)); %repmat(W0,[1 SizF]);
       % pass it through the sigmoid function 
       fx = sigmoid(fx);
    case 'michMenten'
       % 
       fx = michMenten(ff, W, LikParams.Net_X(J,:));
       %
    %case 'michMentenAct'
    %   xp = repmat(ff,[NumGenes 1]); 
    %   fx = xp./(xp + repmat(exp(-W0),[1 SizF]));
    %case 'michMentenRepres'
    %   xp = repmat(ff,[NumGenes 1]); 
    %   fx = 1./(xp + repmat(exp(-W0),[1 SizF])); 
end

% binarize the outputs if necessary 
if LikParams.jointActBin == 1
    fx = round(fx);
end

