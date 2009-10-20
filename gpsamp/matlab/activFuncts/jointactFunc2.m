function fx = jointactFunc2(Par,ff,J)
%
%
%
%

NumGenes = size(J(:),1);

W = Par.W(J,:);
Tausindex = Par.Tausindex(J); 

W0 = Par.W0(J); 
if strcmp(Par.jointAct,'michMenten')
Net_X = Par.Net_X(J,:);
end
switch Par.jointAct
    case 'lin'
       % 
       % this account for delays in gene expression
       for m=1:size(J,2)
         fx(m,:) = W(m,:)*ff(:,Tausindex(m):Tausindex(m)+Par.sizTime-1);
       end
       fx = fx + repmat(W0,[1 Par.sizTime]);
       %
    case 'sigmoid'
       %
       % this account for delayes in gene expression
       for m=1:size(J,2)
         fx(m,:) = W(m,:)*ff(:,Tausindex(m):Tausindex(m)+Par.sizTime-1);
       end
       fx = fx + repmat(W0,[1 Par.sizTime]); 
       % pass it through the sigmoid function 
       fx = sigmoid(fx);
    case 'michMenten'
       ops.Net_X = Net_X; 
       ops.sizTime = Par.sizTime; 
       ops.Tausindex = Tausindex;
       fx = michMenten(ff, W, ops);
    case 'genHill'
       % generalized hill function function  
       for m=1:size(J,2)
         fx(m,:) = W(m,:)*log(ff(:,Tausindex(m):Tausindex(m)+Par.sizTime-1));
       end
       %W0 = sum(W,2).*log(Par.Gammas(J)); 
       %
       fx = fx + repmat(W0,[1 Par.sizTime]); 
       fx = sigmoid(fx); 
    %case 'michMentenAct'
    %   xp = repmat(ff,[NumGenes 1]); 
    %   fx = xp./(xp + repmat(exp(-W0),[1 SizF]));
    %case 'michMentenRepres'
    %   xp = repmat(ff,[NumGenes 1]); 
    %   fx = 1./(xp + repmat(exp(-W0),[1 SizF])); 
end

% binarize the outputs if necessary 
if Par.jointActBin == 1
    fx = round(fx);
end

