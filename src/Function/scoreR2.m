%% scoreR2 Coefficient of Determination
% 
%  From definition of Score: v/u 
% ``v = ((y_true - y_pred)** 2).sum()``
% ``u =((y_true - y_true.mean()) ** 2).sum()``

function score = scoreR2(Model, X, Ytrue)
    [ Ypredict, pp] = predict(Model, X);
    Yaverage = mean(Ytrue);
    u = sum(power(( Ytrue - Yaverage), 2));
    v = sum(power(( Ytrue - Ypredict), 2));
    score = v/u;
end