function [xh,Sxh,K] = estimate_MLE(x,Sx,y,Sy,Q)


% SIMPLE FUSION:  z = x + y

n = size(x,1); % = size(y,1)

A = Sx  + Sy - (Q+Q');
B = Sx - Q';
X = A\B;
K=X';


xh = x + K*(y-x);

Sxh = [eye(n)-K K]*[Sx Q;Q' Sy]*[eye(n)-K' ;  K'];

end

