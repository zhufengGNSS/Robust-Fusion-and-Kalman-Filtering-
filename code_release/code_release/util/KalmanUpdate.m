function [xhat,Sxhat] = KalmanUpdate(C,xhat,Sxhat,y,R)



CSt = (C*Sxhat)';
K = CSt*(C*CSt+R)^(-1);
xhat = xhat + K*(y-C*xhat);
Sxhat = Sxhat - K*CSt';



end

