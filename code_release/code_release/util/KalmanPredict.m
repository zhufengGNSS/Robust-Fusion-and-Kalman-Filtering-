function [xhat,Sxhat] = KalmanPredict(A,B,xhat,Sxhat,Q,u)

% x(t+1) = A x(t) + B (u(t)+w(t))

xhat = A*xhat+ B*u;
Sxhat = A*Sxhat*A' + B*Q*B';



end