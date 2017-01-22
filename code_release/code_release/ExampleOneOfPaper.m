clear all
close all
clc

% fusion of two normal random variables  with unknown correlation

addpath util


n = 2; % dimension of x, y

x = [0;0];
Sx = [5 0;0 5];


y = [0;0]; 
Sy = [3 0;0 7];


%% mymethod
R = zeros(n);
C = eye(n);
D = eye(n);
z = C*x+D*y;


[~,Sxh,K,Q] =mybarrier(x,Sx,C,y,Sy,D,z,zeros(n));
xh = (eye(2)-K)*x+K*y;



%% MLE ESTIMATE
% generate some random correlation
tt =0.8;
Sxy = tt*[3 0;0 5];
SS = [Sx Sxy;Sxy' Sy];

if ~all( eigs(SS) >=0)
    error('Not a valid correlation')
end


A = Sx  + Sy - (Sxy+Sxy');
B = Sx - Sxy';

[xmle,Sxmle,Kmle] = estimate_MLE(x,Sx,y,Sy,Sxy); 




%% CI estimate

[xci,Pxci,f,omopt] = CovarianceIntersection(x,y,z,Sx,Sy,R,C,D);






%% PLOT RESULTS


%% PLOTS RESULTS

figure,
error_ellipse('C',Pxci ,'mu',xci,'style','g'); 
hold on;
error_ellipse('C', Sxh,'mu',xh,'style','r')
error_ellipse('C',Sx ,'mu',x,'style','b'); 
error_ellipse('C',Sxmle ,'mu',xmle,'style','k'); 

plot(xci(1),xci(2),'g+')
plot(xh(1),xh(2),'r+')
plot(x(1),x(2),'b+')
plot(xmle(1),xmle(2),'k+')

hold off;
legend('CI','RF','Initial','MLE')
axis equal

%% DISPLAY RESULTS
display('============ Covariance intersection mean and covariance ==========')
xci
Pxci 

display('============ Robust fusion mean and covariance ==========')
xh
Sxh

display('============ MLE mean and covariance ==========')
xmle
Sxmle

