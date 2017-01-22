clear all
close all
clc

% fusion with partial measurements


addpath util


n = 2; % dimension of x, y

x = [0;0];
Sx = [5 0;0 5];


y = 0.0; % zero mean measurement noise, possibly correlated to state x
Sy = 1;


%% Robust fusion

R = 0;
C = [1 0];
D = 1;
z = C*x+D*y;


[~,Sxh,K,Q] = mybarrier(x,Sx,C,y,Sy,D,z,R);
xh = x+K*(z-C*x);



%% CI estimate
[xci,Pxci,f,omopt] = CovarianceIntersection(x,y,z,Sx,Sy,R,C,D);



%% PLOTS RESULTS

figure,
error_ellipse('C',Pxci ,'mu',xci,'style','g'); 
hold on;
error_ellipse('C', Sxh,'mu',xh,'style','r')
error_ellipse('C',Sx ,'mu',x,'style','b'); 

plot(xci(1),xci(2),'g+')
plot(xh(1),xh(2),'r+')
plot(x(1),x(2),'b+')

hold off;
legend('CI','RF','Initial')
axis equal

%% 
display('============ Covariance intersection mean and covariance ==========')
xci
Pxci 

display('============ Robust fusion mean and covariance ==========')
xh
Sxh
