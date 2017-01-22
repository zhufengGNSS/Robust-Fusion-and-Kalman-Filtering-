clear all
close all
clc

addpath util

% See also mybarrier.m




%% parameters
N = 4; % number of agents
n = 2; % agents move on R^2

if n~= 2, error('Wrong dimension'); end

F = 200; % # of frames
sigma_Q = 0.001; % process noise
sigma_R1 = 1; % gps error 
sigma_Rij = 0.1; %  relative position error






vmax = 0.3; % maximum velocity of agents

Mx = 10;
Mv{1} = .02;
Mv{2} = .02;
Mv{3} = .02;
Mv{4} = .02;

E=[1 2;...
   2 3;...
   2 4;...
   3 4;...
   1 4;...
   2 1;...
   3 1;...
   4 1];% directed edges, inward to 1 in order to test robustness to corr

nE = size(E,1);

%% limits of space where agents move
llim = [-10;-10];
ulim = [10;10];

xlim =[-15 15];
ylim =[-15 15];

%% set up  system matrices
II = [eye(n) eye(n)];
IO = [eye(n) zeros(n)];
OI = [zeros(n) eye(n)];
OIt = [zeros(n); eye(n)]; 
OO = [zeros(n) zeros(n)];

% x(t+1) = A*x(t)+B*w(t), w(t) ~ N(0,Q)
for ia=1:N, A{ia} = [II; OI]; end
for ia=1:N, B{ia} = OIt; end
for ia=1:N, Q{ia} = sigma_Q^2*eye(n); end % process noise
for ia=1:N, sQ{ia} = sqrtm(Q{ia}); end


% measurement noise
R1 = sigma_R1^2*eye(n); % GPS
sqR1 = sqrtm(R1);

Rp= sigma_Rij^2*eye(n); % relative position error
sqRij = sqrtm(Rp);

for ie=1:nE, Rij{ie} = Rp; end


% centralized quantities
At = blkdiag(A{1:N}); % centralized KF
Bt = blkdiag(B{1:N}); % centralized KF

Ct = zeros(n*nE+n,2*n*N);
Ct(1:n,1:n) = eye(n); % GPS 


for ie=1:nE
   
    ia = E(ie,1);
    ib = E(ie,2); % edge ia ---> ib
    
    idxr = (n+(ie-1)*n+1):(n+ie*n);
    idxa = (2*n*(ia-1)+1):(2*n*ia);
    idxb = (2*n*(ib-1)+1):(2*n*ib);
    
    Ct(idxr,idxb) = IO;
    Ct(idxr,idxa) = -IO;
end
  
Qt = blkdiag(Q{1:N});
sqQt = sqrtm(Qt);    

Rt = blkdiag(R1,Rij{1:nE});
    


%% initial conditions
for ia=1:N, v0{ia} =Mv{ia}*randn(n,1); end


x0{1} = [-5;-5];
x0{2} = [-5;5];
x0{3} = [5;-5];
x0{4} = [5;5];


for ia=1:N, x{ia} = [x0{ia};v0{ia}]; end



xt = cat(1,x{:});


%% initial estimations

for ia=1:N, xh{ia} = 0*x{ia}; end % my estimator
for ia=1:N, xci{ia} = 0*x{ia}; end % CI estimator
for ia=1:N, xnf{ia} = 0*x{ia}; end % naive fusion estimator

xcf = xt; % centralized KF estimator


for ia=1:N, Sxh{ia} = blkdiag(Mx^2*eye(n),Mv{ia}^2*eye(n)); end
for ia=1:N, Sci{ia} = blkdiag(Mx^2*eye(n),Mv{ia}^2*eye(n)); end
for ia=1:N, Snf{ia} = blkdiag(Mx^2*eye(n),Mv{ia}^2*eye(n)); end

Scf = blkdiag(Sxh{1:N});


%% initial control inputs
for ia=1:N, 
    u{ia}=zeros(n,1);
end


%% error initialization
serrrf = zeros(F,N);
serrci = zeros(F,N);
serrnf = zeros(F,N);
serrcf = zeros(F,N);

perrrf = zeros(F,N); % position error
perrci = zeros(F,N);
perrnf = zeros(F,N);
perrcf = zeros(F,N);

perrrf3s = zeros(F,N); % position error 3s
perrci3s = zeros(F,N);
perrnf3s = zeros(F,N);
perrcf3s = zeros(F,N);


verrrf = zeros(F,N); % velocity error
verrci = zeros(F,N);
verrnf = zeros(F,N);
verrcf = zeros(F,N);

verrrf3s = zeros(F,N); % velocity error 3s
verrci3s = zeros(F,N);
verrnf3s = zeros(F,N);
verrcf3s = zeros(F,N);



%% main loop here

tpause = 0.001;

figure(30000),
plot(x{1}(1),x{1}(2),'bo');
axis([xlim ylim])
hold on;
for ia=2:N,    plot(x{ia}(1),x{ia}(2),'bo');    end
pause(tpause)




for f = 1 : F

 
    %% generate measurements
    
    y{1} = IO*x{1} + sqR1*randn(n,1); % GPS measurement
    
    for ie=1:nE, 
        ia = E(ie,1);
        ib = E(ie,2);
        y{ie+1} = IO*(x{ib} - x{ia}) + sqRij*randn(n,1); 
    end
    
    yt = cat(1,y{:}); % stacked measurements
    
    %% update step
    
    
    % RKF
    
        xhp = xh; % previous
        Sxhp = Sxh;
    
        [xh{1},Sxh{1}] = KalmanUpdate(IO,xh{1},Sxh{1},y{1},R1);
        
        for ie=1:nE,
            ia = E(ie,1);
            ib = E(ie,2);
            [xh{ib},Sxh{ib},~,~] = mybarrier(xh{ib},Sxh{ib},IO,xhp{ia},Sxhp{ia},-IO,y{ie+1},Rij{ie}); % agent ib updates
        end

    
    
    
   % CI
    
        % covariance intersection
    
        xcip = xci; % previous
        Scip = Sci;
        
         [xci{1},Sci{1}] = KalmanUpdate(IO,xci{1},Sci{1},y{1},R1);
    
         for ie=1:nE,
            ia = E(ie,1);
            ib = E(ie,2);
            [xci{ib},Sci{ib},~,~] = CovarianceIntersection(xci{ib},xcip{ia},y{ie+1},Sci{ib},Scip{ia},Rij{ie},IO,-IO);
         end
         
         
    
        %% naive fusion
    
        
        xnfp = xnf; % previous
        Snfp = Snf;
        
        [xnf{1},Snf{1}] = KalmanUpdate(IO,xnf{1},Snf{1},y{1},R1);
        
        for ie=1:nE,
            ia = E(ie,1);
            ib = E(ie,2);
            [xnf{ib},Snf{ib}] =  KalmanUpdate(IO,xnf{ib},Snf{ib},y{ie+1}+IO*xnfp{ia},Rij{ie}+IO*Snfp{ia}*IO');
        end

        
        
        
        
         %% centralized KF
         [xcf,Scf] = KalmanUpdate(Ct,xcf,Scf,yt,Rt);


    %% prediction step
    % same for all methods, similar to KF prediction step

    % control input
    kv = 0.1;
    for ia=1:N, 
        u{ia} = zeros(n,1);
        isOut{ia} = false;
        for idim=1:n
            if (x{ia}(idim) >=ulim(idim)) || (x{ia}(idim) <=llim(idim)) 
                isOut{ia} = true;
            end
            
            if  isOut{ia} 
                u{ia} = -x{ia}(n+1:2*n)  -kv*x{ia}(1:n)/norm(x{ia}(1:n));
                break;
            end
                     
        end     
        
        
        if norm(x{ia}(n+1:2*n))> vmax && ~isOut{ia}
            u{ia} = -kv*x{ia}(n+1:2*n);
        end
        
    end
    ut = cat(1,u{:}); % stacked inputs
    
    
    
    % my method
    for ia=1:N, [xh{ia},Sxh{ia}] = KalmanPredict(A{ia},B{ia},xh{ia},Sxh{ia},Q{ia},u{ia}); end
    
    % covariance intersection
    for ia=1:N, [xci{ia},Sci{ia}] = KalmanPredict(A{ia},B{ia},xci{ia},Sci{ia},Q{ia},u{ia}); end
    
    % naive fusion
    for ia=1:N, [xnf{ia},Snf{ia}] = KalmanPredict(A{ia},B{ia},xnf{ia},Snf{ia},Q{ia},u{ia}); end
      
    % centralized KF
    [xcf,Scf] = KalmanPredict(At,Bt,xcf,Scf,Qt,ut);
    
    
    
     %% measure some errors  
    for ia=1:N
        serrrf(f,ia) = norm(xh{ia}-x{ia});
        serrci(f,ia) = norm(xci{ia}-x{ia});
        serrnf(f,ia) = norm(xnf{ia}-x{ia});
        serrcf(f,ia) = norm(xcf(2*n*(ia-1)+1:2*n*ia)-x{ia});
        
        perrrf(f,ia) = norm(xh{ia}(1:2)-x{ia}(1:2));
        perrci(f,ia) = norm(xci{ia}(1:2)-x{ia}(1:2));
        perrnf(f,ia) = norm(xnf{ia}(1:2)-x{ia}(1:2));
        perrcf(f,ia) = norm(xcf(2*n*(ia-1)+1:2*n*(ia-1)+2)-x{ia}(1:2));
        
        verrrf(f,ia) = norm(xh{ia}(3:4)-x{ia}(3:4));
        verrci(f,ia) = norm(xci{ia}(3:4)-x{ia}(3:4));
        verrnf(f,ia) = norm(xnf{ia}(3:4)-x{ia}(3:4));
        verrcf(f,ia) = norm(xcf(2*n*(ia-1)+3:2*n*(ia-1)+4)-x{ia}(3:4));
        
        
        perrrf3s(f,ia) = 3*sqrt(trace(Sxh{ia}));
        perrci3s(f,ia) =  3*sqrt(trace(Sci{ia}));
        perrnf3s(f,ia) =  3*sqrt(trace(Snf{ia}));
        perrcf3s(f,ia) =  3*sqrt( trace(  Scf(2*n*(ia-1)+1:2*n*(ia-1)+4,2*n*(ia-1)+1:2*n*(ia-1)+4)  ) );
    end
    
    %% simulation step
    
    
     for ia=1:N,
         w{ia} = sQ{ia}*randn(n,1);
        x{ia} = A{ia}*x{ia} + B{ia}*(u{ia}+w{ia});
     end
    
    
    plot(x{1}(1),x{1}(2),'bo');
    axis([xlim ylim])
    hold on;
    for ia=2:N,    
        plot(x{ia}(1),x{ia}(2),'bo');    
    end
    
    for ia=1:N
        idx = 2*n*(ia-1)+1:2*n*(ia-1)+n;
        error_ellipse('C',9*Scf(idx,idx),'mu',xcf(idx),'style','b--'); 
        plot(xcf( 2*n*(ia-1)+1),xcf( 2*n*(ia-1)+2),'b+')
    end
    
    for ia=1:N
        error_ellipse('C',9*Sxh{ia}(1:2,1:2),'mu',xh{ia}(1:2),'style','r--'); 
        plot(xh{ia}(1),xh{ia}(2),'r+')
    end
    
    for ia=1:N
        error_ellipse('C',9*Sci{ia}(1:2,1:2),'mu',xci{ia}(1:2),'style','g'); 
        plot(xci{ia}(1),xci{ia}(2),'g+')
    end
    
    
    for ia=1:N
        error_ellipse('C',9*Snf{ia}(1:2,1:2),'mu',xnf{ia}(1:2),'style','k'); 
        plot(xnf{ia}(1),xnf{ia}(2),'k+')
    end
    
    hold off;
    title(sprintf('Frame %d',f))
    
    pause(tpause)
    
    
    
        
end








%% plots



figure(1),
for ia=1:N
      subplot(2,ceil(N/2),ia);
      plot(perrrf(:,ia),'r');
      hold on;
      plot(perrci(:,ia),'g')
      plot(perrnf(:,ia),'k')
      plot(perrcf(:,ia),'b')
      hold off;
      title(['Position error, agent: ' sprintf('%d',ia)])
      xlabel('k')
      legend('RF','CI','NF','KF') % Robust Fusion, Covariance Interesection, Naive Fusion, Kalman Filter
       
end






figure(2),
for ia=1:N
      subplot(2,ceil(N/2),ia);
      plot(perrrf(:,ia),'r');
      hold on;
      plot(perrci(:,ia),'g')
      plot(perrcf(:,ia),'b')
      hold off;
      title(['Position error, agent: ' sprintf('%d',ia)])
      xlabel('k')
      legend('RF','CI', 'KF') % Robust Fusion, Covariance Interesection,   Kalman Filter
       
end





figure(3),
for ia=1:N
      subplot(2,ceil(N/2),ia);
      plot(verrrf(:,ia),'r');
      hold on;
      plot(verrci(:,ia),'g')
      plot(verrcf(:,ia),'b')
      hold off;
      title(['Velocity error, agent: ' sprintf('%d',ia)])
      xlabel('k')
      legend('RF','CI', 'KF') % Robust Fusion, Covariance Interesection,   Kalman Filter
       
end






