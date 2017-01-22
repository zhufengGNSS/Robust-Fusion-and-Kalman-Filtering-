function [xh,Sxh,K,Q] =mybarrier(x,Si,Cii,y,Sj,Cij,z,Rij)

% measurement model: z = Cii*x + Cij*D + measurement noise
% Si: error covariance of x
% Sj: error covariance of x
% Rij: measurement noise covariance

 

% barrier method parameters, for outer loop
maxiter =  10; % not really necessary parameter
tol = 10^(-3);
mu = 20;
t_outer = 1;

% Infeasible start newton method parameters, inner loop
a = 0.01;
b = 0.5;
tol_inner = 10^(-3);
maxiter_inner = 200;





nx = size(Si,1);
if size(Si,1)~= size(Si,2)
    error('Input covariance not square')
end

ny = size(Sj,1);
if size(Sj,1)~= size(Sj,2)
    error('Input covariance not square')
end

m = size(Cii,1);

if size(Cii,1)~= size(Cij,1)
   error('Dimension mismatch') 
end

if size(Cii,2)~= nx
   error('Dimension mismatch') 
end

if size(Cij,2)~= ny
   error('Dimension mismatch') 
end







iSx = Si\eye(size(Si));
iSy = Sj\eye(size(Sj));

iSyr = sqrtm(iSy);


E = Cii*Si*Cii' + Cij*Sj*Cij' + Rij;


% useful for vectorization of equations
Pix  = TvecMat(nx,m) ; % Pi*X(:) = X'(:), sparse   
Pinxy = TvecMat(ny,nx);

Inx = eye(nx);
Iny = eye(ny);





% initial values, strictly feasible
X = E \ (Cii*Si); % X = K^T, Kalman gain K is nx x m
Q = zeros(nx,ny); % Sigma_xy ,correlation, nx x ny



for k = 1 : maxiter % outer iteration, barrier method
    
 
    
    % residual progress
    cost = zeros(maxiter_inner,1);


   for l=1:maxiter_inner
    

  
  
    
    iSx_Q = iSx*Q;  
    Qt_iSx_Q = Q'*iSx_Q; 
    Cii_Q_Cijt = Cii*Q*Cij';
    E_Cii_Q_Cijt_Cii_Q_Cijtt =   (E + Cii_Q_Cijt +  Cii_Q_Cijt');
    
     Cijt_X = Cij'*X;
     Ciit_X = Cii'*X;
     Xt_Cij = Cijt_X'; 
     
     double_t_outer = 2*t_outer;
     
    f1Q = iSyr*Qt_iSx_Q*iSyr - Iny;
    if1Q =  f1Q \ Iny;
    
      iSyr_if1Q_iSyr = iSyr*if1Q*iSyr;
       temp1 = iSyr*if1Q*iSyr;
   
   iSx_Q_temp1 =  iSx_Q*temp1;
    Cii_Si = Cii*Si;
    
    
    gradfX = double_t_outer*(  E_Cii_Q_Cijt_Cii_Q_Cijtt*X - (Cii_Si+Cij*Q')    );
    gradfQ = double_t_outer*(   Ciit_X*Xt_Cij  -Xt_Cij ) + 2* (iSx_Q*iSyr_if1Q_iSyr);



    
    
    double_t_outer_E_Cii_Q_Cijt_Cii_Q_Cijtt  = double_t_outer*(E_Cii_Q_Cijt_Cii_Q_Cijtt);
    
    A11 =     kron( Inx , double_t_outer_E_Cii_Q_Cijt_Cii_Q_Cijtt);

    

  Cijt_X_kron_Cii = kron(   Cijt_X', Cii   );
  Ciit_Xt_kron_Cij  = kron(  Ciit_X',Cij);
  Cij_t_kron_Ciit_X =  kron( Cij' , Ciit_X );
   
 
       Inx_kron_Cij  = kron( Inx  , Cij  );
      Cijt_kron_Inx = kron( Cij',Inx );

   A22 = kron(temp1', iSx) - kron( iSx_Q_temp1', iSx_Q*temp1 )*Pinxy - kron(temp1', iSx_Q_temp1*iSx_Q');
  
    A12 = Cijt_X_kron_Cii+ (Ciit_Xt_kron_Cij - Inx_kron_Cij)*Pinxy  ;
    A12 = double_t_outer*A12;
    
   A21 = Cijt_X_kron_Cii' ;

   A21 = A21 +   (Cij_t_kron_Ciit_X- Cijt_kron_Inx)*Pix ;
   
   A21 = double_t_outer*A21;    
   A22 = 2*A22;
       
     b1=  -myvec(gradfX);
     b2 = -myvec(gradfQ);
         


           AA = [A11 A12; ... 
                A21 A22];
           bb = [b1;b2];   
 
       DXQ = AA\bb;
        
         dX = reshape(DXQ(1:m*nx),[m nx]); % Newton step
     dQ = reshape(DXQ(m*nx+1:end),[nx ny]); % Newton step
       

            % backtracking, residual of optimality  
            t = 1;

            X_t =  X + t*dX;
            Q_t =  Q + t*dQ;
            nres = norm([myvec(gradfX);myvec(gradfQ)]);
            
            
   
           while min(real(eig(-iSyr*(Q_t'*iSx*Q_t)*iSyr + eye(ny))))<=0  % strict feasibility check
             
                    t = b*t;
                    X_t =  X + t*dX;
                    Q_t =  Q + t*dQ;     
           end

 
           f1Qt = iSyr*(Q_t'*iSx*Q_t)*iSyr - Iny;
            if1Qt = f1Qt \ Iny;

             Cii_Q_t_Cij_t = Cii*Q_t*Cij';
             gradfXt0 = (E + Cii_Q_t_Cij_t + Cii_Q_t_Cij_t')*X_t - (Cii_Si+Cij*Q_t') ;
            
              gradfQt_a = Cii'*(X_t*X_t')*Cij -X_t'*Cij;
             gradfQt_b = (iSx*Q_t*iSyr*if1Qt*iSyr);

              gradfXt = double_t_outer *gradfXt0;
             gradfQt = double_t_outer *gradfQt_a + 2*gradfQt_b;
  

            nres_t = norm([myvec(gradfXt);myvec(gradfQt)]);
         
            while  (nres_t > (1-a*t) *nres   )
         
                    t = b*t;
                    X_t =  X + t*dX;
                    Q_t =  Q + t*dQ;
                    
                    f1Qt = iSyr*(Q_t'*iSx*Q_t)*iSyr - Iny;

                    if1Qt = f1Qt \ Iny;
                    
                    gradfXt = double_t_outer *(  (E + (Cii*Q_t*Cij') +  (Cii*Q_t*Cij')')*X_t - (Cii*Si+Cij*Q_t')    );
                    gradfQt = double_t_outer *(   Cii'*(X_t*X_t')*Cij -X_t'*Cij  ) + 2*iSx*Q_t*iSyr*if1Qt*iSyr;
                    nres_t = norm([myvec(gradfXt);myvec(gradfQt)]);
                    
            end

  
        X = X_t;
        Q = Q_t;
        nres = nres_t;

        

        
        cost(l) =  nres;

        
        if nres < tol_inner % m/t < eps
            break;
        end
        
     
        
        
   end % inner loop ends
        


    if 1/t_outer < tol % m/t < eps, we have only one inequality
        break;
    else
        t_outer = mu*t_outer;
    end
    
    
end



Sxh = [eye(nx)-X'*Cii  -X'*Cij] * [Si Q;Q' Sj] * [eye(nx)-X'*Cii  -X'*Cij]' + X'*Rij*X;

xh = x + X'*(z-Cii*x-Cij*y);

K = X'; % kalman gain