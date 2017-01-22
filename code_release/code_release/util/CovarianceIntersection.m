function [xci,Pxci,f,om] = CovarianceIntersection(x,y,z,Sx,Sy,Sz,C,D)
% function [xci,Pci,f,om] = CovarianceIntersection(x,y,z,Sx,Sy,Sz,C,D)


om = 0.5; % initial point
maxiter = 200;
tol = 10^(-8);
fprev = 10^8;


iPxx = (Sx)^(-1);
iPzz = C'*((D*Sy*D'+Sz)^(-1))*C;


   
for iter = 1:maxiter
    
    X=  (om*iPxx + (1-om)*iPzz );
    iX = (X)^(-1);
    g = - trace (iX*(iPxx-iPzz)*iX);
    
    if iter >1, 
        fprev = f; 
    end
    
    f = trace(iX);
    
    if norm(g) < tol
        break
    end
    
    if fprev-f <tol
        break
    end
   
    
    b = 0.5;
    s = 1; % step size
    
    omn = om - s*g;
    fn = trace( inv(omn*iPxx + (1-omn)*iPzz) );
    
    while (omn>1) || (omn<0) ||  (fn > f)
        s = s*b;
        omn = om - s*g;
        if omn >=0 && omn <=1
            fn = trace( inv(omn*iPxx + (1-omn)*iPzz) );
        end
    end
     
    om = omn;
    
end



Pxci =   (om*iPxx + (1-om)*iPzz )\eye(size(iPxx));
xci = Pxci*(  om*iPxx*x + (1-om)*C'*(D*Sy*D'+Sz)^(-1)*(z-D*y) );  



if norm(Pxci-Pxci')>10^(-9)
    error('Output covariance not symmetric')
end

if ~all(eig(Pxci)>10^(-9))
    error('Output covariance not strictly positive definite')
end


Pxci = .5*(Pxci+Pxci');




end

