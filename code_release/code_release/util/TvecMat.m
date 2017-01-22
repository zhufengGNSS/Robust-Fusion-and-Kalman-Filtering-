function Tmn = TvecMat(m,n)

% Tmn*vec(A) = vec(A')


d = m*n;


i = 1:d;
rI = 1+m.*(i-1)-(m*n-1).*floor((i-1)./n);

Tmn = sparse(rI,i,ones(1,d),d,d);