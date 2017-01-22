function Pis = MatVech2Vec(n)


ns = n*(n+1)/2;

% vec(S) = Pis * vech(S), where S is a symmetric matrix


idxl = find(tril(ones(n)))';
 


idxs = find(~tril(ones(n)))';

cI = zeros(size(idxs));

for i=1:numel(idxs)
    
    ri = mod(idxs(i)-1,n)+1;
    ci = floor( (idxs(i)-1)/n)+1;
    k= n-ri+1;
    l= n-ci;
    cI(i) = ns - (k*(k-1)/2 + l)   ;% position in vech
    
end


%Pis(sub2ind([n^2 ns],idxs,cI)) = 1; % stricly upper triangular part

Pis = sparse([idxl idxs],[1:ns cI],1,n^2,ns);