function x = myvech(X)

% X must be symmetric

x = X(find(tril(ones(size(X,1)))));



end

