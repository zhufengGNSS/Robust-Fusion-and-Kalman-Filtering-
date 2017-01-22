function X = myivech(x)

% x is n(n+1)/2 x 1

m = numel(x);
n = round((-1+sqrt(1+8*m))/2);

X = zeros(n);

X(find(tril(ones(n)))) = x;

X = X + X';
X(1:n+1:n^2) = X(1:n+1:n^2)/2;


end