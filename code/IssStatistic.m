function [W,info] = IssStatistic(X, X_ko, y,tol)

n = size(X,2);
opts.history = true;
[u,p,info] = fast_iss([X,X_ko], y,tol,opts);
Z = zeros(1,2*n);
for i = length(info.hist_t):-1:1
    Z(info.hist_in(i)) = 1/info.hist_t(i);
end
W = max([Z(1:n) ; Z((n+1):(2*n))]).*sign(Z(1:n) - Z((n+1):(2*n)));
