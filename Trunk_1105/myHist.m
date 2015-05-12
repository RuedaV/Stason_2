function [] = myHist(n, p, q, m, X_t, X, MSPE, MAPE,  x0, nbins, bool)
% Here I visualize my results for simulations.
% Namely, for each parameter I built a histogram
% we consider the following model
% r(t) = mu_ + theta*d(t) + rho*r(t-1)  + delta*h2(t) + e(t)
% e(t) = h(t)*z(t)
% h2(t) = omega + alpha(1)*h2(t-1) + ... + alpha(p)*h2(t-p) + 
%               +  beta(1)*e2(t-1) + ... +  beta(q)*e2(t-q) +
%       + gamma(1)*e2(t-1)*I{e(t-1)>0} + ... + gamma(m)*e2(t-m)*I{e(t-m)>0}
% MyHistEl( y, target, nbins)

MyHistEl(X_t(:,4+n), X(:,4+n), x0(4+n), nbins, '$\frac{\hat\omega-\omega}{\omega}$');

if p == 1
    MyHistEl(X_t(:,5+n), X(:,5+n), x0(5+n), nbins, '$\frac{\hat\beta-\beta}{\beta}$');
else
    j = 1;
    for i = 5+n:4+n+p
        MyHistEl(X_t(:,i), X(:,i), x0(i), nbins, strcat ('Histogram for $\beta_$', num2str(j)));
        j = j + 1;
    end    
end

if q == 1
    MyHistEl(X_t(:,5+n+p), X(:,5+n+p), x0(5+n+p), nbins, '$\frac{\hat\alpha-\alpha}{\alpha}$');    
else
    j = 1;
    for i = 5+n+p: 4+n+p+q
        MyHistEl(X_t(:,i), X(:,i), x0(i), nbins, strcat ('Histogram for $\alpha_$', num2str(j)));        
        j = j + 1;
    end 
end

if m == 1
    MyHistEl(X_t(:,5+n+p+q), X(:,5+n+p+q), x0(5+n+p+q), nbins, '$\frac{\hat\gamma-\gamma}{\gamma}$');    
else  
    j = 1;
    for i = 5+n+p+q:4+n+p+q+m
        MyHistEl(X_t(:,i), X(:,i), x0(i), nbins, strcat ('Histogram for $\gamma_$', num2str(j)));
        j = j + 1;
    end  
end

% MyHistEl(MSPE, 0, nbins, '$\frac{MSPE_f-MSPE_t}{MSPE_t}$');
% 
% MyHistEl(MAPE, 0, nbins, '$\frac{MAPE_f-MAPE_t}{MAPE_t}$');




end

