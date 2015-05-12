function [return_value] = s_logL(model, misspType, r, d, x)
% misspType belongs to {0, 1, 2}
% misspType = 0 for mu  + rho
% misspType = 1 for mu + delta
% misspType = 2 for mu + theta

    
    if (misspType == 0)      
        mu_ = x(1,1);
        rho = x(1,2);
        omega = x(1,3);
        alpha = x(1,4);
        beta  = x(1,5);
        gamma = x(1,6);
        h2 = s_condVariance(model, misspType, r, d, mu_, rho, omega, alpha, beta, gamma);        
        
        return_value = sum( -log(h2(2:end,1)) - (r(2:end,1) - mu_ -  rho*r(1:end-1,1)).^2./h2(2:end,1));    
     
    elseif (misspType == 1)      
        mu_   = x(1,1);
        delta = x(1,2);
        omega = x(1,3);
        alpha = x(1,4);
        beta  = x(1,5);
        gamma = x(1,6);
        h2 = s_condVariance(model, misspType, r, d, mu_, delta, omega, alpha, beta, gamma);        
        
        return_value = sum( -log(h2(1:end,1)) - (r(1:end,1) - mu_ - delta*h2(1:end,1)).^2./h2(1:end,1));  
        
    elseif (misspType == 2)
                
        theta = x(1,1:5);
        omega = x(1,6);
        alpha = x(1,7);
        beta  = x(1,8);
        gamma = x(1,9);
        h2 = s_condVariance(model, misspType, r, d, theta, omega, alpha, beta, gamma);        
        
        w1 = (mod(d,5) == 1);
        w2 = (mod(d,5) == 2);
        w3 = (mod(d,5) == 3);
        w4 = (mod(d,5) == 4);
        w5 = (mod(d,5) == 0);
        
        w = [w1, w2, w3, w4, w5];
        
        return_value = sum( -log(h2(:,1)) - (r(:,1)  - w*theta').^2./h2(:,1));
    end
    
end

