function [return_value] = logL(model, misspType, r, d, x)
% misspType belongs to {0, 1, 2, 3}
% misspType = 0 for mu
% misspType = 1 for rho
% misspType = 2 for delta
% misspType = 3 for theta 

    
    if (misspType == 0)      
        mu_   = x(1,1);
        omega = x(1,2);
        alpha = x(1,3);
        beta  = x(1,4);
        gamma = x(1,5);
        h2 = condVariance(model, misspType, r, d, mu_, omega, alpha, beta, gamma);     
        return_value = sum( -log(h2(:,1)) - (r(:,1) - mu_).^2./h2(:,1));
        
    elseif (misspType == 1)       
        rho = x(1,1);
        omega = x(1,2);
        alpha = x(1,3);
        beta  = x(1,4);
        gamma = x(1,5);
        h2 = condVariance(model, misspType, r, d, rho, omega, alpha, beta, gamma);    
        return_value = sum( -log(h2(2:end,1)) - (r(2:end,1) - rho*r(1:end-1,1)).^2./h2(2:end,1));
        
    elseif (misspType == 2)      
        delta = x(1,1);
        omega = x(1,2);
        alpha = x(1,3);
        beta  = x(1,4);
        gamma = x(1,5);
        h2 = condVariance(model, misspType, r, d, delta, omega, alpha, beta, gamma);        
        
        return_value = sum( -log(h2(:,1)) - (r(:,1) - delta*h2(:,1)).^2./h2(:,1));   
    elseif (misspType == 3)
              
        theta = x(1,1:5);
        omega = x(1,6);
        alpha = x(1,7);
        beta  = x(1,8);
        gamma = x(1,9);
        h2 = condVariance(model, misspType, r, d, theta, omega, alpha, beta, gamma);        
        
        w1 = (mod(d,5) == 1);
        w2 = (mod(d,5) == 2);
        w3 = (mod(d,5) == 3);
        w4 = (mod(d,5) == 4);
        w5 = (mod(d,5) == 0);
        
        w = [w1, w2, w3, w4, w5];
        
        return_value = sum( -log(h2(:,1)) - (r(:,1) - w*theta').^2./h2(:,1) );
    end
    
end

