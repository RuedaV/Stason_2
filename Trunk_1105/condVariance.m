function h2 = condVariance(model, misspType, r, d, param, omega, alpha, beta, gamma)
% misspType belongs to {0, 1, 2, 3}
% misspType = 0 for mu
% misspType = 1 for rho
% misspType = 2 for delta
% misspType = 3 for theta


if strcmp(model, 'gjr') == 1
    sigma2 = mean(r(1:50,1).^2);
    h2 = sigma2*ones(size(r));
    e = zeros(size(r));
    
    if (misspType == 0)   

        mu_ = param;
        e = r - mu_;
        i = 2;
        while i <= length(e)      
            h2(i,1) = omega + alpha*e(i-1,1)^2 + beta*h2(i-1,1) + gamma*(e(i-1,1)<0);
            i = i + 1;
        end
        
    elseif (misspType == 1)

        rho = param;
        e(1,1) = r(1,1) - rho*mean(r(1:5,1));
        e(2:end, 1) = r(2:end, 1) - rho*r(1:end-1, 1);
        i = 2;

        while i <= length(e)      
            h2(i,1) = omega + alpha*e(i-1,1)^2 + beta*h2(i-1,1) + gamma*(e(i-1,1)<0);
            i = i + 1;
        end
        
    elseif (misspType == 2)

        delta = param;
        e(1,1) = r(1,1) - delta*sigma2;
        i = 2;

        while i <= length(e)      
            h2(i,1) = omega + alpha*e(i-1,1)^2 + beta*h2(i-1,1) + gamma*(e(i-1,1)<0);
            e(i,1) = r(i,1) - delta*h2(i-1,1);
            i = i + 1;
        end
        
    elseif (misspType == 3)  
        
        w1 = (mod(d,5) == 1);
        w2 = (mod(d,5) == 2);
        w3 = (mod(d,5) == 3);
        w4 = (mod(d,5) == 4);
        w5 = (mod(d,5) == 0);
        
        w = [w1, w2, w3, w4, w5];
        
        h2 = sigma2*ones(size(r));
        e = r - w*param';
        i = 2;

        while i <= length(e)
            h2(i,1) = omega + alpha*e(i-1,1)^2 + beta*h2(i-1,1) + gamma*(e(i-1,1)<0);
            i = i + 1;
        end
    end   
        
elseif (strcmp(model, 'egarch') == 1)
    sigma2 = mean(r(1:50,1).^2);
    log_h2 = log(sigma2)*ones(size(r));
    h2 = sigma2*ones(size(r));    
    e = zeros(size(r));

    if (misspType == 0)   

        mu_ = param;
        e = r - mu_;
        i = 2;

        while i <= length(e)      
            log_h2(i,1) = omega + alpha*abs(e(i-1,1))/sqrt(h2(i-1,1)) + beta*log_h2(i-1,1) + gamma*e(i-1,1)/sqrt(h2(i-1,1));
            h2(i,1) = exp(log_h2(i,1));
            i = i + 1;
        end
        
    elseif (misspType == 1)

        rho = param;
        e(1,1) = r(1,1) - rho*mean(r(1:5,1));
        e(2:end, 1) = r(2:end, 1) - rho*r(1:end-1, 1);
        i = 2;

        while i <= length(e)      
            log_h2(i,1) = omega + alpha*abs(e(i-1,1))/sqrt(h2(i-1,1)) + beta*log_h2(i-1,1) + gamma*e(i-1,1)/sqrt(h2(i-1,1));
            h2(i,1) = exp(log_h2(i,1));
            i = i + 1;
        end
        
    elseif (misspType == 2)

        delta = param;
        e(1,1) = r(1,1) - delta*sigma2;
        i = 2;

        while i <= length(e)      
            log_h2(i,1) = omega + alpha*abs(e(i-1,1))/sqrt(h2(i-1,1)) + beta*log_h2(i-1,1) + gamma*e(i-1,1)/sqrt(h2(i-1,1));
            h2(i,1) = exp(log_h2(i,1));
            e(i,1) = r(i,1) - delta*log_h2(i-1,1);
            i = i + 1;
        end
        
    elseif (misspType == 3)  

        w1 = (mod(d,5) == 1);
        w2 = (mod(d,5) == 2);
        w3 = (mod(d,5) == 3);
        w4 = (mod(d,5) == 4);
        w5 = (mod(d,5) == 0);
        
        w = [w1, w2, w3, w4, w5];
        
        h2 = sigma2*ones(size(r));
        e = r - w*param';
        i = 2;


        while i <= length(e)            
            log_h2(i,1) = omega + alpha*abs(e(i-1,1))/sqrt(h2(i-1,1)) + beta*log_h2(i-1,1) + gamma*e(i-1,1)/sqrt(h2(i-1,1));
            h2(i,1) = exp(log_h2(i,1));
            i = i + 1;
        end
    end
    
else
    error('Choose an appropriate model!')
end
