function [r, h2] = s_mySimulate(model, misspType, mu_, param, omega, alpha, beta, gamma, T)
% T - число наблюдений в симул€ции

if (strcmp(model, 'gjr') == 1)
    sigma2 = omega/(1 - alpha - beta - gamma/2);
    h2 = sigma2*ones(T, 1);
    z = normrnd(0, 1, T, 1);
    r = zeros(T,1);

    if (misspType == 0)

        rho = param;    
        r(1,1) = mu_ + rho*mu_/(1-rho) + h2(1,1)^(1/2)*z(1,1);

        for i = 2:length(r)
            h2(i,1) = omega + alpha*(h2(i-1,1)^(1/2)*z(i-1,1))^2 + beta*h2(i-1,1) + gamma*(z(i-1,1)<0);
        end
        r(2:end,1) = mu_ + rho*r(1:end-1,1) + h2(2:end,1).^(1/2).*z(2:end,1);

    elseif (misspType == 1)

        delta = param;
        
        for i = 2:length(r)
            h2(i,1) = omega + alpha*(h2(i-1,1)^(1/2)*z(i-1,1))^2 + beta*h2(i-1,1) + gamma*(z(i-1,1)<0);
        end
        r = mu_ + delta*h2 + h2.^(1/2).*z;


    elseif (misspType == 2)    
        d = ones(T, 1);
        for i = 1:T
            d(i, 1) = mod(i, 5);
        end

        w1 = (mod(d,5) == 1);
        w2 = (mod(d,5) == 2);
        w3 = (mod(d,5) == 3);
        w4 = (mod(d,5) == 4);
        w5 = (mod(d,5) == 0);
        
        w = [w1, w2, w3, w4, w5];       

        for i = 2:length(r)
            h2(i,1) = omega + alpha*(h2(i-1,1)^(1/2)*z(1,1))^2 + beta*h2(i-1,1) + gamma*(z(i-1,1)<0);
        end
        r = mu_ + param*w' + h2.^(1/2).*z;
    else
        error('misspType belongs to {0, 1, 2}');
    end
    
elseif strcmp(model, 'egarch') == 1
    sigma2 = exp((omega + alpha/2 + gamma)/(1 - beta));    
    h2 = sigma2*ones(T,1);
    log_h2 = log(sigma2)*ones(T,1);
    z = normrnd(0, 1, [T, 1]);
    r = zeros(T,1);

    if (misspType == 0)

        rho = param;    
        r(1,1) = mu_ + rho*mu_/(1-rho) + h2(1,1)^(1/2)*z(1,1);

        for i = 2:length(r)
            log_h2(i,1) = omega + alpha*abs(z(i-1,1)) + beta*log_h2(i-1,1) + gamma*z(i-1,1);
        end
        h2 = exp(log_h2);
        r(2:end,1) = mu_ + rho*r(1:end-1,1) + h2(2:end,1).^(1/2).*z(2:end,1);
        
    elseif (misspType == 1)

        delta = param;

        for i = 2:length(r)
            log_h2(i,1) = omega + alpha*abs(z(i-1,1)) + beta*log_h2(i-1,1) + gamma*z(i-1,1);
        end
        h2 = exp(log_h2);
        r = mu_ + delta*h2 + h2.^(1/2).*z;

    elseif (misspType == 2)    
        d = ones(T, 1);
        for i = 1:T
            d(i, 1) = mod(i, 5);
            if mod(i, 5) == 0
                d(i, 1) = 5;
            end
        end

        w1 = (mod(d,5) == 1);
        w2 = (mod(d,5) == 2);
        w3 = (mod(d,5) == 3);
        w4 = (mod(d,5) == 4);
        w5 = (mod(d,5) == 0);
        
        w = [w1, w2, w3, w4, w5];
       

        for i = 2:length(r)
            log_h2(i,1) = omega + alpha*abs(z(i-1,1)) + beta*log_h2(i-1,1) + gamma*z(i-1,1);
        end
        h2 = exp(log_h2);
        r = w*param' + h2.^(1/2).*z;
    else
        error('misspType belongs to {0, 1, 2}');
    end
    
else
    error('choose an appropriate model');
end