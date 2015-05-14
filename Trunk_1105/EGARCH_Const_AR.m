classdef EGARCH_Const_AR < EGARCH_BaseModel
    properties
        mu, mu0, rho, rho0
    end

    methods (Access = public)
        function obj = EGARCH_Const_AR(mu0, rho0, omega0, alpha0, beta0, gamma0) % конструктор
            obj.mu0    = mu0;
            obj.rho0   = rho0;
            obj.omega0 = omega0;
            obj.alpha0 = alpha0;
            obj.beta0  =  beta0;
            obj.gamma0 = gamma0;
            obj.mu     = mu0;
            obj.rho    = rho0;
            obj.omega  = omega0;
            obj.alpha  = alpha0;
            obj.beta   =  beta0;
            obj.gamma  = gamma0;
        end
        
        function MyDisplay(self)            
            fprintf('%6s %12s\n\n', 'model', 'EGARCH, Const + AR');
            
            fprintf('%6s %12s %12s\r','param','true', 'estimated');
            fprintf('%6s %12.6f %12.6f\n', 'mu',    self.mu0,    self.mu);
            fprintf('%6s %12.6f %12.6f\n', 'rho',   self.rho0,   self.rho);
            fprintf('%6s %12.6f %12.6f\n', 'omega', self.omega0, self.omega);
            fprintf('%6s %12.6f %12.6f\n', 'alpha', self.alpha0, self.alpha);
            fprintf('%6s %12.6f %12.6f\n', 'beta',  self.beta0,  self.beta);
            fprintf('%6s %12.6f %12.6f\n', 'gamma', self.gamma0, self.gamma);
        end 
        
        function SimData = Simulate(self, T)
            sigma2 = exp((self.omega0 + self.alpha0*sqrt(2/pi))/(1 - self.beta0));
            z      = normrnd(0, 1, T, 1);
            h2     = sigma2*ones(T,1);
            log_h2 = log(sigma2)*ones(T,1);
            
            SimData = zeros(T,1);
            SimData(1,1) = self.mu0 + self.rho0*self.mu0/(1 - self.rho0) + sqrt(h2(1,1))*z(1,1);
                        
            for i = 2:T
                log_h2(i,1) = self.omega0 + self.alpha0*abs(z(i-1,1)) + ...
                              self.beta0*log_h2(i-1,1) + self.gamma0*z(i-1,1);
                h2(i,1) = exp(log_h2(i,1));
                SimData(i,1) = self.mu0 + self.rho0*SimData(i-1,1) + sqrt(h2(i,1))*z(i,1);
            end
            self.sigma2 = h2;
            
        end
        
        function [loss, loss2, VaR_exceeded] = Predict(self, p)
            [h2, e] = self.CondVar();
            log_h2_pred = self.omega + ...
                        self.alpha*( abs(e(end,1))/sqrt(h2(end,1)) ) + ...
                        self.beta*log(h2(end,1)) + ...
                        self.gamma*e(end,1)/sqrt(h2(end,1));
                    
            h2_pred = exp(log_h2_pred);
            
            VaR = self.mu + self.rho*self.data(end,1) + sqrt(h2_pred)*norminv(p,0,1);
            VaR_exceeded = (VaR > self.data_plus(end,1));
            
            loss  = QLIKE(self.sigma2(end,1), h2_pred);
            loss2 = QLIKE2(self.sigma2(end,1), h2_pred);
        end
       
    end

    methods (Access = protected)
        function x0 = GetX0 (self) 
            x0 = [self.mu0, self.rho0, self.omega0, self.alpha0, self.beta0, self.gamma0];
        end

        function Install (self, x, data)
            self.data  = data;
            self.mu    = x(1);
            self.rho   = x(2);
            self.omega = x(3);
            self.alpha = x(4);
            self.beta  = x(5);
            self.gamma = x(6);
        end
        
        function log_ret = LogL(self)
            h2 = self.CondVar();
            log_ret = sum( -0.5*(log(2*pi) + log(h2(2:end,1)) +...
                (self.data(2:end,1) - self.mu - self.rho*self.data(1:end-1, 1)).^2./h2(2:end,1)) );
        end
        
        function size = Param0size(self)
            size = [1, 2];
        end
        
        function [h2, e] = CondVar(self) 
            sigma2 = exp((self.omega0 + self.alpha0*sqrt(2/pi))/(1 - self.beta0));
            log_h2 = log(sigma2)*ones(size(self.data));
            h2 = sigma2*ones(size(self.data));
            e = zeros(size(self.data));
            
            e(1,1)     = self.data(1,1)     - self.mu - self.rho*mean(self.data);
            e(2:end,1) = self.data(2:end,1) - self.mu - self.rho*self.data(1:end-1,1);            
           
            for i = 2:length(self.data)     
                log_h2(i,1) = self.omega + self.alpha*( abs(e(i-1,1))/sqrt(h2(i-1,1))) + ...
                              self.beta*log_h2(i-1,1) + self.gamma*e(i-1,1)/sqrt(h2(i-1,1));
                h2(i,1) = exp(log_h2(i,1));
            end   
        end      
        
    end
end