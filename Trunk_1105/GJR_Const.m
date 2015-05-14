classdef GJR_Const < GJR_BaseModel
    properties
        mu, mu0
    end

    methods (Access = public)
        function obj = GJR_Const(mu0, omega0, alpha0, beta0, gamma0) % конструктор
            obj.mu0    = mu0;
            obj.omega0 = omega0;
            obj.alpha0 = alpha0;
            obj.beta0  =  beta0;
            obj.gamma0 = gamma0;
            obj.mu     = mu0;
            obj.omega  = omega0;
            obj.alpha  = alpha0;
            obj.beta   =  beta0;
            obj.gamma  = gamma0;
        end
        
        function MyDisplay(self)            
            fprintf('%6s %12s\n\n', 'model', 'GJR, Const');
            
            fprintf('%6s %12s %12s\r','param','true', 'estimated');
            fprintf('%6s %12.6f %12.6f\n', 'mu',    self.mu0,    self.mu);
            fprintf('%6s %12.6f %12.6f\n', 'omega', self.omega0, self.omega);
            fprintf('%6s %12.6f %12.6f\n', 'alpha', self.alpha0, self.alpha);
            fprintf('%6s %12.6f %12.6f\n', 'beta',  self.beta0,  self.beta);
            fprintf('%6s %12.6f %12.6f\n', 'gamma', self.gamma0, self.gamma);
        end 
        
         function [loss, loss2, VaR_exceeded] = Predict(self, p)
            [h2, e] = self.CondVar();
            h2_pred = self.omega + self.alpha*e(end,1)^2 ...
                    + self.beta*h2(end,1) + self.gamma*(e(end,1)<0)*e(end,1)^2;
                
            VaR = self.mu + sqrt(h2_pred)*norminv(p,0,1);
            VaR_exceeded = (VaR > self.data_plus(end,1));
            
            loss  = QLIKE(self.sigma2(end,1), h2_pred);
            loss2 = QLIKE2(self.sigma2(end,1), h2_pred);
        end
        
    end

    methods (Access = protected)
        function x0 = GetX0 (self) 
            x0 = [self.mu0, self.omega0, self.alpha0, self.beta0, self.gamma0];
        end

        function Install (self, x, data)
            self.data  = data;
            self.mu    = x(1);
            self.omega = x(2);
            self.alpha = x(3);
            self.beta  = x(4);
            self.gamma = x(5);
        end
        
        function log_ret = LogL (self)
            h2 = self.CondVar();
            log_ret = sum( -0.5*(log(2*pi) + log(h2(:,1)) + (self.data(:,1) - self.mu).^2./h2(:,1)));
        end
        
        function size = Param0size (self)
            size = 1;
        end
        
        function [h2, e] = CondVar(self) 
            sigma2 = sum((self.data - mean(self.data)).^2)/(length(self.data) - 1);
            h2 = sigma2*ones(size(self.data));
            e = self.data - self.mu;
            
            for i = 2:length(self.data)      
                h2(i,1) = self.omega + self.alpha*e(i-1,1)^2 ...
                    + self.beta*h2(i-1,1) + self.gamma*(e(i-1,1)<0)*e(i-1,1)^2;
            end   
        end    
        
    end
end