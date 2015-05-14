classdef EGARCH_Var < EGARCH_BaseModel
    properties
        delta, delta0
    end

    methods (Access = public)
        function obj = EGARCH_Var(delta0, omega0, alpha0, beta0, gamma0) % конструктор
            %obj.data   = data;
            %ƒату запихиваем ручками
            obj.delta0 = delta0;
            obj.omega0 = omega0;
            obj.alpha0 = alpha0;
            obj.beta0  = beta0;      
            obj.gamma0 = gamma0;
            
            obj.delta  = delta0;
            obj.omega  = omega0;
            obj.alpha  = alpha0;
            obj.beta   = beta0;
            obj.gamma  = gamma0;
        end
        
        function MyDisplay(self)            
            fprintf('%6s %12s\n\n', 'model', 'EGARCH, Variance-in-Mean');
            
            fprintf('%6s %12s %12s\r','param','true', 'estimated');
            fprintf('%6s %12.6f %12.6f\n', 'delta', self.delta0, self.delta);
            fprintf('%6s %12.6f %12.6f\n', 'omega', self.omega0, self.omega);
            fprintf('%6s %12.6f %12.6f\n', 'alpha', self.alpha0, self.alpha);
            fprintf('%6s %12.6f %12.6f\n', 'beta',  self.beta0,  self.beta);
            fprintf('%6s %12.6f %12.6f\n', 'gamma', self.gamma0, self.gamma);
        end 
        
        function [loss, VaR_true, VaR_pred] = Predict(self, p)
            [h2, e] = self.CondVar();
            log_h2_pred = self.omega + ...
                        self.alpha*( abs(e(end,1))/sqrt(h2(end,1)) ) + ...
                        self.beta*log(h2(end,1)) + ...
                        self.gamma*e(end,1)/sqrt(h2(end,1));
                    
            h2_pred = exp(log_h2_pred);
            

            VaR_pred = self.delta*h2_pred + sqrt(h2_pred)*norminv(0.05,0,1);
            VaR_true = self.delta0*h2_pred + sqrt(self.sigma2(end,1))*norminv(0.05,0,1);
           
            
            loss  = QLIKE(self.sigma2(end,1), h2_pred);
        end
        
        
    end

    methods (Access = protected)        
        function x0 = GetX0 (self) 
            x0 = [self.delta0, self.omega0, self.alpha0, self.beta0, self.gamma0];
        end

        function Install (self, x, data)
            self.data  = data;
            self.delta = x(1);
            self.omega = x(2);
            self.alpha = x(3);
            self.beta  = x(4);
            self.gamma = x(5);
        end
        
        function log_ret = LogL (self)
            h2 = self.CondVar();
            log_ret = sum( -0.5*(log(2*pi) + log(h2(:,1)) + ...
                (self.data(:,1) - self.delta*h2(:,1)).^2./(h2(:,1))));
        end
        
        function size = Param0size (self)
            size = 1;
        end
        
        function [h2, e] = CondVar(self) 
            sigma2 = sum((self.data - mean(self.data)).^2)/(length(self.data) - 1);
            log_h2 = log(sigma2)*ones(size(self.data));
            h2 = sigma2*ones(size(self.data));
            e = zeros(size(self.data));
            e(1,1) = self.data(1,1) - self.delta*h2(1,1);
            
            for i = 2:length(self.data)     
                log_h2(i,1) = self.omega + ...
                        self.alpha*( abs(e(i-1,1))/sqrt(h2(i-1,1)) ) + ...
                        self.beta*log_h2(i-1,1) + ...
                        self.gamma*e(i-1,1)/sqrt(h2(i-1,1));
                    
                h2(i,1) = exp(log_h2(i,1));
                e(i,1) = self.data(i,1) - self.delta*h2(i,1);
            end  
        end       
 
    end
end















