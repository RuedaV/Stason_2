classdef EGARCH_Ssn < EGARCH_BaseModel
    properties        
        day, day_plus, 
        theta1,  theta2,  theta3,  theta4,  theta5, 
        theta10, theta20, theta30, theta40, theta50
    end

    methods (Access = public)
        function obj = EGARCH_Ssn(theta10, theta20, theta30, theta40, theta50,...
                omega0, alpha0, beta0, gamma0) % конструктор
            obj.omega0  = omega0;
            obj.alpha0  = alpha0;
            obj.beta0   =  beta0;
            obj.gamma0  = gamma0;
            obj.theta10 = theta10;
            obj.theta20 = theta20;
            obj.theta30 = theta30;
            obj.theta40 = theta40;
            obj.theta50 = theta50;
            
            obj.omega   = omega0;
            obj.alpha   = alpha0;
            obj.beta    =  beta0;
            obj.gamma   = gamma0;
            obj.theta1  = theta10;
            obj.theta2  = theta20;
            obj.theta3  = theta30;
            obj.theta4  = theta40;
            obj.theta5  = theta50;
        end
        
        function MyDisplay(self)            
            fprintf('%6s %12s\n\n', 'model', 'EGARCH, Ssn');
            
            fprintf('%6s %12s %12s\r','param','true', 'estimated');
            fprintf('%6s %12.6f %12.6f\n', 'theta1', self.theta10, self.theta1);
            fprintf('%6s %12.6f %12.6f\n', 'theta2', self.theta20, self.theta2);
            fprintf('%6s %12.6f %12.6f\n', 'theta3', self.theta30, self.theta3);
            fprintf('%6s %12.6f %12.6f\n', 'theta4', self.theta40, self.theta4);
            fprintf('%6s %12.6f %12.6f\n', 'theta5', self.theta50, self.theta5);
            fprintf('%6s %12.6f %12.6f\n', 'omega',  self.omega0,  self.omega);
            fprintf('%6s %12.6f %12.6f\n', 'alpha',  self.alpha0,  self.alpha);
            fprintf('%6s %12.6f %12.6f\n', 'beta',   self.beta0,   self.beta);
            fprintf('%6s %12.6f %12.6f\n', 'gamma',  self.gamma0,  self.gamma);
        end 
        
        function [SimData, day] = Simulate(self, T)
            day = zeros(T, 1);
            for i = 1:T
                day(i,1) = mod(i,5);
                if day(i,1) == 0
                    day(i,1) = 5;
                end
            end
            self.day = day;
            
            w1 = (mod(day, 5) == 1);
            w2 = (mod(day, 5) == 2);
            w3 = (mod(day, 5) == 3);
            w4 = (mod(day, 5) == 4);
            w5 = (mod(day, 5) == 0);
            w = [w1, w2, w3, w4, w5];
            
            sigma2 = exp((self.omega0 + self.alpha0*sqrt(2/pi))/(1 - self.beta0));
            z      = normrnd(0, 1, T, 1);
            h2     = sigma2*ones(T,1);
            log_h2 = log(sigma2)*ones(T,1);            
            
            for i = 2:T
                log_h2(i,1) = self.omega0 + self.alpha0*abs(z(i-1,1)) + ...
                              self.beta0*log_h2(i-1,1) + self.gamma0*z(i-1,1);
                h2(i,1) = exp(log_h2(i,1));
            end
            
            SimData =  w*[self.theta10, self.theta20, self.theta30, self.theta40, self.theta50]' + h2.^(0.5).*z;
            self.sigma2 = h2;
        end
        
        
        function setDaily (self, d)                           
            self.day = d;
        end
        
        function [loss, VaR_true, VaR_pred] = Predict(self, p)
            [h2, e] = self.CondVar();
            log_h2_pred = self.omega + ...
                        self.alpha*( abs(e(end,1))/sqrt(h2(end,1)) ) + ...
                        self.beta*log(h2(end,1)) + ...
                        self.gamma*e(end,1)/sqrt(h2(end,1));
                    
            h2_pred = exp(log_h2_pred);
            
            w1 = (mod(self.day, 5) == 1);
            w2 = (mod(self.day, 5) == 2);
            w3 = (mod(self.day, 5) == 3);
            w4 = (mod(self.day, 5) == 4);
            w5 = (mod(self.day, 5) == 0);
            
            VaR_pred = self.theta1*w1(end,1)...
            + self.theta2*w2(end,1) + self.theta3*w3(end,1) + self.theta4*w4(end,1) + self.theta5*w5(end,1) ...
            + sqrt(h2_pred)*norminv(p,0,1);
        
            VaR_true = self.theta10*w1(end,1)...
            + self.theta20*w2(end,1) + self.theta30*w3(end,1) + self.theta40*w4(end,1) + self.theta50*w5(end,1) ...
            + sqrt(self.sigma2(end,1))*norminv(p,0,1);
            
            loss  = QLIKE(self.sigma2(end,1), h2_pred);
        end
        
        
    end

    methods (Access = protected)
        function x0 = GetX0 (self) 
            x0 = [self.theta10, self.theta20, self.theta30, self.theta40, self.theta50,...
                self.omega0, self.alpha0, self.beta0, self.gamma0];
        end

        function Install (self, x, data)
            self.data   = data;
            self.theta1 = x(1);
            self.theta2 = x(2);
            self.theta3 = x(3);
            self.theta4 = x(4);
            self.theta5 = x(5);
            self.omega  = x(6);
            self.alpha  = x(7);
            self.beta   = x(8);
            self.gamma  = x(9);
        end
        
        function log_ret = LogL(self)
            h2 = self.CondVar();
            w1 = (mod(self.day, 5) == 1);
            w2 = (mod(self.day, 5) == 2);
            w3 = (mod(self.day, 5) == 3);
            w4 = (mod(self.day, 5) == 4);
            w5 = (mod(self.day, 5) == 0);
            w = [w1, w2, w3, w4, w5];
            delta_new = [self.theta1, self.theta2, self.theta3, self.theta4, self.theta5];
            
            log_ret = sum( -0.5*(log(2*pi) + log(h2) +...
                (self.data - w*delta_new').^2./h2) );
        end
        
        function size = Param0size(self)
            size = [1, 5];
        end
        
        function [h2, e] = CondVar(self) 
%             sigma2 = self.omega/(1 - self.alpha - self.beta - self.gamma/2);
            w1 = (mod(self.day, 5) == 1);
            w2 = (mod(self.day, 5) == 2);
            w3 = (mod(self.day, 5) == 3);
            w4 = (mod(self.day, 5) == 4);
            w5 = (mod(self.day, 5) == 0);
            
            w = [w1, w2, w3, w4, w5];
            delta_new = [self.theta1, self.theta2, self.theta3, self.theta4, self.theta5];
            
            sigma2 = sum((self.data - w*delta_new' - mean(self.data - w*delta_new')).^2)/(length(self.data) - 1);
            log_h2 = log(sigma2)*ones(size(self.data));
            h2 = sigma2*ones(size(self.data));
            
            e = self.data - w*delta_new';
            
            for i = 2:length(self.data)      
                log_h2(i,1) = self.omega + ...
                        self.alpha*( abs(e(i-1,1))/sqrt(h2(i-1,1)) ) + ...
                        self.beta*log_h2(i-1,1) + ...
                        self.gamma*e(i-1,1)/sqrt(h2(i-1,1));
                    
                h2(i,1) = exp(log_h2(i,1));
            end   
        end       
        
    end
end