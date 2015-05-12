classdef GJR_Ssn < GJR_BaseModel
    properties
        day, day_plus,
        theta1,  theta2,  theta3,  theta4,  theta5, 
        theta10, theta20, theta30, theta40, theta50
    end

    methods (Access = public)
        function obj = GJR_Ssn( theta10, theta20, theta30, ...
            theta40, theta50, omega0, alpha0, beta0, gamma0) % конструктор
            %obj.data   = data;
            %ƒату запихиваем ручками   
            %obj.day     = day;
            
            obj.theta10 = theta10;
            obj.theta20 = theta20;
            obj.theta30 = theta30;
            obj.theta40 = theta40;
            obj.theta50 = theta50;                
            obj.omega0  = omega0;
            obj.alpha0  = alpha0;
            obj.beta0   =  beta0;
            obj.gamma0  = gamma0;
            
            obj.theta1  = theta10;
            obj.theta2  = theta20;
            obj.theta3  = theta30;
            obj.theta4  = theta40;
            obj.theta5  = theta50; 
            obj.omega   = omega0;
            obj.alpha   = alpha0;
            obj.beta    =  beta0;
            obj.gamma   = gamma0;
        end
        
        function MyDisplay(self)            
            fprintf('%6s %12s\n\n', 'model', 'GJR, Seasonality');
            
            fprintf('%6s %12s %12s\r','param','true', 'estimated');
            fprintf('%6s %12.6f %12.6f\n', 'theta1', self.theta10, self.theta1);
            fprintf('%6s %12.6f %12.6f\n', 'theta2', self.theta20, self.theta2);
            fprintf('%6s %12.6f %12.6f\n', 'theta3', self.theta30, self.theta3);
            fprintf('%6s %12.6f %12.6f\n', 'theta4', self.theta40, self.theta4);
            fprintf('%6s %12.6f %12.6f\n', 'theta5', self.theta50, self.theta5);
            fprintf('%6s %12.6f %12.6f\n', 'omega', self.omega0, self.omega);
            fprintf('%6s %12.6f %12.6f\n', 'alpha', self.alpha0, self.alpha);
            fprintf('%6s %12.6f %12.6f\n', 'beta',  self.beta0,  self.beta);
            fprintf('%6s %12.6f %12.6f\n', 'gamma', self.gamma0, self.gamma);
        end 
        
        function res = Compare(self, x, rival)
            if sum(x(1:5).^2) < sum(rival(1:5).^2)
                res = 1;
            else
                res = 0;
            end                
        end   
          
        function [SimData, day] = Simulate(self, T)
            day = zeros(T, 1);
            for i = 1:T
                day(i,1) = mod(i,5);
                if day(i,1) == 0
                    day(i,1) = 5;
                end
            end
            
            w1 = (mod(day, 1) == 1);
            w2 = (mod(day, 1) == 2);
            w3 = (mod(day, 1) == 3);
            w4 = (mod(day, 1) == 4);
            w5 = (mod(day, 1) == 0);
            w = [w1, w2, w3, w4, w5];
            
            sigma2 = self.omega0/(1 - self.alpha0 - self.beta0 - self.gamma0/2);
            z = normrnd(0, 1, T, 1);
            h2 = sigma2*ones(T,1);
            
            for i = 2:T
                h2(i,1) = self.omega0 + self.alpha0*h2(i-1,1)*z(i-1,1)^2 +...
                    self.beta0*h2(i-1,1) + self.gamma0*h2(i-1,1)*z(i-1,1)^2*(z(i-1,1) < 0);
            end     
            
            SimData = w*[self.theta10, self.theta20, self.theta30, self.theta40, self.theta50]' + sqrt(h2).*z;
            
        end
        
        function setDaily (self, d)                           
            self.day = d;
        end
        
        function [loss, loss2, VaR_exceeded] = Predict(self)
            [h2, e] = self.CondVar();
            h2_pred = self.omega + self.alpha*e(end,1)^2 ...
                    + self.beta*h2(end,1) + self.gamma*(e(end,1)<0)*e(end,1)^2;
                
            VaR = (self.data(end, 1) - e(end, 1)) + sqrt(h2(end,1))*norminv(0.05,0,1);
            VaR_exceeded = (VaR > self.data_plus(end,1));
            
            day_tmp = self.day;
            data_temp = self.data;
            self.data = self.data_plus;
            self.day  = self.day_plus;
            self.Switch();
            [h2_plus, e_plus] = self.CondVar();
            h2_proxy = h2_plus(end, 1);
            self.Switch();
            
            self.data = data_temp;
            self.day = day_tmp;
            loss = QLIKE(h2_proxy, h2_pred);
            loss2 = QLIKE2(h2_proxy, h2_pred);
        end
        
    end

    methods (Access = protected)
        
        function x0 = GetX0 (self) 
            x0 = [self.theta10, self.theta20, self.theta30, self.theta40, self.theta50, ...
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
        
        function log_ret = LogL (self)
            w1 = (mod(self.day,5) == 1);
            w2 = (mod(self.day,5) == 2);
            w3 = (mod(self.day,5) == 3);
            w4 = (mod(self.day,5) == 4);
            w5 = (mod(self.day,5) == 0);
            w = [w1, w2, w3, w4, w5];
            
            h2 = self.CondVar();
            log_ret = sum( -0.5*(log(2*pi) + log(h2(:,1)) +...
                (self.data(:,1) - w*[self.theta1, self.theta2, self.theta3, self.theta4, self.theta5]').^2./h2(:,1)));
        end
        
        function size = Param0size (self)
            size = [1, 5];
        end
        
        function [h2, e] = CondVar(self) 
            sigma2 = sum((self.data - mean(self.data)).^2)/(length(self.data) - 1);
            h2 = sigma2*ones(size(self.data));
            w1 = (mod(self.day,5) == 1);
            w2 = (mod(self.day,5) == 2);
            w3 = (mod(self.day,5) == 3);
            w4 = (mod(self.day,5) == 4);
            w5 = (mod(self.day,5) == 0);
            w = [w1, w2, w3, w4, w5];
            e = self.data - w*[self.theta1, self.theta2, self.theta3, self.theta4, self.theta5]';
            
            for i = 2:length(self.data)      
                h2(i,1) = self.omega + self.alpha*e(i-1,1)^2 ...
                    + self.beta*h2(i-1,1) + self.gamma*(e(i-1,1)<0)*e(i-1,1)^2;
            end   
        end       
        
         function Switch(self)
                [self.theta1, self.theta10]  =  self.Switch_pair(self.theta1, self.theta10);
                [self.theta2, self.theta20]  =  self.Switch_pair(self.theta2, self.theta20);
                [self.theta3, self.theta30]  =  self.Switch_pair(self.theta3, self.theta30);
                [self.theta4, self.theta40]  =  self.Switch_pair(self.theta4, self.theta40);
                [self.theta5, self.theta50]  =  self.Switch_pair(self.theta5, self.theta50);            
                [self.omega,  self.omega0]   =  self.Switch_pair(self.omega,  self.omega0);
                [self.alpha,  self.alpha0]   =  self.Switch_pair(self.alpha,  self.alpha0);
                [self.beta,   self.beta0]    =  self.Switch_pair(self.beta,   self.beta0);
                [self.gamma,  self.gamma0]   =  self.Switch_pair(self.gamma,  self.gamma0);
        end
    end
end