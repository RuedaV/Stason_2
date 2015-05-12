classdef GJR_Const_AR < GJR_BaseModel
    properties
        mu, mu0, rho, rho0
    end

    methods (Access = public)
        function obj = GJR_Const_AR(mu0, rho0, omega0, alpha0, beta0, gamma0) % конструктор
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
            fprintf('%6s %12s\n\n', 'model', 'GJR, Const + AR');
            
            fprintf('%6s %12s %12s\r','param','true', 'estimated');
            fprintf('%6s %12.6f %12.6f\n', 'mu',    self.mu0,    self.mu);
            fprintf('%6s %12.6f %12.6f\n', 'rho',   self.rho0,   self.rho);
            fprintf('%6s %12.6f %12.6f\n', 'omega', self.omega0, self.omega);
            fprintf('%6s %12.6f %12.6f\n', 'alpha', self.alpha0, self.alpha);
            fprintf('%6s %12.6f %12.6f\n', 'beta',  self.beta0,  self.beta);
            fprintf('%6s %12.6f %12.6f\n', 'gamma', self.gamma0, self.gamma);
        end 
        
        function SimData = Simulate(self, T)
            sigma2 = self.omega0/(1 - self.alpha0 - self.beta0 - self.gamma0/2);
            z = normrnd(0, 1, T, 1);
            h2 = sigma2*ones(T,1);
            SimData = zeros(T,1);
            SimData(1,1) = self.mu0 + self.rho0*self.mu0/(1 - self.rho0) + sqrt(h2(1,1))*z(1,1);
            
            for i = 2:T
                h2(i,1) = self.omega0 + self.alpha0*h2(i-1,1)*z(i-1,1)^2 +...
                    self.beta0*h2(i-1,1) + self.gamma0*h2(i-1,1)*z(i-1,1)^2*(z(i-1,1) < 0);
                SimData(i,1) = self.mu0 + self.rho0*SimData(i-1,1) + sqrt(h2(i,1))*z(i,1);
            end
            
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
%                 -0.5*(log(2*pi) + log(h2(1,1)) + (self.data(1,1) - self.mu - self.rho*mean(self.data(:, 1))).^2/h2(1,1));
        end
        
        function size = Param0size(self)
            size = [1, 2];
        end
        
        function [h2, e] = CondVar(self) 
%             sigma2 = self.omega/(1 - self.alpha - self.beta - self.gamma/2);
            sigma2 = (1 - self.rho^2)*sum((self.data - mean(self.data)).^2)/(length(self.data) - 1);
            h2 = sigma2*ones(size(self.data));
            e = zeros(size(self.data));
            e(1,1)     = self.data(1,1)     - self.mu - self.rho*mean(self.data);
            e(2:end,1) = self.data(2:end,1) - self.mu - self.rho*self.data(1:end-1,1);
            
            for i = 2:length(self.data)      
                h2(i,1) = self.omega + self.alpha*e(i-1,1)^2 ...
                    + self.beta*h2(i-1,1) + self.gamma*e(i-1,1)^2*(e(i-1,1)<0);
            end   
        end       
        
        function Switch(self)   
            [self.mu,    self.mu0]     =  self.Switch_pair(self.mu,    self.mu0);
            [self.rho,   self.rho0]    =  self.Switch_pair(self.rho,   self.rho0);
            [self.omega, self.omega0]  =  self.Switch_pair(self.omega, self.omega0);
            [self.alpha, self.alpha0]  =  self.Switch_pair(self.alpha, self.alpha0);
            [self.beta,  self.beta0]   =  self.Switch_pair(self.beta,  self.beta0);
            [self.gamma, self.gamma0]  =  self.Switch_pair(self.gamma, self.gamma0);
        end
    end
end