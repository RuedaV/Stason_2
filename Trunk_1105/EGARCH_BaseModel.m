classdef (Abstract) EGARCH_BaseModel < IBasicModel
    properties        
        omega,  alpha,  beta,  gamma, 
        omega0, alpha0, beta0, gamma0
    end
    
    methods (Access = protected)
        function [Aineq, bineq] = GetConditions (self) 
            % beta < 1000
            A1 = [zeros(self.Param0size),     0,      0,     1,      0]; 
            Aineq = A1;
            bineq = 1000;
        end
    end 
    
    methods (Access = public)
        function res = Compare(self, x, rival)
            if abs(x(1)) < abs(rival(1))
                res = 1;
            else
                res = 0;
            end                
        end
        
        function setDaily (self, d)                           
        end
        
        function [loss, loss2, VaR_exceeded] = Predict(self)
            [h2, e] = self.CondVar();
            log_h2_pred = self.omega + ...
                        self.alpha*( abs(e(end,1))/sqrt(h2(end,1)) ) + ...
                        self.beta*log(h2(end,1)) + ...
                        self.gamma*e(end,1)/sqrt(h2(end,1));
                    
            h2_pred = exp(log_h2_pred);
            
            VaR = (self.data(end, 1) - e(end, 1)) + sqrt(h2(end,1))*norminv(0.05,0,1);
            VaR_exceeded = (VaR > self.data_plus(end,1));
            
            data_temp = self.data;            
            self.data = self.data_plus;
            self.Switch();
            [h2_plus, e_plus] = self.CondVar();
            h2_proxy = h2_plus(end, 1);
            self.Switch();
            self.data = data_temp;
            loss  = QLIKE(h2_proxy, h2_pred);
            loss2 = QLIKE2(h2_proxy, h2_pred);
        end
    end
    
    methods (Abstract, Access = protected) 
        Param0size (self)
        Install (self, x)
        LogL(self)
        CondVar(self)
        GetX0 (self)
        Switch(self) % просто меняет местами param0 и param
    end
end