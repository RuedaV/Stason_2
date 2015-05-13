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