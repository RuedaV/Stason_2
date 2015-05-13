%
% Абстрактный класс всех GJR поднобных моделей нужен
% так как estimate и getCondition одинаковы для всех GJR
% также общие параметры omega alpha beta gamma
% стоит лишь переопределять
% фукнции - condVar, 
% getX0, Install, param0Size, LogL()
classdef (Abstract) GJR_BaseModel < IBasicModel
    properties
        omega,  alpha,  beta,  gamma, 
        omega0, alpha0, beta0, gamma0
    end
    
    % Абстрактные виртуальные функции
    % Которые нужно переопределить
    methods (Access = protected)
        function [Aineq, bineq] = GetConditions (self) 
            s = 1e-10;
            A1 = [zeros(self.Param0size),     -1,      0,     0,      0]; % omega > 0
            A2 = [zeros(self.Param0size),      0,     -1,     0,      0]; % alpha >= 0
            A3 = [zeros(self.Param0size),      0,      0,    -1,      0]; % beta >= 0
            A4 = [zeros(self.Param0size),      0,     -1,     0,     -1]; % alpha + gamma >= 0
            
            Aineq = [A1; A2; A3; A4];
            bineq = [-s;  -s;  -s;  -s];
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
        Install(self, x)
        LogL(self)
        CondVar(self)
        GetX0 (self)     
        Switch(self) % просто меняет местами param0 и param
    end
end