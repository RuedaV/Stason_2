classdef (Abstract) IBasicModel < handle
    % Свойства (переменные обьекта)
    properties 
        data, % данные для оценивания
        data_plus, % данные чисто для построения функции потерь, чтобы считать ошибки
        % Данные для колибрования
        % массив векторов для разных индексов
        collibration,
        % error
        errors
    end
    methods (Access = public)
        function [est_x, LogL] = Estimate (self)
            x0 = self.GetX0();
            [Aineq, bineq] = self.GetConditions();
            [est_x, LogL] = fmincon(@(x)(-self.LogLForEstimate(self.data, x)), ...
                    x0,Aineq,bineq,[],[],[],[],[], ...
                    optimset('Algorithm','interior-point','Hessian','bfgs', 'Display', 'off'));
            LogL = - LogL;
        end
        
        function self = IBasicModel()
            self.setCollibration();
        end
        
        function showCollibration (self)
            fprintf('%12s\n', 'MODEL 1')
            fprintf('%12s %12s %12s %12s %12s %12s\n', 'model',    'mu',    'omega',     'alpha',     'beta',    'gamma');
            fprintf('%12s %12.6f %12.6f %12.6f %12.6f %12.6f\n', 'CAC40',     self.collibration{1}(1),    self.collibration{1}(2),    self.collibration{1}(3), self.collibration{1}(4), self.collibration{1}(5));
            fprintf('%12s %12.6f %12.6f %12.6f %12.6f %12.6f\n', 'FTSE100',   self.collibration{2}(1),    self.collibration{2}(2),    self.collibration{2}(3), self.collibration{2}(4), self.collibration{2}(5));
            fprintf('%12s %12.6f %12.6f %12.6f %12.6f %12.6f\n', 'HandSeng',  self.collibration{3}(1),    self.collibration{3}(2),    self.collibration{3}(3), self.collibration{3}(4), self.collibration{3}(5));
            fprintf('%12s %12.6f %12.6f %12.6f %12.6f %12.6f\n', 'Nikkei225', self.collibration{4}(1),    self.collibration{4}(2),    self.collibration{4}(3), self.collibration{4}(4), self.collibration{4}(5));
            fprintf('%12s %12.6f %12.6f %12.6f %12.6f %12.6f\n', 'S&P500',    self.collibration{5}(1),    self.collibration{5}(2),    self.collibration{5}(3), self.collibration{5}(4), self.collibration{5}(5));
        end    
    end

    
    methods (Access = protected) 
        function log = LogLForEstimate (self, data, x) 
            self.Install (x, data);
            log = self.LogL();
        end 
        
       
        function setCollibration (self) 
            self.collibration = {[];[];[];[];[]};
            self.errors = [0, 0, 0, 0, 0];
            % 5 векторов для 5ти индексов
            % и еще 5 на прозапас
        end
        
        function [a_new, b_new] = Switch_pair(self, a, b)
            a_new = b;
            b_new = a;
        end
    end 
     
    % функции виртуальные, которые надо переопределять в каждой реализации
    % данного класса
    methods (Abstract, Access = public) 
        Compare(self, x, rival) % Сравнивает колибрируюий параметр 

    end
    % закрытые виртуальные функции
    methods (Abstract, Access = protected) 
        GetConditions (self) % Для любого estimate нужны граничные условия
        Param0size (self)    % Возвращает размер первого параметра
        Install (self, x)    % Приравнивает параметры модели к вектору x
        LogL(self)           % Считает функцию максимального правдоподобия
        CondVar(self)        % Считает условную вероятность
        GetX0 (self)         % Возвращает первоначальный вектор для fmincon
    end
end
