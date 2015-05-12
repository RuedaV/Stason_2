classdef (Abstract) IBasicModel < handle
    % �������� (���������� �������)
    properties 
        data, % ������ ��� ����������
        data_plus, % ������ ����� ��� ���������� ������� ������, ����� ������� ������
        % ������ ��� ������������
        % ������ �������� ��� ������ ��������
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
            % 5 �������� ��� 5�� ��������
            % � ��� 5 �� ��������
        end
        
        function [a_new, b_new] = Switch_pair(self, a, b)
            a_new = b;
            b_new = a;
        end
    end 
     
    % ������� �����������, ������� ���� �������������� � ������ ����������
    % ������� ������
    methods (Abstract, Access = public) 
        Compare(self, x, rival) % ���������� ������������ �������� 

    end
    % �������� ����������� �������
    methods (Abstract, Access = protected) 
        GetConditions (self) % ��� ������ estimate ����� ��������� �������
        Param0size (self)    % ���������� ������ ������� ���������
        Install (self, x)    % ������������ ��������� ������ � ������� x
        LogL(self)           % ������� ������� ������������� �������������
        CondVar(self)        % ������� �������� �����������
        GetX0 (self)         % ���������� �������������� ������ ��� fmincon
    end
end
