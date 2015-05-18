clear;
close all
filename = {'HandSeng.csv'};
% filename = {'CAC40.csv', 'FTSE100.csv', 'HandSeng.csv', ... 
%     'Nikkei225.csv', 'S&P500.csv'};
filenameSize = 1;

% �������� ����� �������� �������
delta0 = 5.001;
mu0 = 0.0006;
rho0 = -0.006;
omega0 = 0.000002;
alpha0 = 0.03;
beta0 = 0.94;
gamma0 = 0.04;
window = 2000;

% ������� ������ �������
% ���� ����������
models = {
%     GJR_Const(0, 0.0001, 0.006, 0.9, 0.11), ...
%     GJR_AR(0, 0.0001, 0.006, 0.9, 0.11), ...
%     GJR_Var(1, 0.0001, 0.006, 0.9, 0.11), ...
%     GJR_Ssn(0, 0, 0, 0, 0, 0.0001, 0.006, 0.9, 0.11),...
%   EGARCH_Const(0, -0.2, 0.1, 0.9, -0.08), ...
  EGARCH_AR(0, -0.2, 0.1, 0.9, -0.08)};
%   EGARCH_Var(0, 0, 0.1, 0.9, -0.08), ...
%   EGARCH_Ssn(0, 0, 0, 0, 0, -0.2, 0.1, 0.9, -0.08)};
modelsNumber = 1;
      
% ������� ���� �� ���� ������
for fileIndex = 1:filenameSize
    %��������� ������
    datasource = xlsread(filename{fileIndex});
    n = length(datasource(:,1));
    d = datasource(1:n, 1);
    r = datasource(1:n, 2);
    disp (filename{fileIndex})
    window  = length(r);
    
    % ������ �� ���� ���������
    for modelIndex = 1:8
        currentModel = models{modelIndex};
        
        % ������ ��� �������� ������������� ������
        clear RESULT_COLIBRATION_VECTOR;
        
        length(r)
        i = 1;
        while i <= length(r)-window+1  
           
            returns = r(i:i+window-1, 1);
            currentModel.data = returns;
            currentModel.setDaily(d(i:i+window-1, 1));
            % ��� ��� ��������� � try
            % �������������� ������������ ������ � 
            % ���� ���������� ������ �� ��������� �������� � 
            % catch � ��� ������������ ������ � ��������� ������ ������
            % ��������
            try
                [est_x, logL] = currentModel.Estimate();
                % �� ������ ���� ���������� ����� ���������
                % �������������� ������
                if (i == 1)
                    RESULT_COLIBRATION_VECTOR = est_x;
                end 
                
                % ������� ���������, ������ ���� ���������� �� ���� �������
                % ���������� ������� � ���������� ������� ���� ������ �����
                if (currentModel.Compare(RESULT_COLIBRATION_VECTOR, est_x) == 1) 
                    RESULT_COLIBRATION_VECTOR = est_x;
                end
            catch
                disp ('������ ������ ������� estimate. �� �� ��������� ������');
                currentModel.errors(fileIndex) = currentModel.errors(fileIndex) + 1;
            end
            
            disp(RESULT_COLIBRATION_VECTOR);
            i = i + 100;
        end
        disp ('Model done');
        % ��������� � ������� ������ ����������
        if (exist('RESULT_COLIBRATION_VECTOR','var') == 1)
            currentModel.collibration{fileIndex} = RESULT_COLIBRATION_VECTOR;
        end
    end
end

%% ��� �� ������ �������� ����������?
% models{1}.collibration{1} % ��������� ��� ������� ������� ������ ������
% models{1}.collibration{2} % ��������� ��� ������� ������� ������ ������
% models{1}.collibration{3}
% models{1}.collibration{4}
% models{1}.collibration{5}

