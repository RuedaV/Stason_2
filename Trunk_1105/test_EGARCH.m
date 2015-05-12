%%
clear all
close all

MUMUMU = 0.001;
RHORHORHO = 0.15;
DELTADELTADELTA = 3.303;
OMEGA = -0.8;
ALPHA = 0.25;
BETTA = 0.93;
GAMMA = -0.12;
SIMULATION_NUMBER = 200;
NBINS = 200;
%% EGARCH_Const_AR
fprintf('\n%12s\r', 'EGARCH_Const_AR');

mu0 = MUMUMU;
rho0 = RHORHORHO;
omega0 = OMEGA;
alpha0 = ALPHA;
beta0 = BETTA;
gamma0 = GAMMA;

S = SIMULATION_NUMBER;
T = 1000;

%x0 = [mu0, rho0, omega0, alpha0, beta0, gamma0];
omega  = zeros(S, 1);
alpha  = zeros(S, 1);
beta   = zeros(S, 1);
gamma  = zeros(S, 1);
loss   = zeros(S, 1);
loss2  = zeros(S, 1);
VaR_a  = 0;
VaR_b  = 0;
for i = 1:S
    i
    try
    a  = EGARCH_Const_AR(mu0, rho0, omega0, alpha0, beta0, gamma0);
    b  = EGARCH_Const(mu0, omega0, alpha0, beta0, gamma0);
    Data = a.Simulate(T);    
    a.data = Data(1:end-1);
    a.data_plus = Data;    
    b.data = Data(1:end-1);
    b.data_plus = Data;
    % True
    a.Estimate();
    [loss_a, loss2_a, VaR_exceeded_a] = a.Predict();
    VaR_a = VaR_a + VaR_exceeded_a;
%     a.MyDisplay();
    % False    
    b.Estimate();
    [loss_b, loss2_b, VaR_exceeded_b] = b.Predict();
    VaR_b = VaR_b + VaR_exceeded_b;
%     b.MyDisplay();
    omega(i,1) = (a.omega - b.omega)/omega0;
    alpha(i,1) = (a.alpha - b.alpha)/alpha0;
    beta(i,1)  = (a.beta  - b.beta)/beta0;
    gamma(i,1) = (a.gamma - b.gamma)/gamma0;
    loss(i,1)  = (loss_a - loss_b)/loss_a;
    loss2(i,1)  = (loss2_a - loss2_b)/loss2_a;
    i
    catch
        disp ('Ошибка');
    end
end
nbins = NBINS;
MyHistEl( omega, 0, nbins, '$\bf \omega$')
MyHistEl( alpha, 0, nbins, '$\alpha$')
MyHistEl( beta,  0, nbins, '$\beta$')
MyHistEl( gamma, 0, nbins, '$\gamma$')
MyHistEl( loss,  0, nbins, '$L$')
% fprintf('\n%6s %12s %12s\r', 'Model',    'estimated','realized');
% fprintf('%6s %12.6f %12.6f\n', 'True',    h2_pred_a,   h2_proxy_a);
% fprintf('%6s %12.6f %12.6f\n', 'False',   h2_pred_b,   h2_proxy_b);

fprintf('\n%6s %12s \r', 'Model', 'VaR is exceeded');
fprintf('%6s %12.6f \n', 'True',    VaR_a);
fprintf('%6s %12.6f \n', 'False',   VaR_b);

DisplayStats(omega, alpha, beta, gamma, loss, loss2);




%% EGARCH_Const_Var
fprintf('\n%12s\r', 'EGARCH_Const_Var');

mu0 = MUMUMU;
delta0 = DELTADELTADELTA;
omega0 = OMEGA
alpha0 = ALPHA;
beta0 = BETTA;
gamma0 = GAMMA;

S = 2000;
T = 1000;

%x0 = [mu0, delta0, omega0, alpha0, beta0, gamma0];
omega  = zeros(S, 1);
alpha  = zeros(S, 1);
beta   = zeros(S, 1);
gamma  = zeros(S, 1);
loss   = zeros(S, 1);
loss2  = zeros(S, 1);
VaR_a  = 0;
VaR_b  = 0;
for i = 1:S
    i
    try 
    a  = EGARCH_Const_Var(mu0, delta0, omega0, alpha0, beta0, gamma0);
    b  = EGARCH_Const(mu0, omega0, alpha0, beta0, gamma0);
    Data = a.Simulate(T);    
    a.data = Data(1:end-1);
    a.data_plus = Data;    
    b.data = Data(1:end-1);
    b.data_plus = Data;
        % True
    a.Estimate();
    [loss_a, loss2_a, VaR_exceeded_a] = a.Predict();
    VaR_a = VaR_a + VaR_exceeded_a;
%     a.MyDisplay();
    % False    
    b.Estimate();
    [loss_b, loss2_b, VaR_exceeded_b] = b.Predict();
    VaR_b = VaR_b + VaR_exceeded_b;
%     b.MyDisplay();
    omega(i,1) = (a.omega - b.omega)/omega0;
    alpha(i,1) = (a.alpha - b.alpha)/alpha0;
    beta(i,1)  = (a.beta  - b.beta)/beta0;
    gamma(i,1) = (a.gamma - b.gamma)/gamma0;
    loss(i,1)  = (loss_a - loss_b)/loss_a;
    loss2(i,1)  = (loss2_a - loss2_b)/loss2_a;
    catch
        disp ('Ошибка');
    end
end
nbins = NBINS;
MyHistEl( omega, 0, nbins, '$\bf \omega$')
MyHistEl( alpha, 0, nbins, '$\alpha$')
MyHistEl( beta,  0, nbins, '$\beta$')
MyHistEl( gamma, 0, nbins, '$\gamma$')
MyHistEl( loss,  0, nbins, '$L$')
% fprintf('\n%6s %12s %12s\r', 'Model',    'estimated','realized');
% fprintf('%6s %12.6f %12.6f\n', 'True',    h2_pred_a,   h2_proxy_a);
% fprintf('%6s %12.6f %12.6f\n', 'False',   h2_pred_b,   h2_proxy_b);

fprintf('\n%6s %12s \r', 'Model', 'VaR is exceeded');
fprintf('%6s %12.6f \n', 'True',    VaR_a);
fprintf('%6s %12.6f \n', 'False',   VaR_b);

DisplayStats(omega, alpha, beta, gamma, loss, loss2);



%% EGARCH_Ssn

clear all
close all

MUMUMU = 0.001;
RHORHORHO = 0.15;
DELTADELTADELTA = 3.303;
THETHA01 = -0.000583391;
THETHA02 = 0.0014;
THETHA03 = 0.0022;
THETHA04 = 0.0002;
THETHA05 = 0.0016;
OMEGA = -0.8;
ALPHA = 0.25;
BETTA = 0.93;
GAMMA = -0.12;
SIMULATION_NUMBER = 2000;
NBINS = 200;

fprintf('\n%12s\r', 'EGARCH_Ssn');

mu0 = MUMUMU;
theta10 = THETHA01;
theta20 = THETHA02;
theta30 = THETHA03;
theta40 = THETHA04;
theta50 = THETHA05;
omega0 = OMEGA;
alpha0 = ALPHA;
beta0 = BETTA;
gamma0 = GAMMA;


S = SIMULATION_NUMBER;
T = 1000;

%x0 = [theta10, theta20, theta30, theta40, theta50, omega0, alpha0, beta0, gamma0];
omega  = zeros(S, 1);
alpha  = zeros(S, 1);
beta   = zeros(S, 1);
gamma  = zeros(S, 1);
loss   = zeros(S, 1);
loss2  = zeros(S, 1);
VaR_a  = 0;
VaR_b  = 0;
for i = 1:S
    i
%     try
    a  = EGARCH_Ssn(theta10, theta20, theta30, theta40, theta50,...
        omega0, alpha0, beta0, gamma0);
    b  = EGARCH_Const(mu0, omega0, alpha0, beta0, gamma0);
    [Data, Day] = a.Simulate(T);    
    a.data = Data(1:end-1);
    a.data_plus = Data;  
    
    a.day = Day(1:end-1);
    a.day_plus = Day;
    
    b.data = Data(1:end-1);
    b.data_plus = Data;
    % True
    a.Estimate();
    [loss_a, loss2_a, VaR_exceeded_a] = a.Predict();
    VaR_a = VaR_a + VaR_exceeded_a;
%     a.MyDisplay();
    % False    
    b.Estimate();
    [loss_b, loss2_b, VaR_exceeded_b] = b.Predict();
    VaR_b = VaR_b + VaR_exceeded_b;
%     b.MyDisplay();
    omega(i,1) = (a.omega - b.omega)/omega0;
    alpha(i,1) = (a.alpha - b.alpha)/alpha0;
    beta(i,1)  = (a.beta  - b.beta)/beta0;
    gamma(i,1) = (a.gamma - b.gamma)/gamma0;
    loss(i,1)  = (loss_a - loss_b)/loss_a;
    loss2(i,1) = (loss2_a - loss2_b)/loss2_a;
    
%     fprintf('%6s %12s %12s\r','param','true', 'estimated');
%     fprintf('%6s %12.6f %12.6f\n', 'theta1',    a.theta10,    a.theta1);
%     fprintf('%6s %12.6f %12.6f\n', 'theta2',    a.theta20,    a.theta2);
%     fprintf('%6s %12.6f %12.6f\n', 'theta3',    a.theta30,    a.theta3);
%     fprintf('%6s %12.6f %12.6f\n', 'theta4',    a.theta40,    a.theta4);
%     fprintf('%6s %12.6f %12.6f\n', 'theta5',    a.theta50,    a.theta5);
%     fprintf('%6s %12.6f %12.6f\n', 'omega',     a.omega0,     a.omega);
%     fprintf('%6s %12.6f %12.6f\n', 'alpha',     a.alpha0,     a.alpha);
%     fprintf('%6s %12.6f %12.6f\n', 'beta',      a.beta0,      a.beta);
%     fprintf('%6s %12.6f %12.6f\n', 'gamma',     a.gamma0,     a.gamma);
%     catch  
%         disp ('Ошибка');
%     end
end
nbins = 150;
MyHistEl( omega, 0, nbins, '$\bf \omega$')
MyHistEl( alpha, 0, nbins, '$\alpha$')
MyHistEl( beta,  0, nbins, '$\beta$')
MyHistEl( gamma, 0, nbins, '$\gamma$')
MyHistEl( loss,  0, nbins, '$L$')
% fprintf('\n%6s %12s %12s\r', 'Model',    'estimated','realized');
% fprintf('%6s %12.6f %12.6f\n', 'True',    h2_pred_a,   h2_proxy_a);
% fprintf('%6s %12.6f %12.6f\n', 'False',   h2_pred_b,   h2_proxy_b);

fprintf('\n%6s %12s \r', 'Model', 'VaR is exceeded');
fprintf('%6s %12.6f \n', 'True',    VaR_a);
fprintf('%6s %12.6f \n', 'False',   VaR_b);

DisplayStats(omega, alpha, beta, gamma, loss, loss2);



%% EGARCH_General
fprintf('\n%12s\r', 'EGARCH_General');

clear all
close all

MUMUMU = 0.001;
RHORHORHO = 0.15;
DELTADELTADELTA = 3.303;
THETHA01 = -0.000583391;
THETHA02 = 0.0014;
THETHA03 = 0.0022;
THETHA04 = 0.0002;
THETHA05 = 0.0016;
OMEGA = -0.8;
ALPHA = 0.25;
BETTA = 0.93;
GAMMA = -0.12;
SIMULATION_NUMBER = 2000;
NBINS = 200;

mu0 = MUMUMU;
rho0 = RHORHORHO;
delta0 = DELTADELTADELTA;
theta10 = THETHA01;
theta20 = THETHA02;
theta30 = THETHA03;
theta40 = THETHA04;
theta50 = THETHA05;
omega0 = OMEGA;
alpha0 = ALPHA;
beta0 = BETTA;
gamma0 = GAMMA;

S = SIMULATION_NUMBER;
T = 1000;

x0 = [theta10, theta20, theta30, theta40, theta50, omega0, alpha0, beta0, gamma0];
omega  = zeros(S, 1);
alpha  = zeros(S, 1);
beta   = zeros(S, 1);
gamma  = zeros(S, 1);
loss   = zeros(S, 1);
loss2  = zeros(S, 1);
VaR_a  = 0;
VaR_b  = 0;
for i = 1:S
    i
%     try
    a  = EGARCH_General(rho0, delta0, theta10, theta20, theta30, theta40, theta50,...
        omega0, alpha0, beta0, gamma0);
    b  = EGARCH_Const(mu0, omega0, alpha0, beta0, gamma0);
    [Data, Day] = a.Simulate(T);    
    a.data = Data(1:end-1);
    a.data_plus = Data;  
    a.day = Day(1:end-1);
    a.day_plus = Day;
    b.data = Data(1:end-1);
    b.data_plus = Data;
    % True
    a.Estimate();
    [loss_a, loss2_a, VaR_exceeded_a] = a.Predict();
    VaR_a = VaR_a + VaR_exceeded_a;
%     a.MyDisplay();
    % False    
    b.Estimate();
    [loss_b, loss2_b, VaR_exceeded_b] = b.Predict();
    VaR_b = VaR_b + VaR_exceeded_b;
%     b.MyDisplay();
    omega(i,1) = (a.omega - b.omega)/omega0;
    alpha(i,1) = (a.alpha - b.alpha)/alpha0;
    beta(i,1)  = (a.beta  - b.beta)/beta0;
    gamma(i,1) = (a.gamma - b.gamma)/gamma0;
    loss(i,1)  = (loss_a - loss_b)/loss_a;
    loss2(i,1)  = (loss2_a - loss2_b)/loss2_a;
%     catch
%         disp ('Ошибка');
%     end
end
nbins = NBINS;
MyHistEl( omega, 0, nbins, '$\bf \omega$')
MyHistEl( alpha, 0, nbins, '$\alpha$')
MyHistEl( beta,  0, nbins, '$\beta$')
MyHistEl( gamma, 0, nbins, '$\gamma$')
MyHistEl( loss,  0, nbins, '$L$')
% fprintf('\n%6s %12s %12s\r', 'Model',    'estimated','realized');
% fprintf('%6s %12.6f %12.6f\n', 'True',    h2_pred_a,   h2_proxy_a);
% fprintf('%6s %12.6f %12.6f\n', 'False',   h2_pred_b,   h2_proxy_b);

fprintf('\n%6s %12s \r', 'Model', 'VaR is exceeded');
fprintf('%6s %12.6f \n', 'True',    VaR_a);
fprintf('%6s %12.6f \n', 'False',   VaR_b);

DisplayStats(omega, alpha, beta, gamma, loss, loss2);


