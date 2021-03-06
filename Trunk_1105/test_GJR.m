%%
clear all
close all
MUMUMU = 0.001;
RHORHORHO = 0.009;
DELTADELTADELTA = 3.303;
THETHA01 = -0.000583391;
THETHA02 = 0.00143266434255308;
THETHA03 = 0.00224614578474341;
THETHA04 = 0.000154878290132809;
THETHA05 = 0.0015718522946019;
OMEGA = 0.000004;
ALPHA = 0.033;
BETTA = 0.9;
GAMMA = 0.08;
SIMULATION_NUMBER = 10;
NBINS = 50;
p = 0.05;

%% GJR_Const_AR
fprintf('\n%12s\r', 'GJR_Const_AR');

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
VaR_a  = 0;
VaR_b  = 0;
err = 0;

i = 1;
while i <= S
    i
    try
    a  = GJR_Const_AR(mu0, rho0, omega0, alpha0, beta0, gamma0);
    b  = GJR_Const(mu0, omega0, alpha0, beta0, gamma0);
    Data = a.Simulate(T);    
    a.data = Data(1:end-1);
    a.data_plus = Data;    
    b.data = Data(1:end-1);
    b.data_plus = Data;
    b.sigma2 = a.sigma2;
    % True
    a.Estimate();
    [loss_a, loss2_a, VaR_exceeded_a] = a.Predict(p);
    VaR_a = VaR_a + VaR_exceeded_a;
%     a.MyDisplay();
    % False    
    b.Estimate();
    [loss_b, loss2_b, VaR_exceeded_b] = b.Predict(p);
    VaR_b = VaR_b + VaR_exceeded_b;
%     b.MyDisplay();
    omega(i,1) = (a.omega - b.omega)/omega0;
    alpha(i,1) = (a.alpha - b.alpha)/alpha0;
    beta(i,1)  = (a.beta  - b.beta)/beta0;
    gamma(i,1) = (a.gamma - b.gamma)/gamma0;
    loss(i,1)  = (loss_a - loss_b)/QLIKE(a.sigma2(end,1),a.sigma2(end,1));
    if and(isnan(loss(i,1))== 0,(abs(omega(i,1)) < 10))
        i = i + 1;
        err = 0;
    end
    catch
        err = err + 1;
        disp ('������');  
        if err > 10
            i = i + 1;
            err = 0;
        end
    end
end

DisplayStats(omega, alpha, beta, gamma, loss);

nbins = NBINS;
MyHistEl( omega, 0, nbins, '$\frac{\hat\omega - \hat\omega''}{\omega}$')
MyHistEl( alpha, 0, nbins, '$\frac{\hat\alpha - \hat\alpha''}{\alpha}$')
MyHistEl( beta,  0, nbins, '$\frac{\hat\beta - \hat\beta''}{\beta}$')
MyHistEl( gamma, 0, nbins, '$\frac{\hat\gamma - \hat\gamma''}{\gamma}$')
MyHistEl( loss,  0, nbins, '$\frac{\hat L - \hat L''}{L}$')
% fprintf('\n%6s %12s %12s\r', 'Model',    'estimated','realized');
% fprintf('%6s %12.6f %12.6f\n', 'True',    h2_pred_a,   h2_proxy_a);
% fprintf('%6s %12.6f %12.6f\n', 'False',   h2_pred_b,   h2_proxy_b);

p_value_a = Xi_squared(p, SIMULATION_NUMBER,  VaR_exceeded_a);
p_value_b = Xi_squared(p, SIMULATION_NUMBER,  VaR_exceeded_b);

fprintf('\n%6s %12s \r', 'Model', 'p-value');
fprintf('%6s %12.3f \n', 'True',    p_value_a);
fprintf('%6s %12.3f \n', 'False',   p_value_a);




%% GJR_Const_Var
fprintf('\n%12s\r', 'GJR_Const_Var');

mu0 = MUMUMU;
delta0 = DELTADELTADELTA;
omega0 = OMEGA;
alpha0 = ALPHA;
beta0 = BETTA;
gamma0 = GAMMA;

S = SIMULATION_NUMBER;
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

err = 0;

i = 1;
while i <= S
    i
    try 
    a  = GJR_Const_Var(mu0, delta0, omega0, alpha0, beta0, gamma0);
    b  = GJR_Const(mu0, omega0, alpha0, beta0, gamma0);
    Data = a.Simulate(T);    
    a.data = Data(1:end-1);
    a.data_plus = Data;    
    b.data = Data(1:end-1);
    b.data_plus = Data;
    b.sigma2 = a.sigma2;
        % True
    a.Estimate();
    [loss_a, loss2_a, VaR_exceeded_a] = a.Predict(p);
    VaR_a = VaR_a + VaR_exceeded_a;
%     a.MyDisplay();
    % False    
    b.Estimate();
    [loss_b, loss2_b, VaR_exceeded_b] = b.Predict(p);
    VaR_b = VaR_b + VaR_exceeded_b;
%     b.MyDisplay();
    omega(i,1) = (a.omega - b.omega)/omega0;
    alpha(i,1) = (a.alpha - b.alpha)/alpha0;
    beta(i,1)  = (a.beta  - b.beta)/beta0;
    gamma(i,1) = (a.gamma - b.gamma)/gamma0;
    loss(i,1)  = (loss_a - loss_b)/QLIKE(a.sigma2(end,1),a.sigma2(end,1));
    if and(isnan(loss(i,1))== 0,(abs(omega(i,1)) < 10))
        i = i + 1;
        err = 0;
    end
    catch
        err = err + 1;
        disp ('������');  
        if err > 10
            i = i + 1;
            err = 0;
        end
    end
end

DisplayStats(omega, alpha, beta, gamma, loss);

nbins = NBINS;
MyHistEl( omega, 0, nbins, '$\frac{\hat\omega - \hat\omega''}{\omega}$')
MyHistEl( alpha, 0, nbins, '$\frac{\hat\alpha - \hat\alpha''}{\alpha}$')
MyHistEl( beta,  0, nbins, '$\frac{\hat\beta - \hat\beta''}{\beta}$')
MyHistEl( gamma, 0, nbins, '$\frac{\hat\gamma - \hat\gamma''}{\gamma}$')
MyHistEl( loss,  0, nbins, '$\frac{\hat L - \hat L''}{L}$')
% fprintf('\n%6s %12s %12s\r', 'Model',    'estimated','realized');
% fprintf('%6s %12.6f %12.6f\n', 'True',    h2_pred_a,   h2_proxy_a);
% fprintf('%6s %12.6f %12.6f\n', 'False',   h2_pred_b,   h2_proxy_b);

fprintf('\n%6s %12s \r', 'Model', 'VaR is exceeded');
fprintf('%6s %12.6f \n', 'True',    VaR_a);
fprintf('%6s %12.6f \n', 'False',   VaR_b);



%% GJR_Ssn
fprintf('\n%12s\r', 'GJR_Ssn');

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

err = 0;

i = 1;
while i <= S
    i
    try
    a  = GJR_Ssn(theta10, theta20, theta30, theta40, theta50,...
        omega0, alpha0, beta0, gamma0);
    b  = GJR_Const(mu0, omega0, alpha0, beta0, gamma0);
    [Data, Day] = a.Simulate(T);    
    a.data = Data(1:end-1);
    a.data_plus = Data;  
    a.day = Day(1:end-1);
    a.day_plus = Day;
    b.data = Data(1:end-1);
    b.data_plus = Data;
    b.sigma2 = a.sigma2;
        % True
    a.Estimate();
    [loss_a, loss2_a, VaR_exceeded_a] = a.Predict(p);
    VaR_a = VaR_a + VaR_exceeded_a;
%     a.MyDisplay();
    % False    
    b.Estimate();
    [loss_b, loss2_b, VaR_exceeded_b] = b.Predict(p);
    VaR_b = VaR_b + VaR_exceeded_b;
%     b.MyDisplay();
    omega(i,1) = (a.omega - b.omega)/omega0;
    alpha(i,1) = (a.alpha - b.alpha)/alpha0;
    beta(i,1)  = (a.beta  - b.beta)/beta0;
    gamma(i,1) = (a.gamma - b.gamma)/gamma0;
    loss(i,1)  = (loss_a - loss_b)/loss_a;
    loss2(i,1)  = (loss2_a - loss2_b)/loss2_a;
    if and(isnan(loss(i,1))== 0,(abs(omega(i,1)) < 10))
        i = i + 1;
        err = 0;
    end
    catch
        err = err + 1;
        disp ('������');  
        if err > 10
            i = i + 1;
            err = 0;
        end
    end
end
nbins = NBINS;
MyHistEl( omega, 0, nbins, '$\frac{\hat\omega - \hat\omega''}{\omega}$')
MyHistEl( alpha, 0, nbins, '$\frac{\hat\alpha - \hat\alpha''}{\alpha}$')
MyHistEl( beta,  0, nbins, '$\frac{\hat\beta - \hat\beta''}{\beta}$')
MyHistEl( gamma, 0, nbins, '$\frac{\hat\gamma - \hat\gamma''}{\gamma}$')
MyHistEl( loss,  0, nbins, '$\frac{\hat L - \hat L''}{L}$')
MyHistEl( loss2,  0, nbins, '$\frac{\hat L_2 - \hat L_2''}{L_2}$')
% fprintf('\n%6s %12s %12s\r', 'Model',    'estimated','realized');
% fprintf('%6s %12.6f %12.6f\n', 'True',    h2_pred_a,   h2_proxy_a);
% fprintf('%6s %12.6f %12.6f\n', 'False',   h2_pred_b,   h2_proxy_b);

fprintf('\n%6s %12s \r', 'Model', 'VaR is exceeded');
fprintf('%6s %12.6f \n', 'True',    VaR_a);
fprintf('%6s %12.6f \n', 'False',   VaR_b);

DisplayStats(omega, alpha, beta, gamma, loss);

%% GJR_General
fprintf('\n%12s\r', 'GJR_General');

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

err = 0;

i = 1;

while i <= S
    i
    try
    a  = GJR_General(rho0, delta0, theta10, theta20, theta30, theta40, theta50,...
        omega0, alpha0, beta0, gamma0);
    b  = GJR_Const(mu0, omega0, alpha0, beta0, gamma0);
    [Data, Day] = a.Simulate(T);    
    a.data = Data(1:end-1);
    a.data_plus = Data;  
    a.day = Day(1:end-1);
    a.day_plus = Day;
    b.data = Data(1:end-1);
    b.data_plus = Data;
    b.sigma2 = a.sigma2;
        % True
    a.Estimate();
    [loss_a, loss2_a, VaR_exceeded_a] = a.Predict(p);
    VaR_a = VaR_a + VaR_exceeded_a;
%     a.MyDisplay();
    % False    
    b.Estimate();
    [loss_b, loss2_b, VaR_exceeded_b] = b.Predict(p);
    VaR_b = VaR_b + VaR_exceeded_b;
%     b.MyDisplay();
    omega(i,1) = (a.omega - b.omega)/omega0;
    alpha(i,1) = (a.alpha - b.alpha)/alpha0;
    beta(i,1)  = (a.beta  - b.beta)/beta0;
    gamma(i,1) = (a.gamma - b.gamma)/gamma0;
    loss(i,1)  = (loss_a - loss_b)/QLIKE(a.sigma2(end,1),a.sigma2(end,1));
    if and(isnan(loss(i,1))== 0,(abs(omega(i,1)) < 10))
        i = i + 1;
        err = 0;
    end
    catch
        err = err + 1;
        disp ('������');  
        if err > 10
            i = i + 1;
            err = 0;
        end
    end
end

DisplayStats(omega, alpha, beta, gamma, loss);


nbins = NBINS;
MyHistEl( omega, 0, nbins, '$\frac{\hat\omega - \hat\omega''}{\omega}$')
MyHistEl( alpha, 0, nbins, '$\frac{\hat\alpha - \hat\alpha''}{\alpha}$')
MyHistEl( beta,  0, nbins, '$\frac{\hat\beta - \hat\beta''}{\beta}$')
MyHistEl( gamma, 0, nbins, '$\frac{\hat\gamma - \hat\gamma''}{\gamma}$')
MyHistEl( loss,  0, nbins, '$\frac{\hat L - \hat L''}{L}$')
% fprintf('\n%6s %12s %12s\r', 'Model',    'estimated','realized');
% fprintf('%6s %12.6f %12.6f\n', 'True',    h2_pred_a,   h2_proxy_a);
% fprintf('%6s %12.6f %12.6f\n', 'False',   h2_pred_b,   h2_proxy_b);

fprintf('\n%6s %12s \r', 'Model', 'VaR is exceeded');
fprintf('%6s %12.6f \n', 'True',    VaR_a);
fprintf('%6s %12.6f \n', 'False',   VaR_b);



