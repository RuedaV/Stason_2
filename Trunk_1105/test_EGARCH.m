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
SIMULATION_NUMBER = 100;
NBINS = 20;
p = 0.05;
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
VaR    = zeros(S, 1);
VaR_a  = zeros(S, 1);
VaR_b  = zeros(S, 1);
r_last = zeros(S, 1);

err = 0;
i = 1;
while i <= S
    i
    try
    a  = EGARCH_Const_AR(mu0, rho0, omega0, alpha0, beta0, gamma0);
    b  = EGARCH_Const(mu0, omega0, alpha0, beta0, gamma0);
    Data = a.Simulate(T);    
    a.data = Data(1:end-1);
    a.data_plus = Data;    
    b.data = Data(1:end-1);
    b.data_plus = Data;
    b.sigma2 = a.sigma2;
    % True
    a.Estimate();
    [loss_a, VaR_true_a, VaR_pred_a] = a.Predict(p);

    % False    
    b.Estimate();
    [loss_b, VaR_true_b, VaR_pred_b] = b.Predict(p);
    
    VaR_a(i,1)  = VaR_pred_a;
    VaR_b(i,1)  = VaR_pred_b;
    r_last(i,1) = Data(end,1);

    omega(i,1) = (a.omega - b.omega)/omega0;
    alpha(i,1) = (a.alpha - b.alpha)/alpha0;
    beta(i,1)  = (a.beta  - b.beta)/beta0;
    gamma(i,1) = (a.gamma - b.gamma)/gamma0;
    loss(i,1)  = (loss_a - loss_b)/QLIKE(a.sigma2(end,1),a.sigma2(end,1));
    VaR(i,1)   = (VaR_pred_a - VaR_pred_b)/VaR_true_a;

    if and(isnan(loss(i,1))== 0,(abs(omega(i,1)) < 10))
        i = i + 1;
        err = 0;
    end
    catch
        err = err + 1;
        disp ('Îøèáêà');  
        if err > 10
            i = i + 1;
            err = 0;
        end
    end
end

DisplayStats(omega, alpha, beta, gamma, loss);

nbins = NBINS;
MyHistEl( omega, 0, nbins, '$\bf \omega$')
MyHistEl( alpha, 0, nbins, '$\alpha$')
MyHistEl( beta,  0, nbins, '$\beta$')
MyHistEl( gamma, 0, nbins, '$\gamma$')
MyHistEl( loss,  0, nbins, '$L$')
MyHistEl( VaR,   0, nbins, '$VaR$')


%% EGARCH_Const_Var
fprintf('\n%12s\r', 'EGARCH_Const_Var');

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
VaR    = zeros(S, 1);


err = 0;
i = 1;
while i <= S
    i
    try
    a  = EGARCH_Const_Var(mu0, delta0, omega0, alpha0, beta0, gamma0);
    b  = EGARCH_Const(mu0, omega0, alpha0, beta0, gamma0);
    Data = a.Simulate(T);    
    a.data = Data(1:end-1);
    a.data_plus = Data;    
    b.data = Data(1:end-1);
    b.data_plus = Data;
    b.sigma2 = a.sigma2;
    % True
    a.Estimate();
    [loss_a, VaR_true_a, VaR_pred_a] = a.Predict(p);

    % False    
    b.Estimate();
    [loss_b, VaR_true_b, VaR_pred_b] = b.Predict(p);
    
    omega(i,1) = (a.omega - b.omega)/omega0;
    alpha(i,1) = (a.alpha - b.alpha)/alpha0;
    beta(i,1)  = (a.beta  - b.beta)/beta0;
    gamma(i,1) = (a.gamma - b.gamma)/gamma0;
    loss(i,1)  = (loss_a - loss_b)/QLIKE(a.sigma2(end,1),a.sigma2(end,1));
    VaR(i,1)   = (VaR_pred_a - VaR_pred_b)/VaR_true_a;

    if and(isnan(loss(i,1))== 0,(abs(omega(i,1)) < 10))
        i = i + 1;
        err = 0;
    end
    catch
        err = err + 1;
        disp ('Îøèáêà');  
        if err > 10
            i = i + 1;
            err = 0;
        end
    end
end

DisplayStats(omega, alpha, beta, gamma, loss);

nbins = NBINS;
MyHistEl( omega, 0, nbins, '$\bf \omega$')
MyHistEl( alpha, 0, nbins, '$\alpha$')
MyHistEl( beta,  0, nbins, '$\beta$')
MyHistEl( gamma, 0, nbins, '$\gamma$')
MyHistEl( loss,  0, nbins, '$L$')
MyHistEl( VaR,   0, nbins, '$VaR$')



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
VaR    = zeros(S, 1);


err = 0;
i = 1;
while i <= S
    i
    try
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
    b.sigma2 = a.sigma2;
    % True
    a.Estimate();
    [loss_a, VaR_true_a, VaR_pred_a] = a.Predict(p);

    % False    
    b.Estimate();
    [loss_b, VaR_true_b, VaR_pred_b] = b.Predict(p);
    
    omega(i,1) = (a.omega - b.omega)/omega0;
    alpha(i,1) = (a.alpha - b.alpha)/alpha0;
    beta(i,1)  = (a.beta  - b.beta)/beta0;
    gamma(i,1) = (a.gamma - b.gamma)/gamma0;
    loss(i,1)  = (loss_a - loss_b)/QLIKE(a.sigma2(end,1),a.sigma2(end,1));
    VaR(i,1)   = (VaR_pred_a - VaR_pred_b)/VaR_true_a;

    if and(isnan(loss(i,1))== 0,(abs(omega(i,1)) < 10))
        i = i + 1;
        err = 0;
    end
    catch
        err = err + 1;
        disp ('Îøèáêà');  
        if err > 10
            i = i + 1;
            err = 0;
        end
    end
end

DisplayStats(omega, alpha, beta, gamma, loss);

nbins = NBINS;
MyHistEl( omega, 0, nbins, '$\bf \omega$')
MyHistEl( alpha, 0, nbins, '$\alpha$')
MyHistEl( beta,  0, nbins, '$\beta$')
MyHistEl( gamma, 0, nbins, '$\gamma$')
MyHistEl( loss,  0, nbins, '$L$')
MyHistEl( VaR,   0, nbins, '$VaR$')


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
VaR    = zeros(S, 1);


err = 0;
i = 1;
while i <= S
    i
    try
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
    b.sigma2 = a.sigma2;
    % True
    a.Estimate();
    [loss_a, VaR_true_a, VaR_pred_a] = a.Predict(p);

    % False    
    b.Estimate();
    [loss_b, VaR_true_b, VaR_pred_b] = b.Predict(p);
    
    omega(i,1) = (a.omega - b.omega)/omega0;
    alpha(i,1) = (a.alpha - b.alpha)/alpha0;
    beta(i,1)  = (a.beta  - b.beta)/beta0;
    gamma(i,1) = (a.gamma - b.gamma)/gamma0;
    loss(i,1)  = (loss_a - loss_b)/QLIKE(a.sigma2(end,1),a.sigma2(end,1));
    VaR(i,1)   = (VaR_pred_a - VaR_pred_b)/VaR_true_a;

    if and(isnan(loss(i,1))== 0,(abs(omega(i,1)) < 10))
        i = i + 1;
        err = 0;
    end
    catch
        err = err + 1;
        disp ('Îøèáêà');  
        if err > 10
            i = i + 1;
            err = 0;
        end
    end
end

DisplayStats(omega, alpha, beta, gamma, loss);

nbins = NBINS;
MyHistEl( omega, 0, nbins, '$\bf \omega$')
MyHistEl( alpha, 0, nbins, '$\alpha$')
MyHistEl( beta,  0, nbins, '$\beta$')
MyHistEl( gamma, 0, nbins, '$\gamma$')
MyHistEl( loss,  0, nbins, '$L$')
MyHistEl( VaR,   0, nbins, '$VaR$')

