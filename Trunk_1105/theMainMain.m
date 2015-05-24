clear all
close all

NBINS = 200;
p = 0.05;
SIMULATION_NUMBER = 10;
S = SIMULATION_NUMBER;
T = 1000;

[MUMUMU, RHORHORHO, DELTADELTADELTA, THETHA01, ...
    THETHA02, THETHA03, THETHA04, THETHA05, OMEGA, ALPHA, ...
    BETTA, GAMMA, FOLDERNAME] = GetParams(1, 1, 1);
                              % GetParams(IsGJR, IsNormal, Coef), 
                              % GJR Extreme x 10 = GetParams(1, 0, 10)
mu0 = MUMUMU;
rho0 = RHORHORHO;
theta10 = THETHA01;
theta20 = THETHA02;
theta30 = THETHA03;
theta40 = THETHA04;
theta50 = THETHA05;
omega0 = OMEGA;
alpha0 = ALPHA;
beta0 = BETTA;
gamma0 = GAMMA;
delta0 = DELTADELTADELTA;
    
a  = GJR_Const_AR(mu0, rho0, omega0, alpha0, beta0, gamma0);
b  = GJR_Const(mu0, omega0, alpha0, beta0, gamma0);
theMainIteration(a, b, p, S, T, NBINS, 'graphs\', 0, FOLDERNAME);

a  = GJR_Const_Var(mu0, delta0, omega0, alpha0, beta0, gamma0);
theMainIteration(a, b, p, S, T, NBINS, 'graphs\', 0, FOLDERNAME);

a  = GJR_Ssn(theta10, theta20, theta30, theta40, theta50,omega0, alpha0, beta0, gamma0);
theMainIteration(a, b, p, S, T, NBINS, 'graphs\', 1, FOLDERNAME);

a  = GJR_General(rho0, delta0, theta10, theta20, theta30, theta40, theta50,...
    omega0, alpha0, beta0, gamma0);
theMainIteration(a, b, p, S, T, NBINS, 'graphs\', 1, FOLDERNAME);

%% EGARCH

[MUMUMU, RHORHORHO, DELTADELTADELTA, THETHA01, ...
    THETHA02, THETHA03, THETHA04, THETHA05, OMEGA, ALPHA, ...
    BETTA, GAMMA, FOLDERNAME] = GetParams(0, 1, 1);
                              % GetParams(IsGJR, IsNormal, Coef)
mu0 = MUMUMU;
rho0 = RHORHORHO;
omega0 = OMEGA;
alpha0 = ALPHA;
beta0 = BETTA;
gamma0 = GAMMA;
delta0 = DELTADELTADELTA;
theta10 = THETHA01;
theta20 = THETHA02;
theta30 = THETHA03;
theta40 = THETHA04;
theta50 = THETHA05;

b  = EGARCH_Const(mu0, omega0, alpha0, beta0, gamma0);

a  = EGARCH_Const_AR(mu0, rho0, omega0, alpha0, beta0, gamma0);
theMainIteration(a, b, p, S, T, NBINS, 'graphs\', 0, FOLDERNAME);

a  = EGARCH_Const_Var(mu0, delta0, omega0, alpha0, beta0, gamma0);
theMainIteration(a, b, p, S, T, NBINS, 'graphs\', 0, FOLDERNAME);

a  = EGARCH_Ssn(theta10, theta20, theta30, theta40, theta50,omega0, alpha0, beta0, gamma0);
theMainIteration(a, b, p, S, T, NBINS, 'graphs\', 1, FOLDERNAME);

a  = EGARCH_General(rho0, delta0, theta10, theta20, theta30, theta40, theta50,...
    omega0, alpha0, beta0, gamma0);
theMainIteration(a, b, p, S, T, NBINS, 'graphs\', 1, FOLDERNAME);