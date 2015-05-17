clear all
close all

MUMUMU = 0.001;
RHORHORHO = 0.009*10;
DELTADELTADELTA = 3.303*10;
THETHA01 = -0.000583391*10;
THETHA02 = 0.00143266434255308*10;
THETHA03 = 0.00224614578474341*10;
THETHA04 = 0.000154878290132809*10;
THETHA05 = 0.0015718522946019*10;
OMEGA = 0.000004;
ALPHA = 0.033;
BETTA = 0.9;
GAMMA = 0.08;
SIMULATION_NUMBER = 2000;
NBINS = 200;
p = 0.05;

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

S = SIMULATION_NUMBER;
T = 1000;
    
a  = GJR_Const_AR(mu0, rho0, omega0, alpha0, beta0, gamma0);
b  = GJR_Const(mu0, omega0, alpha0, beta0, gamma0);
theMainIteration(a, b, p, S, T, NBINS, 'graphs\', 0);

a  = GJR_Const_Var(mu0, delta0, omega0, alpha0, beta0, gamma0);
theMainIteration(a, b, p, S, T, NBINS, 'graphs\', 0);

a  = GJR_Ssn(theta10, theta20, theta30, theta40, theta50,omega0, alpha0, beta0, gamma0);
theMainIteration(a, b, p, S, T, NBINS, 'graphs\', 1);

a  = GJR_General(rho0, delta0, theta10, theta20, theta30, theta40, theta50,...
    omega0, alpha0, beta0, gamma0);
theMainIteration(a, b, p, S, T, NBINS, 'graphs\', 1);

%% EGARCH
SIMULATION_NUMBER = 2000;
NBINS = 200;
p = 0.05;

MUMUMU = 0.001;
RHORHORHO = 0.15*10;
DELTADELTADELTA = 3.303*10;
OMEGA = -0.8;
ALPHA = 0.25;
BETTA = 0.93;
GAMMA = -0.12;
THETHA01 = -0.000583391*10;
THETHA02 = 0.0014*10;
THETHA03 = 0.0022*10;
THETHA04 = 0.0002*10;
THETHA05 = 0.0016*10;

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

a  = EGARCH_Const_Var(mu0, delta0, omega0, alpha0, beta0, gamma0);
b  = EGARCH_Const(mu0, omega0, alpha0, beta0, gamma0);

theMainIteration(a, b, p, S, T, NBINS, 'graphs\', 0);

a  = EGARCH_Const_AR(mu0, rho0, omega0, alpha0, beta0, gamma0);
theMainIteration(a, b, p, S, T, NBINS, 'graphs\', 0);

a  = EGARCH_Ssn(theta10, theta20, theta30, theta40, theta50,omega0, alpha0, beta0, gamma0);
theMainIteration(a, b, p, S, T, NBINS, 'graphs\', 1);

a  = EGARCH_General(rho0, delta0, theta10, theta20, theta30, theta40, theta50,...
    omega0, alpha0, beta0, gamma0);
theMainIteration(a, b, p, S, T, NBINS, 'graphs\', 1);