clear all
close all

%--------------------------------------------------------------------------
% Terminal case + *10
%--------------------------------------------------------------------------

% MUMUMU = 0.001;
% RHORHORHO = 0.15;
% DELTADELTADELTA = 3.303;
% THETHA01 = -0.000583391;
% THETHA02 = 0.00143266434255308;
% THETHA03 = 0.00224614578474341;
% THETHA04 = 0.000154878290132809;
% THETHA05 = 0.0015718522946019;
% OMEGA = 0.000004;
% ALPHA = 0.033;
% BETTA = 0.9;
% GAMMA = 0.08;
% SIMULATION_NUMBER = 2000;
% NBINS = 200;
% p = 0.05;


%--------------------------------------------------------------------------
% Normal case
%--------------------------------------------------------------------------

MUMUMU = 4.7562e-04;
RHORHORHO = 0.0733;
DELTADELTADELTA = 1.6122;
THETHA01 = -3.1013e-06;
THETHA02 = 2.9867e-04;
THETHA03 = 7.2473e-04;
THETHA04 = 2.9180e-04;
THETHA05 = 0.0011;
OMEGA = 0.000004;
ALPHA = 0.033;
BETTA = 0.9;
GAMMA = 0.08;
SIMULATION_NUMBER = 2000;
NBINS = 200;
p = 0.05;
%--------------------------------------------------------------------------




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

%--------------------------------------------------------------------------
% Terminal case + *10
%--------------------------------------------------------------------------

% MUMUMU = 0.001;
% RHORHORHO = 0.15;
% DELTADELTADELTA = 3.303;
% THETHA01 = -0.000583391;
% THETHA02 = 0.0014;
% THETHA03 = 0.0022;
% THETHA04 = 0.0002;
% THETHA05 = 0.0016;
% OMEGA = -0.8;
% ALPHA = 0.25;
% BETTA = 0.93;
% GAMMA = -0.12;


%--------------------------------------------------------------------------
% Normal case
%--------------------------------------------------------------------------

MUMUMU = 4.3409e-04;
RHORHORHO = 0.0733;
DELTADELTADELTA = 1.6122;
THETHA01 = 1.5964e-05;
THETHA02 = 2.6119e-04;
THETHA03 = 6.9964e-04;
THETHA04 = 2.6202e-04;
THETHA05 = 9.5764e-046;
OMEGA = -0.8;
ALPHA = 0.25;
BETTA = 0.93;
GAMMA = -0.12;


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