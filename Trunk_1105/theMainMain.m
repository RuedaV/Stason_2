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
SIMULATION_NUMBER = 15;
NBINS = 50;
p = 0.05;

mu0 = MUMUMU;
rho0 = RHORHORHO;
omega0 = OMEGA;
alpha0 = ALPHA;
beta0 = BETTA;
gamma0 = GAMMA;

S = SIMULATION_NUMBER;
T = 10;

model1 = GJR_Const_AR(mu0, rho0, omega0, alpha0, beta0, gamma0);
model2 = GJR_Const(mu0, omega0, alpha0, beta0, gamma0);

theMainIteration(model1, model2, 0.05, S, T, NBINS, 'graphs\');
