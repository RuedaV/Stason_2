clear;
close all

%% Calibration JGR-model
model = 'gjr';
%model = 'egarch';
window = 2000;
T = 10000;
d = ones(T, 1);
for i = 1:T
    d(i, 1) = mod(i, 5);
end

% Type 1
S = 10;
nbins = 3;
MisspType = 0;
mu0 = 0.01;
rho0 = 0.01;
delta0 = 0.3;
theta10 = 0;
theta20 = 0;
theta30 = 0;
theta40 = 0;
theta50 = 0;
omega0 = 0.3;
alpha0 = 0.2;
beta0 = 0.7;
gamma0 = 0.06;
x0 = [mu0, rho0, delta0, theta10, theta20, theta30, theta40, theta50, omega0, alpha0, beta0, gamma0];
param0 = rho0;
est_x_t = zeros(S, length(x0));
est_x_f = zeros(S, length(x0));

for i = 1:S
    [r, h2] = s_mySimulate(model, MisspType, mu0, param0, omega0, alpha0, beta0, gamma0, T);
    [est_x_t(i, :), LogL_t] = myEstimate(model, 0, r, d, param0, omega0, alpha0, beta0, gamma0);
    [est_x_f(i, :), LogL_f] = s_myEstimate(model, MisspType, r, d, mu0, param0, omega0, alpha0, beta0, gamma0);
end 

relative = (est_x_t - est_x_f);
Text = 'bla';
MyHistEl( (est_x_t(:,1)-est_x_f(:,1))/x0(1), 0, nbins, Text)
 

% % Type 2
% MisspType = 1;
% mu0 = 0.01;
% rho0 = 0;
% delta0 = 0.3;
% theta10 = 0;
% theta20 = 0;
% theta30 = 0;
% theta40 = 0;
% theta50 = 0;
% omega0 = 0.3;
% alpha0 = 0.2;
% beta0 = 0.7;
% gamma0 = 0.06;
% x0 = [mu0, rho0, delta0, theta10, theta20, theta30, theta40, theta50, omega0, alpha0, beta0, gamma0];
% param0 = delta0;
% [r, h2] = s_mySimulate(model, MisspType, mu0, param0, omega0, alpha0, beta0, gamma0, T);
% [est_x_t, LogL_t] = myEstimate(model, 0, r, d, param0, omega0, alpha0, beta0, gamma0);
% [est_x_f, LogL_f] = s_myEstimate(model, MisspType, r, d, mu0, param0, omega0, alpha0, beta0, gamma0);

% % Type 3
% MisspType = 2;
% mu0 = 0;
% rho0 = 0;
% delta0 = 0;
% theta10 = 0.2;
% theta20 = 0.1;
% theta30 = 0.5;
% theta40 = 0.5;
% theta50 = 0.5;
% omega0 = 0.3;
% alpha0 = 0.2;
% beta0 = 0.7;
% gamma0 = 0.06;
% x0 = [mu0, rho0, delta0, theta10, theta20, theta30, theta40, theta50, omega0, alpha0, beta0, gamma0];
% param0 = [theta10, theta20, theta30, theta40, theta50];
% 
% [r, h2] = s_mySimulate(model, MisspType, mu0, param0, omega0, alpha0, beta0, gamma0, T);
% [est_x_t, LogL_t] = myEstimate(model, 0, r, d, mu0, omega0, alpha0, beta0, gamma0);
% [est_x_f, LogL_f] = myEstimate(model, 3, r, d, param0, omega0, alpha0, beta0, gamma0);



% plot(r)
% plot(h2)
s_myDisplay(model, x0, est_x_t, est_x_f)
