clear;
close all

%% Calibration JGR-model
%model = 'gjr';
% model = 'egarch';
T = 2000;
S = 2000;
d = ones(T, 1);
for i = 1:T
    d(i, 1) = mod(i, 5);
    if mod(i, 5) == 0
        d(i, 1) = 5;
    end
end

%% Type 1, gjr
% model = 'gjr';
% nbins = 150;
% MisspType = 0;
% mu0 = -0.147346;
% rho0 = 0.7;
% delta0 = 0;
% theta10 = 0;
% theta20 = 0;
% theta30 = 0;
% theta40 = 0;
% theta50 = 0;
% omega0 = 0.3;
% alpha0 = 0.7;
% beta0 = 0.2;
% gamma0 = 0.07;
% x0 = [mu0, rho0, delta0, theta10, theta20, theta30, theta40, theta50, omega0, alpha0, beta0, gamma0];
% param0 = rho0;
% est_x_t = zeros(S, length(x0));
% est_x_f = zeros(S, length(x0));
% h2_f_t = zeros(S,1);
% h2_f_f = zeros(S,1);
% QLIKE_t = zeros(S,1);
% QLIKE_f = zeros(S,1);
% QLIKE_  = zeros(S,1);
% 
% for i = 1:S
%     [r, h2] = s_mySimulate(model, MisspType, mu0, param0, omega0, alpha0, beta0, gamma0, T+1);
%     [est_x_f(i, :), LogL_f] = myEstimate(model, 0, r(1:T,1), d, param0, omega0, alpha0, beta0, gamma0);
%     [est_x_t(i, :), LogL_t] = s_myEstimate(model, MisspType, r(1:T,1), d, mu0, param0, omega0, alpha0, beta0, gamma0);
%     h2_f = condVariance(model, MisspType, r, d, est_x_f(i, 1), est_x_f(i, 9), est_x_f(i, 10), est_x_f(i, 11), est_x_f(i, 12));
%     h2_t = s_condVariance(model, MisspType, r, d, est_x_t(i, 1), est_x_t(i, 2), est_x_t(i, 9), est_x_t(i, 10), est_x_t(i, 11), est_x_t(i, 12));
%     h2_forecast_f = myForecast(model, r(1:T,1), h2_f, est_x_f(i, :));
%     h2_forecast_t = myForecast(model, r(1:T,1), h2_t, est_x_t(i, :));
%     h2_forecast = myForecast(model, r(1:T,1), h2_t, x0);
%     QLIKE_f(i,1) = QLIKE(r(T+1,1), h2_forecast_f);
%     QLIKE_t(i,1) = QLIKE(r(T+1,1), h2_forecast_t);
%     QLIKE_(i,1)  = QLIKE(r(T+1,1), h2_forecast);
%     i
% end 
% 
% x_omega = (est_x_t(:,9)-est_x_f(:,9))/omega0;
% x_alpha = (est_x_t(:,10)-est_x_f(:,10))/alpha0;
% x_beta = (est_x_t(:,11)-est_x_f(:,11))/beta0;
% x_gamma = (est_x_t(:,12)-est_x_f(:,12))/gamma0;
% x_QLIKE = (QLIKE_t - QLIKE_f)./QLIKE_;
% myDisplay_2(x_omega, x_alpha, x_beta, x_gamma, x_QLIKE);
% 
% MyHistEl( x_omega, 0, nbins, '$\frac{\hat\omega-\hat\omega''}{\omega}$')
% MyHistEl( x_alpha, 0, nbins, '$\frac{\hat\alpha-\hat\alpha''}{\alpha}$')
% MyHistEl( x_beta,  0, nbins, '$\frac{\hat\beta-\hat\beta''}{\beta}$')
% MyHistEl( x_gamma, 0, nbins, '$\frac{\hat\gamma-\hat\gamma''}{\gamma}$')
% MyHistEl( x_QLIKE, 0, nbins, '$\frac{QLIKE_t-QLIKE_f}{QLIKE}$')
% 
%% Type 1, EGARCH
% model = 'egarch';
% nbins = 150;
% MisspType = 0;
% mu0 = 0.03;
% rho0 = 0.7;
% delta0 = 0;
% theta10 = 0;
% theta20 = 0;
% theta30 = 0;
% theta40 = 0;
% theta50 = 0;
% omega0 = 0.3;
% alpha0 = 0.5;
% beta0 = 0.55;
% gamma0 = -0.26;
% x0 = [mu0, rho0, delta0, theta10, theta20, theta30, theta40, theta50, omega0, alpha0, beta0, gamma0];
% param0 = rho0;
% est_x_t = zeros(S, length(x0));
% est_x_f = zeros(S, length(x0));
% h2_f_t = zeros(S,1);
% h2_f_f = zeros(S,1);
% QLIKE_t = zeros(S,1);
% QLIKE_f = zeros(S,1);
% QLIKE_  = zeros(S,1);
% 
% for i = 1:S
%     [r, h2] = s_mySimulate(model, MisspType, mu0, param0, omega0, alpha0, beta0, gamma0, T+1);
%     [est_x_f(i, :), LogL_f] = myEstimate(model, 0, r(1:T,1), d, param0, omega0, alpha0, beta0, gamma0);
%     [est_x_t(i, :), LogL_t] = s_myEstimate(model, MisspType, r(1:T,1), d, mu0, param0, omega0, alpha0, beta0, gamma0);
%     h2_f = condVariance(model, MisspType, r, d, est_x_f(i, 1), est_x_f(i, 9), est_x_f(i, 10), est_x_f(i, 11), est_x_f(i, 12));
%     h2_t = s_condVariance(model, MisspType, r, d, est_x_t(i, 1), est_x_t(i, 2), est_x_t(i, 9), est_x_t(i, 10), est_x_t(i, 11), est_x_t(i, 12));
%     h2_forecast_f = myForecast(model, r(1:T,1), h2_f, est_x_f(i, :));
%     h2_forecast_t = myForecast(model, r(1:T,1), h2_t, est_x_t(i, :));
%     h2_forecast = myForecast(model, r(1:T,1), h2_t, x0);
%     QLIKE_f(i,1) = QLIKE(r(T+1,1), h2_forecast_f);
%     QLIKE_t(i,1) = QLIKE(r(T+1,1), h2_forecast_t);
%     QLIKE_(i,1)  = QLIKE(r(T+1,1), h2_forecast);
%     if mod(i,10) == 0
%         i
%     end
% end 
% 
% x_omega = (est_x_t(:,9)-est_x_f(:,9))/omega0;
% x_alpha = (est_x_t(:,10)-est_x_f(:,10))/alpha0;
% x_beta = (est_x_t(:,11)-est_x_f(:,11))/beta0;
% x_gamma = (est_x_t(:,12)-est_x_f(:,12))/gamma0;
% x_QLIKE = (QLIKE_t - QLIKE_f)./QLIKE_;
% myDisplay_2(x_omega, x_alpha, x_beta, x_gamma, x_QLIKE);
% 
% MyHistEl( x_omega, 0, nbins, '$\frac{\hat\omega-\hat\omega''}{\omega}$')
% MyHistEl( x_alpha, 0, nbins, '$\frac{\hat\alpha-\hat\alpha''}{\alpha}$')
% MyHistEl( x_beta,  0, nbins, '$\frac{\hat\beta-\hat\beta''}{\beta}$')
% MyHistEl( x_gamma, 0, nbins, '$\frac{\hat\gamma-\hat\gamma''}{\gamma}$')
% MyHistEl( x_QLIKE, 0, nbins, '$\frac{QLIKE_t-QLIKE_f}{QLIKE}$')

%% Type 2, gjr
% model = 'gjr';
% nbins = 150;
% S = 2000;
% MisspType = 1;
% mu0 = -0.147346;
% rho0 = 0;
% delta0 = 2;
% theta10 = 0;
% theta20 = 0;
% theta30 = 0;
% theta40 = 0;
% theta50 = 0;
% omega0 = 0.3;
% alpha0 = 0.7;
% beta0 = 0.2;
% gamma0 = 0.07;
% x0 = [mu0, rho0, delta0, theta10, theta20, theta30, theta40, theta50, omega0, alpha0, beta0, gamma0];
% param0 = rho0;
% est_x_t = zeros(S, length(x0));
% est_x_f = zeros(S, length(x0));
% h2_f_t = zeros(S,1);
% h2_f_f = zeros(S,1);
% QLIKE_t = zeros(S,1);
% QLIKE_f = zeros(S,1);
% QLIKE_  = zeros(S,1);
% 
% for i = 1:S
%     [r, h2] = s_mySimulate(model, MisspType, mu0, param0, omega0, alpha0, beta0, gamma0, T+1);
%     [est_x_f(i, :), LogL_f] = myEstimate(model, 0, r(1:T,1), d, param0, omega0, alpha0, beta0, gamma0);
%     [est_x_t(i, :), LogL_t] = s_myEstimate(model, MisspType, r(1:T,1), d, mu0, param0, omega0, alpha0, beta0, gamma0);
%     h2_f = condVariance(model, MisspType, r, d, est_x_f(i, 1), est_x_f(i, 9), est_x_f(i, 10), est_x_f(i, 11), est_x_f(i, 12));
%     h2_t = s_condVariance(model, MisspType, r, d, est_x_t(i, 1), est_x_t(i, 2), est_x_t(i, 9), est_x_t(i, 10), est_x_t(i, 11), est_x_t(i, 12));
%     h2_forecast_f = myForecast(model, r(1:T,1), h2_f, est_x_f(i, :));
%     h2_forecast_t = myForecast(model, r(1:T,1), h2_t, est_x_t(i, :));
%     h2_forecast = myForecast(model, r(1:T,1), h2_t, x0);
%     QLIKE_f(i,1) = QLIKE(r(T+1,1), h2_forecast_f);
%     QLIKE_t(i,1) = QLIKE(r(T+1,1), h2_forecast_t);
%     QLIKE_(i,1)  = QLIKE(r(T+1,1), h2_forecast);
%     if (mod(i, 10) == 0)
%         i
%     end
% end 
% 
% x_omega = (est_x_t(:,9)-est_x_f(:,9))/omega0;
% x_alpha = (est_x_t(:,10)-est_x_f(:,10))/alpha0;
% x_beta = (est_x_t(:,11)-est_x_f(:,11))/beta0;
% x_gamma = (est_x_t(:,12)-est_x_f(:,12))/gamma0;
% x_QLIKE = (QLIKE_t - QLIKE_f)./QLIKE_;
% myDisplay_2(x_omega, x_alpha, x_beta, x_gamma, x_QLIKE);
% 
% MyHistEl( x_omega, 0, nbins, '$\frac{\hat\omega-\hat\omega''}{\omega}$')
% MyHistEl( x_alpha, 0, nbins, '$\frac{\hat\alpha-\hat\alpha''}{\alpha}$')
% MyHistEl( x_beta,  0, nbins, '$\frac{\hat\beta-\hat\beta''}{\beta}$')
% MyHistEl( x_gamma, 0, nbins, '$\frac{\hat\gamma-\hat\gamma''}{\gamma}$')
% MyHistEl( x_QLIKE, 0, nbins, '$\frac{QLIKE_t-QLIKE_f}{QLIKE}$')

%% Type 2, EGARCH
model = 'egarch';
nbins = 150;
S = 2000;
MisspType = 1;
mu0 = 0.01;
rho0 = 0;
delta0 = 2;
theta10 = 0;
theta20 = 0;
theta30 = 0;
theta40 = 0;
theta50 = 0;
omega0 = 0.3;
alpha0 = 0.7;
beta0 = 0.55;
gamma0 = -0.26;
x0 = [mu0, rho0, delta0, theta10, theta20, theta30, theta40, theta50, omega0, alpha0, beta0, gamma0];
param0 = delta0;
est_x_t = zeros(S, length(x0));
est_x_f = zeros(S, length(x0));
h2_f_t = zeros(S,1);
h2_f_f = zeros(S,1);
QLIKE_t = zeros(S,1);
QLIKE_f = zeros(S,1);
QLIKE_  = zeros(S,1);

for i = 1:S
    i
    [r, h2] = s_mySimulate(model, MisspType, mu0, param0, omega0, alpha0, beta0, gamma0, T+1);
    [est_x_f(i, :), LogL_f] = myEstimate(model, 0, r(1:T,1), d, param0, omega0, alpha0, beta0, gamma0);
    [est_x_t(i, :), LogL_t] = s_myEstimate(model, MisspType, r(1:T,1), d, mu0, param0, omega0, alpha0, beta0, gamma0);
    h2_f = condVariance(model, MisspType, r, d, est_x_f(i, 1), est_x_f(i, 9), est_x_f(i, 10), est_x_f(i, 11), est_x_f(i, 12));
    h2_t = s_condVariance(model, MisspType, r, d, est_x_t(i, 1), est_x_t(i, 2), est_x_t(i, 9), est_x_t(i, 10), est_x_t(i, 11), est_x_t(i, 12));
    h2_forecast_f = myForecast(model, r(1:T,1), h2_f, est_x_f(i, :));
    h2_forecast_t = myForecast(model, r(1:T,1), h2_t, est_x_t(i, :));
    h2_forecast = myForecast(model, r(1:T,1), h2_t, x0);
    QLIKE_f(i,1) = QLIKE(r(T+1,1), h2_forecast_f);
    QLIKE_t(i,1) = QLIKE(r(T+1,1), h2_forecast_t);
    QLIKE_(i,1)  = QLIKE(r(T+1,1), h2_forecast);
end 


s_myDisplay(model, x0, est_x_t(end, :), est_x_f(end, :))

x_omega = (est_x_t(:,9)-est_x_f(:,9))/omega0;
x_alpha = (est_x_t(:,10)-est_x_f(:,10))/alpha0;
x_beta = (est_x_t(:,11)-est_x_f(:,11))/beta0;
x_gamma = (est_x_t(:,12)-est_x_f(:,12))/gamma0;
x_QLIKE = (QLIKE_t - QLIKE_f)./QLIKE_;
myDisplay_2(x_omega, x_alpha, x_beta, x_gamma, x_QLIKE);

MyHistEl( x_omega, 0, nbins, '$\frac{\hat\omega-\hat\omega''}{\omega}$')
MyHistEl( x_alpha, 0, nbins, '$\frac{\hat\alpha-\hat\alpha''}{\alpha}$')
MyHistEl( x_beta,  0, nbins, '$\frac{\hat\beta-\hat\beta''}{\beta}$')
MyHistEl( x_gamma, 0, nbins, '$\frac{\hat\gamma-\hat\gamma''}{\gamma}$')
MyHistEl( x_QLIKE, 0, nbins, '$\frac{QLIKE_t-QLIKE_f}{QLIKE}$')

%% Type 3
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
% est_x_t = zeros(S, length(x0));
% est_x_f = zeros(S, length(x0));
% h2_f_t = zeros(S,1);
% h2_f_f = zeros(S,1);
% QLIKE_t = zeros(S,1);
% QLIKE_f = zeros(S,1);
% QLIKE_  = zeros(S,1);
% 
% for i = 1:S
%     [r, h2] = s_mySimulate(model, MisspType, mu0, param0, omega0, alpha0, beta0, gamma0, T+1);
%     [est_x_f(i, :), LogL_f] = myEstimate(model, 0, r(1:T,1), d, param0, omega0, alpha0, beta0, gamma0);
%     [est_x_t(i, :), LogL_t] = s_myEstimate(model, MisspType, r(1:T,1), d, mu0, param0, omega0, alpha0, beta0, gamma0);
%     h2_f = condVariance(model, MisspType, r, d, est_x_f(i, 1), est_x_f(i, 9), est_x_f(i, 10), est_x_f(i, 11), est_x_f(i, 12));
%     h2_t = s_condVariance(model, MisspType, r, d, est_x_t(i, 1), est_x_t(i, 2), est_x_t(i, 9), est_x_t(i, 10), est_x_t(i, 11), est_x_t(i, 12));
%     h2_forecast_f = myForecast(model, r(1:T,1), h2_f, est_x_f(i, :));
%     h2_forecast_t = myForecast(model, r(1:T,1), h2_t, est_x_t(i, :));
%     h2_forecast = myForecast(model, r(1:T,1), h2_t, x0);
%     QLIKE_f(i,1) = QLIKE(r(T+1,1), h2_forecast_f);
%     QLIKE_t(i,1) = QLIKE(r(T+1,1), h2_forecast_t);
%     QLIKE_(i,1)  = QLIKE(r(T+1,1), h2_forecast);
% end 
% 
% 
% 
% MyHistEl( (est_x_t(:,9)-est_x_f(:,9))/omega0, 0, nbins, '$\omega$')
% MyHistEl( (est_x_t(:,10)-est_x_f(:,10))/alpha0, 0, nbins, '$\alpha$')
% MyHistEl( (est_x_t(:,11)-est_x_f(:,11))/beta0, 0, nbins, '$\beta$')
% MyHistEl( (est_x_t(:,12)-est_x_f(:,12))/gamma0, 0, nbins, '$\gamma$')
% MyHistEl( (QLIKE_t - QLIKE_f)./QLIKE_, 0, nbins, '$loss$')
