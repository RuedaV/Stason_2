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

mu0 = MUMUMU;
rho0 = RHORHORHO;
omega0 = OMEGA;
alpha0 = ALPHA;
beta0 = BETTA;
gamma0 = GAMMA;

S = SIMULATION_NUMBER;
T = 1000;

omega  = zeros(S, 1);
alpha  = zeros(S, 1);
beta   = zeros(S, 1);
gamma  = zeros(S, 1);
loss   = zeros(S, 1);
VaR_a  = 0;
VaR_b  = 0;
err = 0;

MODEL1TEMPLATE = GJR_Const_AR(mu0, rho0, omega0, alpha0, beta0, gamma0);
MODEL2TEMPLATE = GJR_Const(mu0, omega0, alpha0, beta0, gamma0);

iteration = 1;
while iteration <= S
    disp (iteration);
    try
        % Making a copy from templates
        a  = copy (MODEL1TEMPLATE);
        b  = copy (MODEL2TEMPLATE);

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
        %a.MyDisplay();
        
        % False    
        b.Estimate();
        [loss_b, loss2_b, VaR_exceeded_b] = b.Predict(p);
        VaR_b = VaR_b + VaR_exceeded_b;
        %b.MyDisplay();
        omega(iteration,1) = (a.omega - b.omega)/omega0;
        alpha(iteration,1) = (a.alpha - b.alpha)/alpha0;
        beta(iteration,1)  = (a.beta  - b.beta)/beta0;
        gamma(iteration,1) = (a.gamma - b.gamma)/gamma0;
        loss(iteration,1)  = (loss_a - loss_b)/QLIKE(a.sigma2(end,1),a.sigma2(end,1));

        if and(isnan(loss(iteration,1))== 0, (abs(omega(iteration,1)) < 10))
            iteration = iteration + 1;
            err = 0;
        end
    catch
        err = err + 1;
        disp ('Îøèáêà');  
        if err > 10
            iteration = iteration + 1;
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
