function [] = theMainIteration (MODEL1, MODEL2, p, S, T, NBINS, PATH)

omega  = zeros(S, 1);
alpha  = zeros(S, 1);
beta   = zeros(S, 1);
gamma  = zeros(S, 1);
loss   = zeros(S, 1);
VaR_a  = 0;
VaR_b  = 0;
err = 0;

omega0 = MODEL1.omega;
alpha0 = MODEL1.alpha;
beta0  = MODEL1.beta;
gamma0 = MODEL1.gamma;

iteration_name = class(MODEL1);

iteration = 1;
while iteration <= S
    disp (strcat (iteration_name, int2str(iteration)));
        % Making a copy from templates
    a  = copy (MODEL1);
    b  = copy (MODEL2);

    try
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
    catch
        err = err + 1;
        disp ('������');  
        if err > 10
            iteration = iteration + 1;
            err = 0;
        end
        continue
    end
    
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
end

DisplayStats(omega, alpha, beta, gamma, loss);

nbins = NBINS;
MainHist( omega, nbins, '$\frac{\hat\omega - \hat\omega''}{\omega}$', strcat (PATH, class(MODEL1), '_omega'))
MainHist( alpha, nbins, '$\frac{\hat\alpha - \hat\alpha''}{\alpha}$', strcat (PATH, class(MODEL1), '_alpha'))
MainHist( beta,  nbins, '$\frac{\hat\beta - \hat\beta''}{\beta}$', strcat (PATH, class(MODEL1), '_beta'))
MainHist( gamma, nbins, '$\frac{\hat\gamma - \hat\gamma''}{\gamma}$', strcat (PATH, class(MODEL1), '_gamma'))
MainHist( loss,  nbins, '$\frac{\hat L - \hat L''}{L}$', strcat (PATH, class(MODEL1), '_loss'))

p_value_a = Xi_squared(p, S,  VaR_exceeded_a);
p_value_b = Xi_squared(p, S,  VaR_exceeded_b);

fprintf('\n%6s %12s \r', 'Model', 'p-value');
fprintf('%6s %12.3f \n', 'True',    p_value_a);
fprintf('%6s %12.3f \n', 'False',   p_value_a);
