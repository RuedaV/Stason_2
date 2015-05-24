function [] = theMainIteration (MODEL1, MODEL2, p, S, T, NBINS, PATH, season)

omega0 = MODEL1.omega;
alpha0 = MODEL1.alpha;
beta0 = MODEL1.beta;
gamma0 = MODEL1.gamma;

omega  = zeros(S, 1);
alpha  = zeros(S, 1);
beta   = zeros(S, 1);
gamma  = zeros(S, 1);
loss   = zeros(S, 1);
VaR    = zeros(S, 1);
VaR_a  = zeros(S, 1);
VaR_b  = zeros(S, 1);
r_last = zeros(S, 1);

iteration_name = class(MODEL1);

err = 0;
i = 1;
while i <= S
    disp (strcat (iteration_name, int2str(i)));
    try
        a  = copy(MODEL1);
        b  = copy(MODEL2);
        if season
            [Data,Day] = a.Simulate(T);    
        else
            Data = a.Simulate(T);
        end
        a.data = Data(1:end-1);
        a.data_plus = Data;
        if season
            a.day = Day(1:end-1);
            a.day_plus = Day;
        end
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

        if and(isnan(loss(i,1))== 0,(abs(omega(i,1)) < 20))
            i = i + 1;
            err = 0;
        end
    catch
        err = err + 1;
        if err > 10
            i = i + 1;
            err = 0;
        end
    end
end
try 
    DisplayStats(omega, alpha, beta, gamma, loss);

    nbins = NBINS;
    MainHist( omega, nbins, '$\omega$', strcat (PATH, class(MODEL1), '_omega'))
    MainHist( alpha, nbins, '$\alpha$', strcat (PATH, class(MODEL1), '_alpha'))
    MainHist( beta,  nbins, '$\beta$', strcat (PATH, class(MODEL1), '_beta'))
    MainHist( gamma, nbins, '$\gamma$', strcat (PATH, class(MODEL1), '_gamma'))
    MainHist( loss,  nbins, '$L$', strcat (PATH, class(MODEL1), '_loss'))
    MainHist( VaR,   nbins, '$VaR$',strcat (PATH, class(MODEL1), '_var'))
catch
end
try 
    QOmega = quantile(omega, 0.5) / (quantile(omega, 0.75) - quantile(omega, 0.25));
    QAlpha = quantile(alpha, 0.5) / (quantile(alpha, 0.75) - quantile(alpha, 0.25));
    QBeta = quantile(beta, 0.5) / (quantile(beta, 0.75) - quantile(beta, 0.25));
    QGamma = quantile(gamma, 0.5) / (quantile(gamma, 0.75) - quantile(gamma, 0.25));
    QLoss = quantile(loss, 0.5) / (quantile(loss, 0.75) - quantile(loss, 0.25));
    QVar = quantile(VaR, 0.5) / (quantile(VaR, 0.75) - quantile(VaR, 0.25));
catch
end   

    save (strcat ('workspace\', iteration_name));

