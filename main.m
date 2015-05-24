clear all
modelList = {'EGARCH_Const_AR', 'EGARCH_Const_Var', ...
             'EGARCH_General','EGARCH_Ssn', ...
             'GJR_Const_AR','GJR_Const_Var', ...
             'GJR_General','GJR_Ssn'};

modelIndex = 1;
modelNumber = 8;

while (modelIndex <= modelNumber)
    W = load (strcat ('JoinedWorkspace\', modelList{modelIndex}));
    filename = strcat ('graphs\', modelList{modelIndex});
    filename = strcat (filename, '_');
     MainJoinedHist (W.alpha1, W.alpha2, W.alpha3, '$\alpha$', strcat (filename, 'alpha'), 1);
     MainJoinedHist (W.beta1, W.beta2, W.beta3, '$\beta$', strcat (filename, 'beta'), 1);
    MainJoinedHist (W.gamma1, W.gamma2, W.gamma3, '$\gamma$', strcat (filename, 'gamma'), 1);
    MainJoinedHist (W.omega1, W.omega2, W.omega3, '$\omega$', strcat (filename, 'omega'), 1);
    MainJoinedHist (W.VaR1, W.VaR2, W.VaR3, '$VaR$', strcat (filename, 'VaR'), 1);
    MainJoinedHist (W.loss1, W.loss2, W.loss3, '$L$', strcat (filename, 'loss'), 1);
    
    MainJoinedHist (W.alpha1, W.alpha2, W.alpha3, '$\alpha$', strcat (filename, 'alpha_2'), 0);
    MainJoinedHist (W.beta1, W.beta2, W.beta3, '$\beta$', strcat (filename, 'beta_2'), 0);
    MainJoinedHist (W.gamma1, W.gamma2, W.gamma3, '$\gamma$', strcat (filename, 'gamma_2'), 0);
    MainJoinedHist (W.omega1, W.omega2, W.omega3, '$\omega$', strcat (filename, 'omega_2'), 0);
    MainJoinedHist (W.VaR1, W.VaR2, W.VaR3, '$VaR$', strcat (filename, 'VaR_2'), 0);
    MainJoinedHist (W.loss1, W.loss2, W.loss3, '$L$', strcat (filename, 'loss_2'), 0);
    modelIndex = modelIndex + 1;
    close all;
end