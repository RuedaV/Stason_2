clear all
workspaceList = {'EGARCH_Const_AR', 'EGARCH_Const_Var', ...
             'EGARCH_General','EGARCH_Ssn', ...
             'GJR_Const_AR','GJR_Const_Var', ...
             'GJR_General','GJR_Ssn'};
workspaceIndex = 1;
workspaceNumber = 8;

while (workspaceIndex <= workspaceNumber)
    W1 = load (strcat ('workspaces\workspace_1\', workspaceList{workspaceIndex}));
    W2 = load (strcat ('workspaces\workspace_2\', workspaceList{workspaceIndex}));
    W3 = load (strcat ('workspaces\workspace_3\', workspaceList{workspaceIndex}));
    
    alpha1 = W1.('alpha');
    alpha2 = W2.('alpha');
    alpha3 = W3.('alpha');
    
    beta1 = W1.('beta');
    beta2 = W2.('beta');
    beta3 = W3.('beta');
    
    gamma1 = W1.('gamma');
    gamma2 = W2.('gamma');
    gamma3 = W3.('gamma');
    
    omega1 = W1.('omega');
    omega2 = W2.('omega');
    omega3 = W3.('omega');
    
    loss1 = W1.('loss');
    loss2 = W2.('loss');
    loss3 = W3.('loss');
    
    VaR1 = W1.('VaR');
    VaR2 = W2.('VaR');
    VaR3 = W3.('VaR');
    
    save (strcat ('JoinedWorkspace\', workspaceList{workspaceIndex}));
    workspaceIndex = workspaceIndex + 1;
end