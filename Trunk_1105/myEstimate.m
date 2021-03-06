function [est_x, LogL] = myEstimate(model, misspType, r, d, param0, omega0, alpha0, beta0, gamma0)

if (strcmp(model, 'gjr') == 1)
    s = 1e-5;
    x0 = [param0, omega0, alpha0, beta0, gamma0];
    A1 = [zeros(size(param0)),     -1,      0,     0,      0]; % omega > 0
    A2 = [zeros(size(param0)),      0,     -1,     0,      0]; % alpha >= 0
    A3 = [zeros(size(param0)),      0,      0,    -1,      0]; % beta >= 0
    A4 = [zeros(size(param0)),      0,     -1,     0,     -1]; % alpha + gamma >= 0

    Aineq = [A1; A2; A3; A4];
    bineq = [-s;  0;  0;  0];


    if (misspType == 0)

        [temp, LogL] = fmincon(@(x)(-logL(model, misspType, r, d, x)),x0,Aineq,bineq,[],[],[],[],[], optimset('Algorithm','interior-point','Hessian','bfgs', 'Display', 'off'));
        est_x = [temp(1), 0, 0, 0, 0, 0, 0, 0, temp(2), temp(3), temp(4), temp(5)];

    elseif (misspType == 1)

        [temp, LogL] = fmincon(@(x)(-logL(model, misspType, r, d, x)),x0,Aineq,bineq,[],[],[],[],[], optimset('Algorithm','interior-point','Hessian','bfgs', 'Display', 'off'));
        est_x = [0, temp(1), 0, 0, 0, 0, 0, 0, temp(2), temp(3), temp(4), temp(5)];

    elseif (misspType == 2)

        [temp, LogL] = fmincon(@(x)(-logL(model, misspType, r, d, x)),x0,Aineq,bineq,[],[],[],[],[], optimset('Algorithm','interior-point','Hessian','bfgs', 'Display', 'off'));
        est_x = [0, 0, temp(1), 0, 0, 0, 0, 0, temp(2), temp(3), temp(4), temp(5)];

    elseif (misspType == 3)    

        [temp, LogL] = fmincon(@(x)(-logL(model, misspType, r, d, x)),x0,Aineq,bineq,[],[],[],[],[], optimset('Algorithm','interior-point','Hessian','bfgs', 'Display', 'off'));
        %       mu  rho  delta   theta1   theta2   theta3   theta4   theta5   omega     alpha     beta     gamma
        est_x = [0,   0,     0, temp(1), temp(2), temp(3), temp(4), temp(5), temp(6), temp(7), temp(8), temp(9)];
    end
    LogL = - LogL;
    
elseif (strcmp(model, 'egarch') == 1)
    x0 = [param0, omega0, alpha0, beta0, gamma0];
    A1 = [zeros(size(param0)),     -1,      0,     0,      0]; % omega > -100


    Aineq = A1;
    bineq = 100;

    if (misspType == 0)

        [temp, LogL] = fmincon(@(x)(-logL(model, misspType, r, d, x)),x0,Aineq,bineq,[],[],[],[],[], optimset('Algorithm','interior-point','Hessian','bfgs', 'Display', 'off'));
        est_x = [temp(1), 0, 0, 0, 0, 0, 0, 0, temp(2), temp(3), temp(4), temp(5)];

    elseif (misspType == 1)

        [temp, LogL] = fmincon(@(x)(-logL(model, misspType, r, d, x)),x0,Aineq,bineq,[],[],[],[],[], optimset('Algorithm','interior-point','Hessian','bfgs', 'Display', 'off'));
        est_x = [0, temp(1), 0, 0, 0, 0, 0, 0, temp(2), temp(3), temp(4), temp(5)];

    elseif (misspType == 2)

        [temp, LogL] = fmincon(@(x)(-logL(model, misspType, r, d, x)),x0,Aineq,bineq,[],[],[],[],[], optimset('Algorithm','interior-point','Hessian','bfgs', 'Display', 'off'));
        est_x = [0, 0, temp(1), 0, 0, 0, 0, 0, temp(2), temp(3), temp(4), temp(5)];

    elseif (misspType == 3)    

        [temp, LogL] = fmincon(@(x)(-logL(model, misspType, r, d, x)),x0,Aineq,bineq,[],[],[],[],[], optimset('Algorithm','interior-point','Hessian','bfgs', 'Display', 'off'));
        %       mu  rho  delta   theta1   theta2   theta3   theta4   theta5   omega     alpha     beta     gamma
        est_x = [0,   0,     0, temp(1), temp(2), temp(3), temp(4), temp(5), temp(6), temp(7), temp(8), temp(9)];
    end

    LogL = - LogL;
else
    error('choose an appropriate model!');
end


end

