function [MUMUMU, RHORHORHO, DELTADELTADELTA, THETHA01, THETHA02, THETHA03, THETHA04, THETHA05, OMEGA, ALPHA, ...
    BETTA, GAMMA, FOLDERNAME] = GetParams(IsGJR, IsNormal, Coef)

if (IsGJR == 1) 
    if (IsNormal == 1)
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
        FOLDERNAME = 'workspace\workspace_1\';
    else 
        MUMUMU = 0.001;
        RHORHORHO = 0.5;
        DELTADELTADELTA = 3.303*Coef;
        THETHA01 = -0.000583391*Coef;
        THETHA02 = 0.00143266434255308*Coef;
        THETHA03 = 0.00224614578474341*Coef;
        THETHA04 = 0.000154878290132809*Coef;
        THETHA05 = 0.0015718522946019*Coef;
        OMEGA = 0.000004;
        ALPHA = 0.033;
        BETTA = 0.9;
        GAMMA = 0.08;

        if (Coef == 1)
            FOLDERNAME = 'workspace\workspace_2\';
        else
            FOLDERNAME = 'workspace\workspace_3\';
        end
    end
    
else
    % EGARCH
     if (IsNormal == 1)
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

        FOLDERNAME = 'workspace\workspace_1\';
     else
         
        MUMUMU = 0.001;
        RHORHORHO = 0.5;
        DELTADELTADELTA = 3.303*Coef;
        THETHA01 = -0.000583391*Coef;
        THETHA02 = 0.0014*Coef;
        THETHA03 = 0.0022*Coef;
        THETHA04 = 0.0002*Coef;
        THETHA05 = 0.0016*Coef;
        OMEGA = -0.8;
        ALPHA = 0.25;
        BETTA = 0.93;
        GAMMA = -0.12;
         
        if (Coef == 1)
            FOLDERNAME = 'workspace\workspace_2\';
        else
            FOLDERNAME = 'workspace\workspace_3\';
        end
     end
end