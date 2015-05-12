function [] = s_myDisplay(model, x0, est_x_t, est_x_f)

if strcmp(model, 'gjr') == 1
    fprintf('%6s %12s\n\n', 'model', 'GJR-GARCH');
elseif strcmp(model, 'egarch') == 1
    fprintf('%6s %12s\n\n', 'model', 'EGARCH');
else
    error('Wrong name of the model! Try gjr or egarch');
end


fprintf('%6s %12s %12s %12s\r','param','true', 'from true', 'from false');

fprintf('%6s %12.6f %12.6f %12.6f\n', 'mu',     x0(1), est_x_t(1), est_x_f(1));
fprintf('%6s %12.6f %12.6f %12.6f\n', 'rho',    x0(2), est_x_t(2), est_x_f(2));
fprintf('%6s %12.6f %12.6f %12.6f\n', 'delta',  x0(3), est_x_t(3), est_x_f(3));
fprintf('%6s %12.6f %12.6f %12.6f\n', 'theta1', x0(4), est_x_t(4), est_x_f(4));
fprintf('%6s %12.6f %12.6f %12.6f\n', 'theta2', x0(5), est_x_t(5), est_x_f(5));
fprintf('%6s %12.6f %12.6f %12.6f\n', 'theta3', x0(6), est_x_t(6), est_x_f(6));
fprintf('%6s %12.6f %12.6f %12.6f\n', 'theta4', x0(7), est_x_t(7), est_x_f(7));
fprintf('%6s %12.6f %12.6f %12.6f\n', 'theta5', x0(8), est_x_t(8), est_x_f(8));
fprintf('%6s %12.6f %12.6f %12.6f\n', 'omega',  x0(9), est_x_t(9), est_x_f(9));
fprintf('%6s %12.6f %12.6f %12.6f\n', 'alpha',  x0(10), est_x_t(10), est_x_f(10));
fprintf('%6s %12.6f %12.6f %12.6f\n', 'beta',   x0(11), est_x_t(11), est_x_f(11));
fprintf('%6s %12.6f %12.6f %12.6f\n', 'gamma',  x0(12), est_x_t(12), est_x_f(12));


end

