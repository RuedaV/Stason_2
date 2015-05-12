function [ h2_f ] = myForecast(model, r, h2, x)
mu_ = x(1,1);
rho = x(1,2);
delta = x(1,3);
theta1 = x(1,4);
theta2 = x(1,5);
theta3 = x(1,6);
theta4 = x(1,7);
theta5 = x(1,8);
omega = x(1,9);
alpha = x(1,10);
beta = x(1,11);
gamma = x(1,12);
d = mod(length(x),5);
w1 = (d-2)*(d-3)*(d-4)*(d-5)/((1-2)*(1-3)*(1-4)*(1-5));
w2 = (d-1)*(d-3)*(d-4)*(d-5)/((2-1)*(2-3)*(2-4)*(2-5));
w3 = (d-1)*(d-2)*(d-4)*(d-5)/((3-1)*(3-2)*(3-4)*(3-5));
w4 = (d-1)*(d-2)*(d-3)*(d-5)/((4-1)*(4-2)*(4-3)*(4-5));
w5 = (d-1)*(d-2)*(d-3)*(d-4)/((5-1)*(5-2)*(5-3)*(5-4));
e = r(end, 1) - mu_ - rho*r(end-1, 1) - delta*h2(end,1) - theta1*w1 - theta2*w2 - theta3*w3 - theta4*w4 - theta5*w5;

if (strcmp(model, 'gjr') == 1)    
    h2_f = omega + alpha*e^2 + beta*h2(end,1) + gamma*(e < 0);
elseif (strcmp(model, 'egarch') == 1) 
    h2_f = exp(omega + alpha*abs(e)/sqrt(h2(end,1)) + beta*log(h2(end,1)) + gamma*e/sqrt(h2(end,1)));
else
    error('choose an appropriate model!')
end

end

