function p_value = Xi_squared(p, SIMULATION_NUMBER,  VaR_exceeded)

x = VaR_exceeded/SIMULATION_NUMBER
LR = 2*log((x/p)^VaR_exceeded*((1-x)/(1-p))^(SIMULATION_NUMBER-VaR_exceeded))
p_value = 1-chi2cdf(LR,1);
end

