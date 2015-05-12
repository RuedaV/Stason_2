function [Aeq, beq]  = buildEqConstr(mu_, rho, delta, omega, alpha, beta, gamma)

x = [mu_, rho, delta, omega, alpha, beta, gamma];
ind = zeros(size(x));

for i = 1:length(x)
    if x(i) == 0
        ind(i) = 1;
    end
end

m = sum(ind);
Aeq = zeros(m, length(x));
beq = zeros(m, 1);

j = 1;
for i = 1:length(x)
    if ind(i) == 1
        Aeq(j, i) = 1;
        j = j+1;
    end
end       

end