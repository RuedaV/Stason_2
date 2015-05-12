function y = statistics(x)

q05      = quantile(x, 0.05);
q25      = quantile(x, 0.25);
q50      = quantile(x, 0.50);
q75      = quantile(x, 0.75);
q95      = quantile(x, 0.95);
y = zeros(1, 5);
y = [q05, q25, q50, q75, q95];

end

