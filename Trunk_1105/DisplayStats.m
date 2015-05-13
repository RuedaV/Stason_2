function DisplayStats(omega, alpha, beta, gamma, loss)
stats = zeros(6, 5);
stats(1,:) = statistics(omega);
stats(2,:) = statistics(alpha);
stats(3,:) = statistics(beta);
stats(4,:) = statistics(gamma);
stats(5,:) = statistics(loss);

fprintf('%6s %10s %10s %10s %10s %10s\n',  'param', 'q05', 'q25', 'median', 'q75', 'q95');
fprintf('%6s %10.3f %10.3f %10.3f %10.3f %10.3f\n', 'omega',stats(1,:));
fprintf('%6s %10.3f %10.3f %10.3f %10.3f %10.3f\n', 'alpha',stats(2,:));
fprintf('%6s %10.3f %10.3f %10.3f %10.3f %10.3f\n', 'beta', stats(3,:));
fprintf('%6s %10.3f %10.3f %10.3f %10.3f %10.3f\n', 'gamma',stats(4,:));
fprintf('%6s %10.3f %10.3f %10.3f %10.3f %10.3f\n', 'loss', stats(5,:));

end

