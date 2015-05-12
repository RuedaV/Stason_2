function [] = myDisplay_2(x1, x2, x3, x4, x5)

[m_omega, sd_omega, sk_omega, k_omega] = statistics(x1);
[m_alpha, sd_alpha, sk_alpha, k_alpha] = statistics(x2);
[m_beta, sd_beta, sk_beta, k_beta] = statistics(x3);
[m_gamma, sd_gamma, sk_gamma, k_gamma] = statistics(x4);
[m_QLIKE, sd_QLIKE, sk_QLIKE, k_QLIKE] = statistics(x5);

fprintf('%50s\n','-----------------------------------------------------');
fprintf('%9s %10s %10s %10s %10s\n', 'parameter', 'mean', 'sd', 'skewness', 'kurtosis');
fprintf('%50s\n','-----------------------------------------------------');
fprintf('%9s %10.6f %10.6f %10.6f %10.6f\n', 'omega', m_omega, sd_omega, sk_omega, k_omega);
fprintf('%9s %10.6f %10.6f %10.6f %10.6f\n', 'alpha', m_alpha, sd_alpha, sk_alpha, k_alpha);
fprintf('%9s %10.6f %10.6f %10.6f %10.6f\n', 'beta',  m_beta, sd_beta, sk_beta, k_beta);
fprintf('%9s %10.6f %10.6f %10.6f %10.6f\n', 'gamma', m_gamma, sd_gamma, sk_gamma, k_gamma);
fprintf('%9s %10.6f %10.6f %10.6f %10.6f\n', 'QLIKE', m_QLIKE, sd_QLIKE, sk_QLIKE, k_QLIKE);
fprintf('%50s\n','-----------------------------------------------------');

end

