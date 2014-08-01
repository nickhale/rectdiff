% Test multiplication of rectangular diff mat against exp(x).

clc, clear all, close all

%#ok<*SAGROW>
NN = round(2.^(4:0.125:10));
op = @(x) exp(x);

k = 0;
for N = NN
    
    disp(N)
    k = k + 1;

    % Chebyshev points:
    t = chebpts(N, 2); tau = chebpts(N - 1, 1);
    % Evaluate f and f':
    f = op(t); fp = op(tau);
    
    % Square:
    D = diffmat(N);         err_sqr(k) = norm(D*f - f, inf);
    % Explicit:
    D = rectdiff_exp(N);    err_exp(k) = norm(D*f - fp, inf);
    % Barycentric:
    D = rectdiff_bary(N);   err_bary(k) = norm(D*f - fp, inf);
    % Aliasing:
    D = rectdiff_alias(N);  err_alias(k) = norm(D*f - fp, inf);
    % Coefficient
    D = rectdiff_coeff(N);  err_coeff(k) = norm(D*f - fp, inf);
    
end

%%

figure(1)
LW = 'LineWidth'; lw = 2;
loglog(NN, err_exp, '-r', NN, err_bary, '-g', NN, err_alias, '-m', ..., 
    NN, err_coeff, '-c', NN, err_sqr, '-b', LW, lw), hold on
loglog(NN, eps*NN.^2, '--k', 'LineWidth', 2); hold off
legend('explicit', 'barycentric', 'aliasing', 'coefficients', 'square', ...
    'O(n^2)', 'location', 'nw');
axis([10 1400 1e-14 1e-9])
set(gca, 'FontSize', 18), legend(gca, 'boxoff')

print -depsc test_Df