% Test multiplication of rectangular diffmat for integrating exp(x).

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
    % BC and RHS:
    BC = eye(1, N);    
    rhs = [exp(-1) ; op(tau)];
    
    % Square:
    D = diffmat(N); D = [BC ; D(2:end,:)];
    err_sqr(k) = norm(D\[exp(-1); f(2:end)]-f, inf);    
    % Explicit:
    D = [BC ; rectdiff_exp(N)];     err_exp(k) = norm(D\rhs-f, inf);
    % Barycentric:
    D = [BC ; rectdiff_bary(N)];    err_bary(k) = norm(D\rhs-f, inf);    
    % Aliasing:
    D = [BC ; rectdiff_alias(N)];   err_alias(k) = norm(D\rhs-f, inf);
    % Coefficient:
    D = [BC ; rectdiff_coeff(N)];   err_coeff(k) = norm(D\rhs-f, inf);
    
end

%%

LW = 'LineWidth'; lw = 2;

loglog(NN, err_exp, '-r', NN, err_bary, '-g', NN, err_alias, '-m', ...
    NN, err_coeff, '-c', NN, err_sqr, '-b',LW, lw), hold on
loglog(NN, eps*NN, '--k', 'LineWidth', 2); hold off
legend('explicit', 'barycentric', 'aliasing', 'coefficient', 'square', ...
    'O(n)', 'location', 'nw');
axis([10, 1400, 1e-15, 1e-12])
set(gca, 'FontSize', 18), legend(gca, 'boxoff')

print -depsc test_invDf
