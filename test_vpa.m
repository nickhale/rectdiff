% VPA computation of rectangular diffmat. 

clc, clear all, close all

%#ok<*SAGROW>
load vpa_diffmats

NN = 10:100;
MM = NN/NN(1);
k = 0;

for N = NN
    
    disp(N)
    k = k + 1;
    
    % VPA
    D_vpa = D_vpa_store{N};   PD_vpa = PD_vpa_store{N};
    % Square
    D = diffmat(N);
    err = D_vpa - D;          err_sqr(k) = norm(err(:), inf);     % abs err
    % Explicit formula
    PD_exp = rectdiff_exp(N);
    err = PD_vpa - PD_exp;    err_exp(k) = norm(err, inf);        % abs err
    err = err./PD_vpa;        err_exp_rel(k) = norm(err, inf);    % rel err
    % Barycentric formula
    PD_bary = rectdiff_bary(N);
    err = PD_vpa - PD_bary;   err_bary(k) = norm(err, inf);       % abs err
    err = err./PD_vpa;        err_bary_rel(k) = norm(err, inf);   % rel err
    % Aliasing
    PD_alias = rectdiff_alias(N);
    err = PD_vpa - PD_alias;  err_alias(k) = norm(err, inf);      % abs err
    err = err./PD_vpa;        err_alias_rel(k) = norm(err, inf);  % rel err
    % Coefficients
    PD_coeff = rectdiff_coeff(N);
    err = PD_vpa - PD_coeff;  err_coeff(k) = norm(err, inf);      % abs err
    err = err./PD_vpa;        err_coeff_rel(k) = norm(err, inf);  % rel err
    
end

%% Absolute errors

figure(1)
LW = 'LineWidth'; lw = 2; FS = 'FontSize';
loglog(NN, err_exp, 'r', NN, err_bary, 'g', NN, err_alias, 'm', ...
    NN, err_coeff, 'c', NN, err_sqr, 'b', NN, err_alias, 'm',  LW, lw); hold on
c = 1e-14;
loglog(NN, c*MM.^2, '--k', NN, c*MM.^4, '--k', LW, lw); hold off
legend('explicit', 'barycentric', 'aliasing', 'coefficient','square',  ...
    'location', 'nw')
axis([10, NN(end), 1e-15, 1e-8])
set(gca, 'FontSize', 18), legend(gca, 'boxoff')

set(gcf, 'position', get(gcf, 'position')-[900 0 0 0])
print -depsc abs_err

%% Relative errors

figure(2)
loglog(NN, err_exp_rel, 'r', NN, err_bary, 'g', ...
    NN, err_alias_rel, 'm', NN, err_coeff_rel, 'c', LW, lw); hold on
c1 = 6e-15; c2 = 2e-14; c3 = 1e-13;
loglog(NN, c1*MM.^1.5, '--k', NN, c2*MM.^2.5, '--k', NN, c3*MM.^3, '--k', ...
    LW, lw); hold off
legend('explicit', 'barycentric', 'aliasing', 'coefficient', 'location', 'nw')
axis([10, NN(end), 1e-15, 1e-8])
set(gca, 'FontSize', 18), legend(gca, 'boxoff')

print -depsc rel_err
