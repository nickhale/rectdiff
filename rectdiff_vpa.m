function [PD, D] = rectdiff_vpa(m, n)
%RECTDIFF_VPA  VPA computation of rectangular diffmat.

% parse input:
if ( nargin == 1 )
    n = m;
    m = n - 1;
end

% For convenience:
n = n - 1; m = m - 1;

% Cast to VPA:
digits = 30;
p = vpa('pi', digits);
N = vpa(n, digits);
M = vpa(m, digits);

% DIFFMAT (VPA):
x = -cos(p*(0:N)/N).'; 
c = ones(n+1, 1); c(2:2:end) = -1; c([1, end]) = 2*c([1, end]);
X = repmat(x, 1, n+1); dX = X - X.';
D  = (c*(1./c).') ./ (dX + eye(n+1));
D  = D - diag(sum(D, 2));

% BARYMAT (VPA):
y = -cos(p*(2*(0:M)+1)/(2*(M+1))).';  w = 1./c.';
Y = repmat(y, 1, n + 1); X = repmat(x, 1, n); W = repmat(w, n, 1);
P = W./(Y - X.');
C = repmat(sum(P, 2), 1, n+1);
P = P./C;

% Rectangular diffmat:
PD = P*D;

% Cast to double:
D = double(D);
PD = double(PD);
    
end
