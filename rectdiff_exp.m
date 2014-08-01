function D = rectdiff_exp(m, n)
%RECTDIFF_COEFF  Explicit constrcution rectangular differentiation matrix.

% Parse input:
if ( nargin == 1 )
    n = m;
    m = n - 1;
end

nm1 = n - 1;                    % For convenience.
cm1 = nm1 - m;                  % Difference between dimensions:
t = chebpts(n).';               % Second-kind grid.
tau = chebpts(m, 1);            % First-kind grid.
T = (nm1:-1:0)*pi/nm1;          % Second-kind grid (angles).
TAU = ((m-1:-1:0).'+.5)*pi/m;   % First-kind grid (angles).

% Explicit expression:
denom = 2*bsxfun( @(u,v) sin((v+u)/2) .* sin((v-u)/2), T, TAU );
numer = bsxfun( @times, 1 - tau*t, cos(cm1*TAU)./sin(TAU) );

sgn = (-1)^cm1;

if ( cm1 == 0 )
    D = numer ./ denom.^2 / nm1;
else
    D = repmat(sin(cm1*TAU), 1, n)./denom + numer./denom.^2 / nm1;
    D = sgn*D;
end
D(:,[1,n]) = .5*D(:,[1,n]);     % Scaling for first and last columns.

% Flipping trick:
ii = logical(rot90(tril(ones(m, n)), 2));
rot90D = rot90(D,2);
D(ii) = sgn*rot90D(ii);

% Sign:
D(1:2:end,1:2:end) = -D(1:2:end,1:2:end);
D(2:2:end,2:2:end) = -D(2:2:end,2:2:end);

% Negative sum trick:
[~, idx] = min(abs(denom), [], 2); 
idx = sub2ind([m n], 1:m, idx.');
D(idx) = 0; D(idx) = -sum(D, 2);

if ( cm1 == 0 )
    % Fix corner values:
    D(1) = -.25/(nm1*sin(pi/(2*m))*sin(pi/(4*m))^2); D(end) = -D(1);
    % Negative sum trick for corner entries:
    D(1,2) = -sum(D(1,[1 3:end])); D(end,end-1) = -D(1,2);
end

end