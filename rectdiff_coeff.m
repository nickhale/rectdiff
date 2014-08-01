function PD = rectdiff_coeff(m, n)
%RECTDIFF_COEFF  Compute rectangular differentiation matrix in coefficient space

% Parse input:
if ( nargin == 1 )
    n = m;
    m = n - 1;
end

% US differentiation and conversion matrices:
dg = .5*ones(n - 1, 1);
S0 = spdiags( [ 1 0 ; .5 0 ; dg -dg ], [0, 2], n, n );
D = spdiags( (0:(n-1))', 1, n, n );

% Discrete Chebyshev transform and its inverse:
dct = @(u) chebtech1.coeffs2vals(flipud(u));
idct = @(u) flipud(chebtech2.vals2coeffs(u));

% Aliasing matrix:
A = aliasmat(m, n); 

% Rectangular differentiation matrix:
PD = dct(A*(S0\(D*idct(eye(n)))));

end