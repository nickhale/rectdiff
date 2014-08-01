function D = rectdiff_alias(m, n)
%RECTDIFF_ALIAS  Compute rectangular differentiation matrix via aliasing.

% Parse input:
if ( nargin == 1 )
    n = m;
    m = n - 1;
end

% Square matrix:
D = colloc2.diffmat(n);

% Convert to coefficient space:
coeffs = chebtech2.vals2coeffs(D);

% Apply aliasing:
% coeffs_aliased = chebtech1.alias(coeffs, m)
coeffs_aliased = rot90(aliasmat(m,n), 2)*coeffs;

% Revert to physical space:
D = chebtech1.coeffs2vals(coeffs_aliased);

end