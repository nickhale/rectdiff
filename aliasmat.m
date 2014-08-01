function A = aliasmat(m, n)
%ALIASMAT  Aliasing matrix.

% Deal with the trivial case:
if ( nargin == 1 )
    n = m; m = n - 1;
    A = speye(m, n);
    return
end

A = [];
A0 = [speye(m), zeros(m,1), -flipud(speye(m, m-1))];
sgn = 1;
while ( size(A, 2) < n )
    A = [A, sgn*A0]; %#ok<AGROW>
    sgn = -sgn;
end
A = A(:,1:n);

end






