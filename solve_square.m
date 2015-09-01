% Returns x = A^-1 b with derivatives.
function [x, mul_dx_dA, mul_dx_db] = solve_square(A, b)

[m, n] = size(A);
if m ~= n
  error('assume square');
end

[L, U] = lu(A);
% A x = b
% L U x = b
% x = U \ (L \ b)
x = U \ (L \ b);

% dx/db = A^-1
% (dx/db) v = A^-1 v
mul_dx_db = @(v) U \ (L \ v);

% dx/dA = dx/dC dC/dA
% c = vec(C), a = vec(A)
% dx/dc = d/dc{C b} = d/dc{(b' kron I) c} = b' kron I
% dc/da = d{vec A^-1}/d{vec A} = - A^-T kron A^-1
% dx/da = dx/dc dc/da
%       = (b' kron I) (-A^-T kron A^-1)
%       = - (A^-1 b)' kron A^-1

% dx_da = -kron(x', C);
% dx_dA = reshape(dx_da, [n, n, n]);

% (dx/da) * v = - ((A^-1 b)' kron A^-1) vec(V)
%             = vec(-A^-1 V A^-1 b)
%             = -C V x
% U is n x n.
mul_dx_dA = @(V) -U \ (L \ (V * x));

end
