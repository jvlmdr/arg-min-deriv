function [x, mul_dx_dA, mul_dx_db] = solve_square(A, b)
% solve_square returns x = A^-1 b and operators that compute products with
% derivatives with respect to A and b. It uses an LU decomposition of A.
%
% Parameters:
% A has size [n, n] and rank(A) = n.
% b has size [n, 1].
%
% Returns:
% x has size [n, 1].
% v = mul_dx_dA(U)
%   u has size [n, n].
%   v has size [n, 1].
% v = mul_dx_db(u)
%   u has size [n, 1].
%   v has size [n, 1].

[m, n] = size(A);
if m ~= n
  error('not square');
end

% Old version using explicit inverse.
% C = inv(A);
% x = C * b;
% dx_db = C;
% dx_da = -kron(x', C);
% mul_dx_db = @(v) dx_db * v;
% mul_dx_dA = @(V) dx_da * V(:);

[L, U] = lu(A);
% A x = b
% L U x = b
% x = U \ (L \ b)
x = U \ (L \ b);
mul_dx_db = @(v) U \ (L \ v);
mul_dx_dA = @(V) -U \ (L \ (V * x));

end
