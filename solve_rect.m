function [x, mul_dx_dA, mul_dx_db] = solve_rect(A, b)
% solve_rect returns x that minimises |A x - b|^2
%   x = (A' A)^-1 A' b
% and operators that compute derivatives with respect to A and b.
% It uses a QR decomposition of A.
%
% Parameters:
% A has size [m, n] with m >= n and rank(A) = n.
% b has size [m, 1].
%
% Returns:
% x has size [n, 1].
% v = mul_dx_dA(U)
%   u has size [m, n].
%   v has size [n, 1].
% v = mul_dx_db(u)
%   u has size [m, 1].
%   v has size [n, 1].

% Old method using explicit inverse.
% G = A' * A;
% C = inv(G);
% x = C * (A' * b);
% mul_dx_db = @(v) C * (A' * v);
% mul_dx_dA = @(V) C * (V'*(b-A*x) - A'*(V*x));

[Q, R] = qr(A, 0);
% A = Q R, Q' Q = I (but not Q Q' = I)
% A' A x = A' b
% R' R x = R' Q' b
% R x = Q' b
x = R \ (Q' * b);
mul_dx_dA = @(V) R \ (R'\(V'*(b-A*x)) - Q'*V*x);
mul_dx_db = @(v) R \ (Q' * v);

end
