% Solves A x = b where A is mxn with m >= n and rank(A) = n.
% A' A is nxn and is full rank.
% Returns x = (A' A)^-1 A' b with derivatives.
% Computes a QR decomposition of A.
function [x, mul_dx_dA, mul_dx_db] = solve_rect(A, b)

% G = A' * A;
% C = inv(G);
% x = C * (A' * b);
% mul_dx_dA = @(V) C * (V'*(b-A*x) - A'*(V*x));
% mul_dx_db = @(v) C * (A' * v);

% A = Q R, Q' Q = I (but not Q Q')
[Q, R] = qr(A, 0);
% A' A x = A' b
% R' R x = R' Q' b
% R x = Q' b
x = R \ (Q' * b);
mul_dx_dA = @(V) R \ (R'\(V'*(b-A*x)) - Q'*V*x);
mul_dx_db = @(v) R \ (Q' * v);

end
