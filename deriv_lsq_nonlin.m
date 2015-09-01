function dy_dx = deriv_lsq_nonlin(rfun, x, y_hat, m, n)
% deriv_lsq_nonlin estimates the derivative of
%   y(x) = argmin_y |r(x, y)|^2 .
%
% Parameters:
% [r, dr_dx, dr_dy, d2r_dxdy] = rfun(x, y)
%   x has size [m, 1].
%   y has size [n, 1].
%   r has size [p, 1].
%   dr_dx has size [p, m] and dr_dx(i, j) = dr(i) / dx(j).
%   dr_dy has size [p, n] and dr_dy(i, j) = dr(i) / dy(j).
%   d2r_dxdy has size [p, m, n] and dr2_dxdy(i, j, k) = d^2 r(i) / dx(j) dy(k).
% y_hat is an estimate of y(x).
%
% Returns:
% dy_dx has size (n, m).

[r, dr_dx, dr_dy, d2r_dxdy] = rfun(x, y_hat);
A = dr_dy;
b = -r;
[~, mul_dy_dA, mul_dy_db] = solve_rect(A, b);
% (p, m, n) -> (p, n, m)
dA_dx = permute(d2r_dxdy, [1, 3, 2]);
dy_dx = zeros(n, m);
for i = 1:m
  dy_dx(:,i) = dy_dx(:,i) + mul_dy_dA(dA_dx(:,:,i));
end
db_dx = -dr_dx;
dy_dx = dy_dx + mul_dy_db(db_dx);

end
