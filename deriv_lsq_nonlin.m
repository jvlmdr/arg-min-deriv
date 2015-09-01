% y(x) = argmin_y |r(x, y)|^2
% using approximation
% |r(x, y)|^2 ~ |r(x, y_hat+dy)|^2
%
% Returns:
% dy/dx(x)
%
% Parameters:
% [r, dr_dx, dr_dy, d2r_dxdy] = rfun(x, y)
%   d2r_dxdy(i, j, k) = d2r(i) / dx(j) dy(k).
% y_hat = argmin_y |r(x, y)|^2 (roughly)
% x has length m, y has length n.
function dy_dx = deriv_lsq_nonlin(rfun, x, y_hat, m, n)

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
