function [r, dr_dx, dr_dy, d2r_dxdy] = quadratic_form(a, x, y, m, n, p)
% size(x) is [m, 1]
% size(y) is [n, 1]
% size(a) is [p, m+n+1, m+n+1]

z = [x; y; 1];

i_x = 1:m;
i_y = m+(1:n);

% r = reshape(a, [p, (m+n+1)^2]) * kron(z, z);
% 
% % d/dz kron(z, z)
% % = d/dz vec(z z')
% % = (z kron I) dz/dz + (I kron z) dz'/dz
% dr_dz = reshape(a, [p, (m+n+1)^2]) * ...
%   (kron(z, eye(m+n+1)) + kron(eye(m+n+1), z));

u = reshape(reshape(a, [p*(m+n+1), m+n+1]) * z, [p, m+n+1]);
r = u * z;

at = permute(a, [1, 3, 2]);
ut = reshape(reshape(at, [p*(m+n+1), m+n+1]) * z, [p, m+n+1]);

% d/dz U(z) z = (z' kron I_p) du/dz + (I_1 kron U(z)) dz/dz
%             = (z' kron I_p) du/dz + U(z)
dr_dz = u + ut;

dr_dx = dr_dz(:, i_x);
dr_dy = dr_dz(:, i_y);
d2r_dxdy = a(:, i_x, i_y) + at(:, i_x, i_y);

end
