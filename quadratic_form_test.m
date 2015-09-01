function pass = quadratic_form_test()

m = 3;
n = 4;
p = 10;
mag = 1e-9;

sample_x = @() randn(m, 1);
sample_dx = @() mag * randn_unit([m, 1]);
sample_y = @() randn(n, 1);
sample_dy = @() mag * randn_unit([n, 1]);

sample_fy = @() sample_fun_x(@(a, x, y) fy(a, x, y, m, n, p), m, n, p);
do_trials(sample_fy, sample_x, sample_dx, 'd/dx fy');

sample_fx = @() sample_fun_y(@(a, x, y) fx(a, x, y, m, n, p), m, n, p);
do_trials(sample_fx, sample_y, sample_dy, 'd/dy fx');

sample_dfy_dx = @() sample_fun_x(@(a, x, y) dfy_dx(a, x, y, m, n, p), m, n, p);
do_trials(sample_dfy_dx, sample_x, sample_dx, 'd/dx d/dy fx');

sample_dfx_dy = @() sample_fun_y(@(a, x, y) dfx_dy(a, x, y, m, n, p), m, n, p);
do_trials(sample_dfx_dy, sample_y, sample_dy, 'd/dy d/dx fy');

end

function do_trials(sample_fun, sample_x, sample_dx, name)
abs_tol = 1e-9;
rel_tol = 1e-5;
for t = 1:100
  fun = sample_fun();
  x = sample_x();
  dx = sample_dx();
  pass = finite_diff(fun, x, dx, abs_tol, rel_tol);
  if ~pass
    fprintf('FAIL: %s\n', name);
    return
  end
end
fprintf('PASS: %s\n', name);
end

function fun = sample_fun_x(g, m, n, p)
a = randn(p, m+n+1, m+n+1);
y = randn(n, 1);
fun = @(x) g(a, x, y);
end

function fun = sample_fun_y(g, m, n, p)
a = randn(p, m+n+1, m+n+1);
x = randn(m, 1);
fun = @(y) g(a, x, y);
end

function [f, mul_df_dx] = fy(a, x, y, m, n, p)
[f, df_dx, ~, ~] = quadratic_form(a, x, y, m, n, p);
mul_df_dx = @(u) df_dx * u;
end

function [f, mul_df_dy] = fx(a, x, y, m, n, p)
[f, ~, df_dy, ~] = quadratic_form(a, x, y, m, n, p);
mul_df_dy = @(u) df_dy * u;
end

function [g, mul_dg_dx] = dfy_dx(a, x, y, m, n, p)
[~, ~, g, dg_dx] = quadratic_form(a, x, y, m, n, p);
g = g(:);
% (p, m, n) -> ((p, n), m)
dg_dx = reshape(permute(dg_dx, [1, 3, 2]), [p*n, m]);
mul_dg_dx = @(u) dg_dx * u;
end

function [g, mul_dg_dy] = dfx_dy(a, x, y, m, n, p)
[~, g, ~, dg_dy] = quadratic_form(a, x, y, m, n, p);
g = g(:);
% (p, m, n) -> ((p, m), n)
dg_dy = reshape(dg_dy, [p*m, n]);
mul_dg_dy = @(u) dg_dy * u;
end
