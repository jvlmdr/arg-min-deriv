function pass = deriv_lsq_nonlin_test()

mag = 1e-6;

do_trials(...
    @() @(x) rosenbrock_y(x, 1, 100), ...
    @() 10*(2*rand()-1), ...
    @() mag * randn_unit([1, 1]), ...
    'rosenbrock');

m = 3;
n = 4;
p = 10;
do_trials(...
    @() sample_quadratic_form_minimiser(m, n, p), ...
    @() randn(m, 1), ...
    @() mag * randn_unit([m, 1]), ...
    'quadratic form minimiser');

end

function do_trials(sample_fun, sample_x, sample_dx, name)
abs_tol = 1e-9;
rel_tol = 1e-5;
num = 100;
for t = 1:num
  fprintf('trial %d/%d\n', t, num);
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

function fun = sample_quadratic_form_minimiser(m, n, p)
a = randn(p, m+n+1, m+n+1);
rfun = @(x, y) quadratic_form(a, x, y, m, n, p);
y_guess = randn(n, 1);
fun = @(x) rfun_y_minimiser(rfun, x, y_guess, m, n);
end

function [y, mul_dy_dx] = rfun_y_minimiser(rfun, x, y_guess, m, n)
rfun_y = @(y) rfun(x, y);
y = lsqnonlin(@(y) lsqnonlin_fun(rfun_y, y), y_guess);
dy_dx = deriv_lsq_nonlin(rfun, x, y, m, n);
mul_dy_dx = @(u) dy_dx * u;
end

function [f, df_dy] = lsqnonlin_fun(fun, y)
[f, ~, df_dy, ~] = fun(y);
end

function [r, dr_dx, dr_dy, d2r_dxdy] = rosenbrock(x, y, a, b)
sqrt_b = sqrt(b);
r = [a-x; sqrt_b*(y-x^2)];
dr_dx = [-1; -sqrt_b*2*x];
dr_dy = [0; sqrt_b];
d2r_dxdy = [0; 0];
end

function [y, mul_dy_dx] = rosenbrock_y(x, a, b)
y = x^2;
rfun = @(x, y) rosenbrock(x, y, a, b);
dy_dx = deriv_lsq_nonlin(rfun, x, y, 1, 1);
mul_dy_dx = @(u) dy_dx * u;
end
