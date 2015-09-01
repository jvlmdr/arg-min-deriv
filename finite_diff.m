function pass = finite_diff(fun, x, d, abs_tol, rel_tol)

[f_x, mul_dfdx_x] = fun(x);
% d = rand_unit(size(x));
% y = x + mag * d;
y = x + d;
mag = norm(d);
[f_y, ~] = fun(y);
g = (f_y - f_x) / mag;
h = 1/mag * mul_dfdx_x(d);
if ~equal(g, h, abs_tol, rel_tol)
  fprintf('FAIL:\n');
  fprintf('finite-diff: %s\n', mat2str(g, 4));
  fprintf('derivative:  %s\n', mat2str(h, 4));
  fprintf('abs norm of difference: %.4g\n', norm(g-h));
  fprintf('rel norm of difference: %.4g\n', norm(g-h)/max(norm(g), norm(h)));
  pass = false;
  return
end
pass = true;

end

function eq = equal(a, b, abs_tol, rel_tol)
eq = norm(a-b) <= abs_tol || ...
    (norm(a) == 0 && norm(b) == 0) || ...
    norm(a-b)/max(norm(a), norm(b)) <= rel_tol;
end
