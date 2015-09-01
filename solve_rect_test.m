function pass = solve_rect_test()

m = 64;
n = 32;
mag = 1e-6;
abs_tol = 1e-9;
rel_tol = 1e-6;

for t = 1:100
  A = randn(m, n);
  b = randn(m, 1);
  x = args_to_vec(A, b);
  dA = randn(m, n);
  db = randn(m, 1);
  d = args_to_vec(dA, db);
  d = mag * d / norm(d);
  pass = finite_diff(@(x) solve_vec(x, m, n), x, d, abs_tol, rel_tol);
  if ~pass
    return
  end
end
fprintf('PASS\n');
pass = true;

end

function [f, df_dx] = solve_vec(x, m, n)
[A, b] = vec_to_args(x, m, n);
[f, df_dA, df_db] = solve_rect(A, b);
df_dx = @(u) mul_deriv_vec(df_dA, df_db, u, m, n);
end

% df_dA and df_db are functions that multiply by the matrix.
function g = mul_deriv_vec(df_dA, df_db, u, m, n)
[u_A, u_b] = vec_to_args(u, m, n);
g = df_dA(u_A) + df_db(u_b);
end

function x = args_to_vec(A, b)
x = [A(:); b];
end

function [A, b] = vec_to_args(x, m, n)
A = reshape(x(1 : m*n), [m, n]);
b = reshape(x(m*n + (1:m)), [m, 1]);
end
