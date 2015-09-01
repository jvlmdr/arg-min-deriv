function pass = solve_square_test()

n = 32;
mag = 1e-6;
abs_tol = 1e-9;
rel_tol = 1e-3;

for t = 1:100
  A = randn(n, n);
  b = randn(n, 1);
  x = args_to_vec(A, b);
  dA = randn(n, n);
  db = randn(n, 1);
  d = args_to_vec(dA, db);
  d = mag * d / norm(d);
  pass = finite_diff(@(x) solve_vec(x, n), x, d, abs_tol, rel_tol);
  if ~pass
    return
  end
end
fprintf('PASS\n');
pass = true;

end

function [f, df_dx] = solve_vec(x, n)
[A, b] = vec_to_args(x, n);
[f, df_dA, df_db] = solve_square(A, b);
df_dx = @(u) mul_deriv_vec(df_dA, df_db, u, n);
end

% df_dA and df_db are functions that multiply by the matrix.
function g = mul_deriv_vec(df_dA, df_db, u, n)
[u_A, u_b] = vec_to_args(u, n);
g = df_dA(u_A) + df_db(u_b);
end

function x = args_to_vec(A, b)
x = [A(:); b];
end

function [A, b] = vec_to_args(x, n)
A = reshape(x(1 : n*n), [n, n]);
b = reshape(x(n*n + (1:n)), [n, 1]);
end
