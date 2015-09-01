function [f, dfdx] = square(x)
f = x .* x;
dfdx = diag(sparse(2 * x));
end
