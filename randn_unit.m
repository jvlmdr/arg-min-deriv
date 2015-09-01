function x = randn_unit(sz)
x = randn(sz);
x = x / norm(x(:));
end
