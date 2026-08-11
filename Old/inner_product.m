function result = inner_product(f1, f2, t)
    % Compute the inner product of two functions over the time grid
    dt = t(2) - t(1);
    result = sum(f1 .* f2) * dt;
end
