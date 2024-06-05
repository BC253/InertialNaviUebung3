function y = RungeKutta3(func, y0, x0, t)
    k1 = func(y0, x0(:, 1), x0(:, 2)); 
    k2 = func(y0 + t/2 * k1, x0(:, 1), x0(:, 2));
    k3 = func(y0 - t * k1 + 2 * t * k2, x0(:, 1), x0(:, 2));
    y = y0 + t/6 * (k1 + 4 * k2 + k3);
    y(7:10) = compact(normalize(quaternion(y(7:10)')));  
end

