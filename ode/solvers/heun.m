function y = heun(t, y0, f)
%HEUN Solve y' = f(t,y) with Heun's explicit second-order method.

    t = t(:).';
    if numel(t) < 2 || any(diff(t) <= 0)
        error('t must contain at least two strictly increasing time nodes.');
    end

    y0 = y0(:);
    y = zeros(numel(y0), numel(t));
    y(:,1) = y0;

    for k = 1:numel(t)-1
        h = t(k+1) - t(k);
        k1 = f(t(k), y(:,k));
        k2 = f(t(k+1), y(:,k) + h*k1);
        y(:,k+1) = y(:,k) + 0.5*h*(k1 + k2);
    end
end
