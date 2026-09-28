function y = rk4_classical(t, y0, f)
%RK4_CLASSICAL Solve y' = f(t,y) with the classical fourth-order RK method.

    t = t(:).';
    if numel(t) < 2 || any(diff(t) <= 0)
        error('t must contain at least two strictly increasing time nodes.');
    end

    y0 = y0(:);
    y = zeros(numel(y0), numel(t));
    y(:,1) = y0;

    for k = 1:numel(t)-1
        h = t(k+1) - t(k);
        tk = t(k);
        yk = y(:,k);

        k1 = f(tk, yk);
        k2 = f(tk + h/2, yk + h*k1/2);
        k3 = f(tk + h/2, yk + h*k2/2);
        k4 = f(tk + h, yk + h*k3);

        y(:,k+1) = yk + h*(k1 + 2*k2 + 2*k3 + k4)/6;
    end
end
