function y = rk4_38(t, y0, f)
%RK4_38 Solve y' = f(t,y) with the fourth-order Runge-Kutta 3/8 rule.

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
        k2 = f(tk + h/3, yk + h*k1/3);
        k3 = f(tk + 2*h/3, yk - h*k1/3 + h*k2);
        k4 = f(tk + h, yk + h*k1 - h*k2 + h*k3);

        y(:,k+1) = yk + h*(k1 + 3*k2 + 3*k3 + k4)/8;
    end
end
