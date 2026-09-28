function y = explicit_euler(t, y0, f)
%EXPLICIT_EULER Solve y' = f(t,y) with the explicit Euler method.
%
%   y = explicit_euler(t, y0, f)
%
% Inputs
%   t  - strictly increasing vector of time nodes
%   y0 - scalar or column vector initial condition
%   f  - function handle f(t,y)
%
% Output
%   y  - solution values; each column corresponds to one time node

    t = t(:).';
    if numel(t) < 2 || any(diff(t) <= 0)
        error('t must contain at least two strictly increasing time nodes.');
    end

    y0 = y0(:);
    y = zeros(numel(y0), numel(t));
    y(:,1) = y0;

    for k = 1:numel(t)-1
        h = t(k+1) - t(k);
        y(:,k+1) = y(:,k) + h * f(t(k), y(:,k));
    end
end
