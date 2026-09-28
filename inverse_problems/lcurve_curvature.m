function curvature = lcurve_curvature(lambdas, residual_norms, solution_norms)
%LCURVE_CURVATURE Estimate geometric curvature of the log-log L-curve.
%
% The curve is parameterized by tau = log(lambda):
%   x(tau) = log ||Ax_lambda-b||_2,
%   y(tau) = log ||x_lambda||_2.

    tau = log(lambdas(:));
    x = log(residual_norms(:));
    y = log(solution_norms(:));

    dx = gradient(x,tau);
    dy = gradient(y,tau);
    d2x = gradient(dx,tau);
    d2y = gradient(dy,tau);

    denominator = (dx.^2 + dy.^2).^(3/2);
    curvature = abs(dx.*d2y - dy.*d2x) ./ max(denominator, eps);

    % Endpoint derivatives are less reliable for a finite sampled curve.
    curvature([1,end]) = 0;
end
