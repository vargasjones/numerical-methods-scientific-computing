%CONVERGENCE_STUDY Empirical convergence of four explicit ODE solvers.
%
% Model: radiative cooling-type nonlinear IVP
%   theta' = k_r (theta^4 - gamma),  theta(0) = 1200.
%
% A high-accuracy ODE113 solution is used as a numerical reference. It is
% not an analytic/exact solution.

clearvars; close all; clc;

this_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(this_dir, 'solvers'));

kr = -2.20e-12;
gamma = 81e8;
theta0 = 1200;
T = 60;
rhs = @(t,theta) kr*(theta.^4 - gamma);

N_values = [20, 40, 80, 160];
h_values = T ./ N_values;

reference_options = odeset('RelTol',1e-12,'AbsTol',1e-13);
[~, theta_reference] = ode113(rhs, [0,T], theta0, reference_options);
theta_T_reference = theta_reference(end);

methods = {@explicit_euler, @heun, @rk4_classical, @rk4_38};
method_names = {'Explicit Euler','Heun','Classical RK4','RK4 3/8'};
errors = zeros(numel(methods), numel(N_values));

for m = 1:numel(methods)
    for j = 1:numel(N_values)
        t = linspace(0, T, N_values(j)+1);
        y = methods{m}(t, theta0, rhs);
        errors(m,j) = abs(y(end) - theta_T_reference);
    end
end

estimated_orders = zeros(numel(methods),1);
for m = 1:numel(methods)
    fit_coefficients = polyfit(log(h_values), log(errors(m,:)), 1);
    estimated_orders(m) = fit_coefficients(1);
end

figure;
for m = 1:numel(methods)
    loglog(h_values, errors(m,:), '-o', 'LineWidth', 1.2, ...
        'DisplayName', sprintf('%s (p \\approx %.2f)', method_names{m}, estimated_orders(m)));
    hold on;
end
grid on;
xlabel('Step size h');
ylabel('Endpoint error |y_h(T)-y_{ref}(T)|');
title('Empirical convergence of explicit ODE methods');
legend('Location','northwest');

fprintf('Estimated convergence orders:\n');
for m = 1:numel(methods)
    fprintf('  %-16s %.3f\n', method_names{m}, estimated_orders(m));
end
