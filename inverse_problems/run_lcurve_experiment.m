%RUN_LCURVE_EXPERIMENT Tikhonov regularization for an ill-conditioned system.
%
% The experiment constructs A with exponentially decaying singular values,
% perturbs exact data by controlled relative Gaussian noise, and selects a
% regularization parameter from the maximum geometric curvature of the
% log-log L-curve.

clearvars; close all; clc;
rng(7);

n = 400;
m = 200;

C_left = test_matrix(n);
C_right = test_matrix(m);
[U,~,~] = svd(C_left,'econ');
[~,~,V] = svd(C_right,'econ');

singular_values = exp(-0.4*(0:m-1)).';
A = U(:,1:m) * diag(singular_values) * V';

x_true = ones(m,1);
b_exact = A*x_true;

noise_levels = [1e-4, 1e-3, 1e-2];
lambdas = logspace(-8,0,120);
selected_lambdas = zeros(size(noise_levels));

figure;
for i = 1:numel(noise_levels)
    eta = randn(n,1);
    eta = eta/norm(eta);
    b = b_exact + noise_levels(i)*norm(b_exact)*eta;

    X = tikhonov_path(A,b,lambdas);
    residual_norms = vecnorm(A*X-b,2,1).';
    solution_norms = vecnorm(X,2,1).';
    curvature = lcurve_curvature(lambdas,residual_norms,solution_norms);

    [~,corner_index] = max(curvature);
    selected_lambdas(i) = lambdas(corner_index);

    subplot(1,numel(noise_levels),i);
    loglog(residual_norms,solution_norms,'-'); hold on;
    loglog(residual_norms(corner_index),solution_norms(corner_index),'o', ...
        'MarkerSize',8,'LineWidth',1.5);
    xlabel('||Ax_\lambda-b||_2');
    ylabel('||x_\lambda||_2');
    title(sprintf('noise=%.0e, \\lambda*=%.2e', ...
        noise_levels(i),selected_lambdas(i)));
    grid on;
end
sgtitle('L-curve parameter selection for Tikhonov regularization');

fprintf('Selected regularization parameters:\n');
for i = 1:numel(noise_levels)
    fprintf('  relative noise %.0e -> lambda* = %.4e\n', ...
        noise_levels(i),selected_lambdas(i));
end

function C = test_matrix(N)
    j = (1:N).';
    k = 1:N;
    C = exp((pi*(2*j-1)/(4*N-2)) .* cos(pi*(2*k-1)/(2*N-1)));
end
