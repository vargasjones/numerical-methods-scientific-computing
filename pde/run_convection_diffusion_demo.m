%RUN_CONVECTION_DIFFUSION_DEMO Diffusion-transport test with known solution.
%
% The test problem satisfies
%   -Delta u + alpha*u_x = 0,  (x,y) in (0,1)^2,
% with Dirichlet data taken from the exact one-dimensional profile.
%
% Two meshes illustrate the effect of the mesh Peclet number on a centered
% finite-difference discretization.

clearvars; close all; clc;

alpha = 40;
exact = @(x,y) (exp(alpha*(x-1))-1) ./ (exp(-alpha)-1) + 0*y;
source = @(x,y) 0*x + 0*y;
boundary = exact;

mesh_sizes = [10, 80];

figure;
for k = 1:numel(mesh_sizes)
    n = mesh_sizes(k);
    [x,y,u,rel_error,Pe] = convection_diffusion_2d( ...
        0,1,0,1,n,n,source,boundary,alpha,exact);

    subplot(1,2,k);
    surf(x,y,u,'EdgeColor','none');
    xlabel('x'); ylabel('y'); zlabel('u_h');
    title(sprintf('n=%d, Pe=%.2f, rel. error=%.2e', n, Pe, rel_error));
    view(35,30);
end
sgtitle('Centered finite differences for a convection-diffusion problem');
