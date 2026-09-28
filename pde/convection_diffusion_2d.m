function [x, y, u, rel_error, Pe] = convection_diffusion_2d( ...
    xmin, xmax, ymin, ymax, nx, ny, source, boundary, alpha, exact_solution)
%CONVECTION_DIFFUSION_2D Central finite differences for
%
%       -Delta u + alpha * du/dx = f
%
% on a rectangular domain with Dirichlet boundary data.
%
% nx, ny denote the number of interior grid points in x and y.
% The diffusion coefficient is one. The mesh Peclet number is therefore
% Pe = |alpha|*hx/2.

    if nargin < 10
        exact_solution = [];
    end
    if nx < 1 || ny < 1
        error('nx and ny must be positive integers.');
    end

    hx = (xmax-xmin)/(nx+1);
    hy = (ymax-ymin)/(ny+1);
    x = linspace(xmin, xmax, nx+2);
    y = linspace(ymin, ymax, ny+2);

    n_unknowns = nx*ny;
    A = spalloc(n_unknowns, n_unknowns, 5*n_unknowns);
    rhs = zeros(n_unknowns,1);

    center = 2/hx^2 + 2/hy^2;
    left   = -1/hx^2 - alpha/(2*hx);
    right  = -1/hx^2 + alpha/(2*hx);
    down   = -1/hy^2;
    up     = -1/hy^2;

    index = @(i,j) i + (j-1)*nx;

    for j = 1:ny
        yj = y(j+1);
        for i = 1:nx
            xi = x(i+1);
            row = index(i,j);
            A(row,row) = center;
            rhs(row) = source(xi,yj);

            if i > 1
                A(row,index(i-1,j)) = left;
            else
                rhs(row) = rhs(row) - left*boundary(xmin,yj);
            end

            if i < nx
                A(row,index(i+1,j)) = right;
            else
                rhs(row) = rhs(row) - right*boundary(xmax,yj);
            end

            if j > 1
                A(row,index(i,j-1)) = down;
            else
                rhs(row) = rhs(row) - down*boundary(xi,ymin);
            end

            if j < ny
                A(row,index(i,j+1)) = up;
            else
                rhs(row) = rhs(row) - up*boundary(xi,ymax);
            end
        end
    end

    u_interior = A \ rhs;

    u = zeros(ny+2,nx+2);
    u(1,:)   = arrayfun(@(xx) boundary(xx,ymin), x);
    u(end,:) = arrayfun(@(xx) boundary(xx,ymax), x);
    u(:,1)   = arrayfun(@(yy) boundary(xmin,yy), y).';
    u(:,end) = arrayfun(@(yy) boundary(xmax,yy), y).';

    for j = 1:ny
        for i = 1:nx
            u(j+1,i+1) = u_interior(index(i,j));
        end
    end

    Pe = abs(alpha)*hx/2;
    if Pe > 1
        warning(['Mesh Peclet number Pe = %.3g > 1. Central differences ' ...
            'may produce non-physical oscillations in the transport-dominated regime.'], Pe);
    end

    rel_error = [];
    if ~isempty(exact_solution)
        [X,Y] = meshgrid(x,y);
        u_exact = arrayfun(exact_solution, X, Y);
        max_exact = max(abs(u_exact(:)));
        max_error = max(abs(u(:)-u_exact(:)));
        if max_exact > 0
            rel_error = max_error/max_exact;
        else
            rel_error = max_error;
        end
    end
end
