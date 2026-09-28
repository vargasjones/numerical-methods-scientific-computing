%STABILITY_REGIONS Absolute-stability regions for the implemented methods.

clearvars; close all; clc;

[x,y] = meshgrid(linspace(-4,2,700), linspace(-4,4,700));
z = x + 1i*y;

R = {
    1 + z, ...
    1 + z + z.^2/2, ...
    1 + z + z.^2/2 + z.^3/6 + z.^4/24 ...
};

names = {'Explicit Euler','Heun','Fourth-order RK'};

figure;
for k = 1:numel(R)
    subplot(1,3,k);
    contourf(x, y, abs(R{k}), [0,1], 'LineStyle','none');
    hold on;
    contour(x, y, abs(R{k}), [1,1], 'k', 'LineWidth', 1.4);
    plot([min(x(:)),max(x(:))],[0,0],'k-');
    plot([0,0],[min(y(:)),max(y(:))],'k-');
    axis equal tight;
    xlabel('Re(z)'); ylabel('Im(z)');
    title(names{k});
end
sgtitle('Absolute-stability regions: |R(z)| \leq 1');
