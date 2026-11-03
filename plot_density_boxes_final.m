% plot_density_boxes.m
% Density f(x,y) = 1 + cos(2 pi x) cos(pi y) on the domain
%   {(x,y): -1 <= x <= 1, -4/3 <= y <= min(s(x), 2/3)},  s(x) = x^2/3 + x^3/8 - x^4/6,
% with the boundary curve Gamma = (x, s(x)) and the 3 x 3 boxes obtained by
% dividing [-1,1] x [-4/3, 2/3] (box side 2/3). The box centered at 0 - i/3,
% i.e. [-1/3,1/3] x [-2/3,0], is drawn in red; all other box edges in blue.

clear; close all;

ffun = @(x,y) 1 + cos(2*pi*x).*cos(pi*y);
sfun = @(x) x.^2/3 + x.^3/8 - x.^4/6;
yb   = -4/3;                 % bottom of D
yTop = 2/3;                  % plot f only up to min(s(x), yTop)

% ---------- curve-following grid on D ----------
nx = 401; ny = 400;
x  = linspace(-1, 1, nx);
t  = linspace(0, 1, ny);
[T, XX] = meshgrid(t, x);                    % nx x ny
YY = yb + (min(sfun(XX), yTop) - yb) .* T;   % y from -4/3 to min(s(x), 2/3)
FF = ffun(XX, YY);

% ---------- 3 x 3 boxes ----------
xE = [-1, -1/3, 1/3, 1];                     % x edges
yE = [-4/3, -2/3, 0, 2/3];                   % y edges
cBox = [-1/3 1/3; -2/3 0];                   % red box: [x1 x2; y1 y2], center 0 - i/3

% ---------- plot ----------
fig = figure('Units', 'centimeters', 'Position', [2 2 16 12], 'Color', 'w');
ax = axes(fig);
pcolor(ax, XX, YY, FF); shading(ax, 'interp');
colormap(ax, jet); caxis(ax, [0 2]);   % f ranges over [0, 2]
hold(ax, 'on');

xs = linspace(-1, 1, 2000);  xs = xs(sfun(xs) <= yTop + 0.05);   % part of Gamma in view
plot(ax, xs, sfun(xs), 'k-', 'LineWidth', 1.5);                         % Gamma

for xe = xE                                                             % blue grid
    plot(ax, [xe xe], [yE(1) yE(end)], 'b-', 'LineWidth', 1.2);
end
for ye = yE
    plot(ax, [xE(1) xE(end)], [ye ye], 'b-', 'LineWidth', 1.2);
end
plot(ax, cBox(1,[1 2 2 1 1]), cBox(2,[1 1 2 2 1]), 'r-', 'LineWidth', 2.2); % red box (on top)
plot(ax, 0, -1/3, 'r+', 'MarkerSize', 8, 'LineWidth', 1.2);              % its center

axis(ax, 'equal');
xlim(ax, [-1.05 1.05]); ylim(ax, [yb - 0.05, yTop + 0.05]);
box(ax, 'on'); set(ax, 'Layer', 'top', 'FontSize', 12);
xlabel(ax, '$x$', 'Interpreter', 'latex');
ylabel(ax, '$y$', 'Interpreter', 'latex');
title(ax, '$f(x,y) = 1 + \cos(2\pi x)\cos(\pi y)$', 'Interpreter', 'latex');
cb = colorbar(ax); cb.Label.String = 'f(x,y)';

exportgraphics(fig, 'density_boxes.pdf', 'ContentType', 'vector');
