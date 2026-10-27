function plot_poly_direct_vs_taylor_deriv(dirDirect, dirTaylor, clims)
% PLOT_POLY_DIRECT_VS_TAYLOR_DERIV  Error distributions of the closed-form
% derivative phi'(w) = Int f/(P-w)^2 dx (double precision) for two
% Mathematica runs (both made with doDerivative = True):
%   left : direct root sum           (data in folder dirDirect)
%   right: pair-Taylor + Newton      (data in folder dirTaylor)
% Each folder must contain the files written by poly_closed_form*.wl:
%   zdataxy_closed.mat : X, Y (nx x ny), Nx, Ny
%   zerror_closed.mat  : logErrDMach (nx x ny, log10|err phi'|), maxErrDMach
%
% Usage:
%   plot_poly_direct_vs_taylor_deriv                         % default folder names
%   plot_poly_direct_vs_taylor_deriv('out_direct','out_taylor')
%   plot_poly_direct_vs_taylor_deriv('out_direct','out_taylor',[-16 -11])
% Default colour limits [-16 -11]: the direct sum for phi' reaches ~4e-12
% (log10 ~ -11.4) near w_c = -0.6981i, so the potential's [-16 -12] would clip.

if nargin < 1 || isempty(dirDirect), dirDirect = 'poly_closed_out_Sum'; end
if nargin < 2 || isempty(dirTaylor), dirTaylor = 'poly_closed_out_Taylor'; end
if nargin < 3 || isempty(clims),     clims     = [-16 -11]; end

[X1, Y1, L1, m1] = readRun(dirDirect);
[X2, Y2, L2, m2] = readRun(dirTaylor);
assert(isequal(size(X1), size(X2)) && max(abs(X1(:)-X2(:))) < 1e-12 && ...
       max(abs(Y1(:)-Y2(:))) < 1e-12, 'The two runs use different grids.');

sfun = @(x) 5*x.^4 + 2*x.^2;                       % boundary Gamma: y = s(x)
xb = linspace(min(X1(:)), max(X1(:)), 3000); yb = sfun(xb);
wc = -0.6981i;                                     % critical value below Gamma

fig = figure('Units','inch','Position',[0 0 10.5 5.0],'Color','w');
t = tiledlayout(fig, 1, 2, 'TileSpacing','compact', 'Padding','loose');
runs = {X1, Y1, L1, m1, 'Direct root sum'; ...
        X2, Y2, L2, m2, 'Taylor Expansion-based Method'};
ax = gobjects(2,1);
for k = 1:2
    [X, Y, L, mx, ttl] = runs{k,:};
    ax(k) = nexttile(t);
    surf(ax(k), X, Y, zeros(size(X)), L, 'EdgeColor','none'); view(ax(k), 2);
    hold(ax(k),'on');
    plot(ax(k), xb, yb, 'k-', 'LineWidth', 2);
    plot(ax(k), real(wc), imag(wc), 'kx', 'MarkerSize', 9, 'LineWidth', 1.5);
    hold(ax(k),'off');
    axis(ax(k), [min(X(:)) max(X(:)) min(Y(:))-0.04 max(yb)]); axis(ax(k), 'square');
    colormap(ax(k), jet); setclim(ax(k), clims);
    title(ax(k), [ttl, ', $\phi''$'], 'Interpreter','latex');
    subtitle(ax(k), sprintf('$\\|err\\|_\\infty = %.2e$', mx), 'Interpreter','latex');
    xlabel(ax(k), '$x$', 'Interpreter','latex');
    set(ax(k), 'FontSize', 13, 'Box','on', 'Layer','top');
end
ylabel(ax(1), '$y$', 'Interpreter','latex');

cb = colorbar(ax(2));                              % one colorbar, same limits
cb.Layout.Tile = 'east';
cb.Label.String = '$\log_{10}|err\,\phi''|$'; cb.Label.Interpreter = 'latex'; cb.FontSize = 13;

exportgraphics(fig, 'LaplaceDLP_poly_direct_vs_taylor_deriv.pdf', 'ContentType','vector');
end

% ---------------------------------------------------------------------
function [X, Y, L, mx] = readRun(folder)
fx = fullfile(folder, 'zdataxy_closed.mat');
fe = fullfile(folder, 'zerror_closed.mat');
assert(isfile(fx) && isfile(fe), 'Missing zdataxy_closed.mat or zerror_closed.mat in %s', folder);
dload = load(fx); eload = load(fe);
d=dload.Expression1; e=eload.Expression1;
assert(isfield(e,'logErrDMach') && isfield(e,'maxErrDMach') && e.maxErrDMach > 0, ...
       'No derivative data in %s (run the .wl script with doDerivative = True).', folder);
X = d.X; Y = d.Y; L = e.logErrDMach;
assert(isequal(size(X), [d.Nx d.Ny]) && isequal(size(L), size(X)), 'Size mismatch in %s', folder);
L(~isfinite(L)) = -16;
mx = e.maxErrDMach;
end

function setclim(ax, lims)
if exist('clim','file') || exist('clim','builtin'), clim(ax, lims); else, caxis(ax, lims); end %#ok<CAXIS>
end
