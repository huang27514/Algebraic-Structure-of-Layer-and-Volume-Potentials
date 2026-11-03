% plot_vol_lowrank_errors.m
% Error of the volume potential V = vol_SLP,top (closed form) + vol_branch(a,b)
% + vol_low-rank (Taylor, K terms, single center w0 = -i/3) against the direct reference,
% on the layer-potential grid below Gamma. Reads the output of vol_closed_lowrank.wl.
%
% Errors are checked and plotted ONLY on the region
%   R = [-1/3, 1/3] x [-2/3, min(0, s(x))];
% nodes outside R are masked (NaN). The mask is recomputed here from X, Y and
% combined with inRegion from the MAT file when present.
% Figure 1: log10|error| maps on R for K = 6, 12, 20 (shared color scale), with
%           the boundary of R, the curve Gamma and the expansion center w0.
% Figure 2: max and median |error| over R versus K, together with the
%           asymptotic rate C rho^K, rho = max_R|w-w0| / dist(w0, nearest corner of D)
%           (for s = x^2/3 + x^3/8 - x^4/6 the nearest corner is -1 + i s(-1); rho ~ 0.44).
%           (maxErrK, medErrK are computed over R in vol_closed_lowrank.wl.)

clear; close all;
dataDir = 'vol_closed_out';
Gload = load(fullfile(dataDir, 'zdataxy_vol.mat'));
Eload = load(fullfile(dataDir, 'zerror_vol.mat'));
G=Gload.Expression1; E=Eload.Expression1;

X = G.X;  Y = G.Y;
sfun = @(x) x.^2/3 + x.^3/8 - x.^4/6;      % boundary curve s(x)
Klist   = double(E.Klist(:)).';
maxErrK = double(E.maxErrK(:)).';
medErrK = double(E.medErrK(:)).';
Kplot   = [6 12 20];
w0 = double(G.w0(:)).';  rho = sqrt(2)/3;

% ---------- region R and mask ----------
xReg = [-1/3 1/3];  yRegMin = -2/3;  tolR = 1e-12;
inR = (X >= xReg(1) - tolR) & (X <= xReg(2) + tolR) & (Y >= yRegMin - tolR) & ...
      (Y <= min(0, sfun(X)) + tolR);
boxC = [0, -1/3];  boxL = 2/3;                    % box centered at 0 - i/3
bxv = boxC(1) + boxL/2*[-1 1 1 -1 -1];
byv = boxC(2) + boxL/2*[-1 -1 1 1 -1];
if isfield(G, 'inRegion')
    M = logical(double(G.inRegion));
    if ~isequal(size(M), size(X)) && isequal(size(M.'), size(X)), M = M.'; end
    if isequal(size(M), size(X))
        if any(M(:) ~= inR(:))
            warning('inRegion from the MAT file differs from the mask computed here at %d nodes; using the intersection.', nnz(M ~= inR));
        end
        inR = inR & M;
    end
end
fprintf('region R: %d of %d nodes\n', nnz(inR), numel(inR));

% ---------- Figure 1: error maps ----------
% Use the 2D variables logErrK6/K12/K20 (always plain matrices in the MAT file)
% instead of slicing the 3D array logErrAll.
fprintf('size(logErrAll) = %s, numel(Klist) = %d, size(X) = %s\n', ...
    mat2str(size(E.logErrAll)), numel(Klist), mat2str(size(X)));
names = {'logErrK6', 'logErrK12', 'logErrK20'};
L = cell(1, numel(Kplot));
for p = 1:numel(Kplot)
    A = double(E.(names{p}));
    if ~isequal(size(A), size(X)) && isequal(size(A.'), size(X)), A = A.'; end
    A(~inR) = NaN;                                  % only the region R
    L{p} = A;
    fprintf('%s: size %s, finite values %d of %d\n', names{p}, mat2str(size(A)), ...
        nnz(isfinite(A)), nnz(inR));
end
allL = cell2mat(cellfun(@(A) A(:), L, 'UniformOutput', false));
allL = allL(isfinite(allL));
if isempty(allL)
    error('No finite error values in %s - check the Mathematica run (NaN/empty results).', dataDir);
end
cl = [floor(min(allL)), ceil(max(allL))];
if cl(2) <= cl(1), cl(2) = cl(1) + 1; end      % all values equal -> widen

fig1 = figure('Units', 'centimeters', 'Position', [2 2 30 10]);
tl = tiledlayout(fig1, 1, numel(Kplot)+1, 'TileSpacing', 'compact', 'Padding', 'compact');
xs = linspace(xReg(1), xReg(2), 400);
for p = 1:numel(Kplot)
    axpp(p) = nexttile(tl,p);
    ax=axpp(p);
    pcolor(ax, X, Y, L{p}); shading(ax, 'interp');
    colormap(ax, jet); caxis(ax, cl);   % (clim in R2022a+)
    hold(ax, 'on');
    plot(ax, xs, sfun(xs), 'k-', 'LineWidth', 1.2);              % Gamma
    plot(ax, [xReg(1) xReg(2) xReg(2) xReg(1) xReg(1)], ...      % boundary of R
         [yRegMin yRegMin 0 0 yRegMin], 'w-', 'LineWidth', 3.0);   % (min(0,s) = 0 here)
    plot(ax, bxv, byv, 'b-', 'LineWidth', 1.5);                   % box: center 0 - i/3, side 2/3
    plot(ax, w0(1), w0(2), 'kp', 'MarkerSize', 11, 'MarkerFaceColor', 'w');   % center
    axis(ax, 'equal'); xlim(ax, xReg); ylim(ax, [yRegMin, 0.05]); box(ax, 'on');
    xlabel(ax, 'x'); ylabel(ax, 'y');
    k = find(Klist == Kplot(p));
    title(ax, sprintf('K = %d,  max error = %.1e', Kplot(p), maxErrK(k)));
end

linkprop(axpp(1:3), 'CLim'); 
% Place the colorbar at the bottom ('south') of the left 3-plot group
cb = colorbar(axpp(3)); cb.Location = 'southoutside'; cb.Layout.Tile = 'south'; 
cb.Label.String = 'log_{10} |error|';
% rr = colorbar(sub(end));
% rr.Location = 'southoutside';
cb.Label.Interpreter = 'latex';
cb.Label.String      = '';
cb.FontSize          = 14;

axpp(4)=nexttile(tl,4);
ax=axpp(4);
semilogy(ax,Klist, maxErrK, 'o-', 'LineWidth', 1.2, 'MarkerSize', 5); hold on;
%semilogy(Klist, maxErrK(end) * rho.^(Klist - Klist(end)), 'k--', 'LineWidth', 1);  
%for K = Kplot
%    xline(K, ':', sprintf('K = %d', K), 'LabelVerticalAlignment', 'bottom');
%end
grid on; box on;axis square;
xlabel(ax,'$K$', 'Interpreter','latex');
ylabel(ax,'$\|err\|_{\infty}$', 'Interpreter','latex');
%legend({'max over R', 'median over R', sprintf('\\rho^K, \\rho = %.3f', rho)}, 'Location', 'southwest');
title(ax,'Convergence in K');




exportgraphics(fig1, 'vol_lowrank_error_maps.pdf', 'ContentType', 'vector');
