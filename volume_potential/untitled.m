function plot_VolumePotential_GlobalVerticalQBX_Ks_30x30
% Plot convergence over K for:
%   vol_globalVertical_recursiveQBX_30x30_Ks_beta125.mat

clear; close all; clc;

fname = 'vol_globalVertical_recursiveQBX_30x30_Ks_beta125.mat';

D0 = load(fname);
if isfield(D0,'Expression1')
    D = D0.Expression1;
else
    D = D0;
end

Nx = double(D.Nx);
Ny = double(D.Ny);

KList = double(D.KList(:));
x = double(D.x(:));
y = double(D.y(:));
ratio = double(D.ratio(:));

directReal = double(D.directReal(:));
realErrMat = double(D.realErrMat);
complexErrMat = double(D.complexErrMat);
approxRealMat = double(D.approxRealMat);

sentinel = double(D.sentinel);

bad = realErrMat == sentinel | ~isfinite(realErrMat);
realErrMat(bad) = NaN;

badC = complexErrMat == sentinel | ~isfinite(complexErrMat);
complexErrMat(badC) = NaN;

maxRealErr = max(realErrMat, [], 2, 'omitnan');
meanRealErr = mean(realErrMat, 2, 'omitnan');
medianRealErr = median(realErrMat, 2, 'omitnan');

maxComplexErr = max(complexErrMat, [], 2, 'omitnan');
meanComplexErr = mean(complexErrMat, 2, 'omitnan');
medianComplexErr = median(complexErrMat, 2, 'omitnan');

fprintf('====================================\n');
fprintf('Global vertical recursive QBX: K sweep\n');
fprintf('Grid size: %d x %d\n', Nx, Ny);
fprintf('Number of points: %d\n', numel(x));
fprintf('KList:\n');
disp(KList.');

fprintf('K    maxRealErr        meanRealErr       medianRealErr\n');
for j = 1:length(KList)
    fprintf('%2d   %.6e   %.6e   %.6e\n', ...
        KList(j), maxRealErr(j), meanRealErr(j), medianRealErr(j));
end
fprintf('====================================\n');

outDir = 'figures_GlobalVerticalQBX_Ks_30x30';
if ~exist(outDir,'dir')
    mkdir(outDir);
end

%% Error vs K

figure('Units','inch','Position',[0 0 6.2 4.8],'Color','w');
semilogy(KList, maxRealErr, 'o-', 'LineWidth', 1.8); hold on;
semilogy(KList, meanRealErr, 's-', 'LineWidth', 1.8);
semilogy(KList, medianRealErr, '^-', 'LineWidth', 1.8);
grid on;
xlabel('$K_{\rm local}$','Interpreter','latex');
ylabel('real error','Interpreter','latex');
legend({'max','mean','median'}, 'Location','southwest');
title('Global Vertical Recursive QBX: Error vs $K$', 'Interpreter','latex');
set(gca,'FontSize',14,'Box','on');
exportgraphics(gcf, fullfile(outDir,'error_vs_K_real.pdf'), 'Resolution',300);

figure('Units','inch','Position',[0 0 6.2 4.8],'Color','w');
plot(KList, log10(maxRealErr), 'o-', 'LineWidth', 1.8); hold on;
plot(KList, log10(meanRealErr), 's-', 'LineWidth', 1.8);
plot(KList, log10(medianRealErr), '^-', 'LineWidth', 1.8);
grid on;
xlabel('$K_{\rm local}$','Interpreter','latex');
ylabel('$\log_{10}$ real error','Interpreter','latex');
legend({'max','mean','median'}, 'Location','southwest');
title('Global Vertical Recursive QBX: $\log_{10}$ Error vs $K$', 'Interpreter','latex');
set(gca,'FontSize',14,'Box','on');
exportgraphics(gcf, fullfile(outDir,'log_error_vs_K_real.pdf'), 'Resolution',300);

%% Area plot for selected K values

s = @(x) 5*x.^4 + 2*x.^2;
xb = linspace(min(x), max(x), 3000);
yb = s(xb);

xx = reshape(x,[Ny,Nx]).';
yy = reshape(y,[Ny,Nx]).';

xCols = xx(:,1);
yBoundary = s(xCols);

Xarea = [xCols, xx];
Yarea = [yBoundary, yy];

Kshow = [KList(1), KList(round(end/2)), KList(end)];

for a = 1:length(Kshow)
    K = Kshow(a);
    idx = find(KList == K, 1);

    Err = reshape(realErrMat(idx,:), [Ny,Nx]).';
    LogErr = log10(Err + eps);
    LogErr(~isfinite(LogErr)) = -16;

    LogErrArea = [LogErr(:,1), LogErr];

    figure('Units','inch','Position',[0 0 5.4 4.8],'Color','w');
    surf(Xarea, Yarea, LogErrArea, 'EdgeColor','none');
    view(2);
    hold on;
    plot(xb, yb, 'k-', 'LineWidth', 2);
    hold off;
    axis tight; axis square;
    colormap('jet');
    caxis([-16 0]);
    colorbar('southoutside');
    title(sprintf('$\\log_{10}$ error, $K=%d$', K), 'Interpreter','latex');
    ylabel('Global Vertical Recursive QBX', 'Interpreter','latex');
    set(gca,'FontSize',14,'Box','on');

    exportgraphics(gcf, fullfile(outDir, ...
        sprintf('area_log_error_K%d.pdf', K)), 'Resolution',300);
end

%% Approx/direct area plot for largest K

idx = length(KList);
approxReal = approxRealMat(idx,:).';

Direct = reshape(directReal,[Ny,Nx]).';
Approx = reshape(approxReal,[Ny,Nx]).';

DirectArea = [Direct(:,1), Direct];
ApproxArea = [Approx(:,1), Approx];

figure('Units','inch','Position',[0 0 10.8 4.8],'Color','w');
t = tiledlayout(1,2,'TileSpacing','Compact','Padding','loose');

ax1 = nexttile(t,1);
surf(Xarea,Yarea,DirectArea,'EdgeColor','none');
view(2); hold on; plot(xb,yb,'k-','LineWidth',2); hold off;
axis tight; axis square;
title('Direct','Interpreter','latex');
set(ax1,'FontSize',14,'Box','on');
colorbar(ax1,'southoutside');

ax2 = nexttile(t,2);
surf(Xarea,Yarea,ApproxArea,'EdgeColor','none');
view(2); hold on; plot(xb,yb,'k-','LineWidth',2); hold off;
axis tight; axis square;
title(sprintf('Recursive QBX, $K=%d$', KList(end)), 'Interpreter','latex');
set(ax2,'FontSize',14,'Box','on');
colorbar(ax2,'southoutside');

exportgraphics(gcf, fullfile(outDir,'direct_vs_approx_area_largestK.pdf'), 'Resolution',300);

fprintf('\nFigures saved to folder:\n%s\n', outDir);

end