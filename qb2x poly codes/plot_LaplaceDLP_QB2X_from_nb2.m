function plot_LaplaceDLP_QB2X_from_nb2
% 使用 nb2 导出的 zdataxy.mat / zerror.mat
% 左三幅：K = 6,8,12 的 log10|err| 场
% 右一幅：收敛曲线，K = 0:20
% 底部 colorbar 覆盖整个横向长度
%
% 修改版：
%   1) 左三幅误差图画成 square
%   2) 右侧 convergence 图也画成 square
%   3) colorbar 下移并避免遮挡 xlabel
%   4) 整体格式与 poly 图统一

clear; close all; clc;

%% 1. 读入 nb2 数据
tmp = load('zdataxy.mat');
evalPoints = tmp.Expression1;          %#ok<NASGU> % 10100x1 complex

tmp = load('zerror.mat');
logErr_allK = tmp.Expression1;         % 41x10100, 行: K=0..40, 值: log10|err|

Kvalues = 0:(size(logErr_allK,1)-1);   % 0..40
nK      = numel(Kvalues);

%% 2. 构造 warped 网格
nx = 101; ny = 100;

x  = -1/3 : 1/150 : 1/3;               % 101 点
yp = -1/150 : -1/150 : -2/3;           % 100 点 (parameter y')

[yy_param, xx] = meshgrid(yp, x);      % yy_param: 参数坐标 y'
curve = 2 * xx.^2 + 5 * xx.^4;   % sc(x)
yy    = yy_param .* (1 + 3/2*curve) + curve;  % 实际 warped y

% 边界曲线
s  = @(x) 2 * x.^2 + 5 * x.^4;
xb = linspace(min(x), max(x), 3000);
yb = s(xb);

% 物理内部区域 y < sc(x)
maskInside = yy < curve;
maskInsideFlat = reshape(maskInside.', [], 1);   % 10100x1 logical

%% 3. 准备 K = 6,8,12 的误差场
Kshow = [6 8 12];
Z = cell(3,1);

for idx = 1:3
    K = Kshow(idx);

    % 对应行是 K+1（因为第一行是 K=0）
    row = logErr_allK(K+1,:);          % 长度 10100，行优先
    row(~isfinite(row)) = -16;

    % 先 reshape 为 [ny,nx]，再转置 -> [nx,ny]
    Zmat = reshape(row, [ny nx]).';
    Zmat(~maskInside) = NaN;

    Z{idx} = Zmat;
end

%% 4. 计算收敛曲线 max|err| (只在内部区域)
convMax = zeros(nK,1);

for k = 1:nK
    row = logErr_allK(k,:).';
    row(~isfinite(row)) = -16;
    maxLog = max(row(maskInsideFlat));
    convMax(k) = 10.^maxLog;
end

% 收敛曲线：只画 K = 0..20
convMask    = (Kvalues <= 20);
K_conv      = Kvalues(convMask);
convMaxPlot = convMax(convMask);

%% 5. 画图
fig = figure('Units','inch','Position',[0 0 13.5 5.0],'Color','w');
t   = tiledlayout(1,4,'TileSpacing','Compact','Padding','loose');

clims = [-16 0];
colormap('jet');

sub = gobjects(3,1);

% ---- 左三幅误差图 ----
for i = 1:3
    sub(i) = nexttile(t,i);

    surf(xx, yy, Z{i}, 'EdgeColor','none');
    view(2);
    hold on;
    plot(xb, yb, 'k-', 'LineWidth', 2);
    hold off;

    axis tight;
    axis square;

    title(sprintf('$K = %d$', Kshow(i)), 'Interpreter','latex');

    % 对应 K 的最大误差
    maxErrK = max(convMax(Kvalues == Kshow(i)));
    subtitle(sprintf('$\\log_{10}(\\|err\\|_\\infty) = %.3f$', ...
        log10(maxErrK)), 'Interpreter','latex');

    caxis(clims);
    set(sub(i), 'FontSize', 14, 'Box','on');
end

ylabel(sub(1), 'QB2X Method (Laplace DLP)', 'Interpreter','latex');

% ---- 右侧收敛曲线：K = 0..20 ----
p = nexttile(t,4);
semilogy(K_conv, convMaxPlot, 'o-', 'LineWidth', 2, 'MarkerSize', 6);
grid on;
xlabel('$K$', 'Interpreter','latex');
ylabel('$\|err\|_{\infty}$', 'Interpreter','latex');
title('QB2X Convergence (Laplace DLP)', 'Interpreter','latex');
set(p, 'FontSize', 14, 'Box','on');

axis square;

%% 6. 底部 colorbar
rr = colorbar(sub(end));
rr.Location = 'southoutside';
rr.Label.Interpreter = 'latex';
rr.Label.String      = '';
rr.FontSize          = 14;

drawnow;

% 手动把 colorbar 拉成一条长条
pos = rr.Position;   % [left bottom width height]
pos(1) = 0.10;
pos(2) = 0.055;      % 下移，但不裁掉数字
pos(3) = 0.80;
rr.Position = pos;

end