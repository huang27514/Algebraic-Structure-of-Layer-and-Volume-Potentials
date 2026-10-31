function plot_LaplaceDLP_QBX_K0_20
% Classical QBX (F11-style), K = 0..20
% 左三幅:  K = 6,12,20 的 log10|err| 场
% 右一幅:  K = 0..20 的 max|err| 收敛曲线
% 数据来自新的 QBX_*.mat (K=0..20)
%
% 修改版：
%   1) 左三幅误差图画成 square
%   2) 右侧 convergence 图也画成 square
%   3) colorbar 下移并避免遮挡 xlabel
%   4) 整体格式与 QB2X / poly 图统一
%   5) 左三幅加 box: 中心 0 - i/3, 边长 2/3 (x in [-1/3,1/3], y in [-2/3,0])
%   6) max|err| 由 logErr_allK 重新计算: (a) 曲线 (x,s(x)) 以下的所有点,
%      (b) box 内的点; 右图同时画 (a) 和 (b)

clear; close all; clc;

%% 1. 读入数据
tmp = load('dQBX_dataxy.mat');
evalPoints = tmp.Expression1;          %#ok<NASGU> % 10100 x 1 complex

tmp = load('dQBX_error.mat');
logErr_allK = tmp.Expression1;         % 21 x 10100, log10|err|

tmp = load('dQBX_Nconv.mat');
convMax = tmp.Expression1(:);          % 21 x 1, max|err| (整个 warp 区域)

tmp = load('dQBX_KList.mat');
Kvalues = tmp.Expression1(:).';        % 0:20
nK      = numel(Kvalues);              %#ok<NASGU>

%% 2. 重建 warp 网格 (和 F11/F12 完全一致)
nx = 101; ny = 100;

x_vec   = -1/3 : 1/150 : 1/3;
y_param = -1/150 : -1/150 : -2/3;

[yy_param, xx] = meshgrid(y_param, x_vec);
curve = 2 * xx.^2 + 5 * xx.^4;
yy    = yy_param .* (1 + 3/2*curve) + curve;

s  = @(x) 2 * x.^2 + 5 * x.^4;
xb = linspace(min(x_vec), max(x_vec), 3000);
yb = s(xb);

maskInside = yy < curve; % only check points inside the curve.

% box centered at 0 - i/3, length = width = 2/3
boxC = [0, -1/3];  boxL = 2/3;  tolB = 1e-12;
bx = boxC(1) + boxL/2*[-1 1 1 -1 -1];
by = boxC(2) + boxL/2*[-1 -1 1 1 -1];
maskBox = maskInside & abs(xx - boxC(1)) <= boxL/2 + tolB & ...
                       abs(yy - boxC(2)) <= boxL/2 + tolB;

%% 2b. max|err| for each K: (a) below the curve, (b) inside the box
errAll = 10.^logErr_allK;               % 21 x 10100; -Inf -> 0
errAll(isnan(errAll)) = 0;
mIn  = reshape(maskInside.', 1, []);    % same ordering as the rows of logErr_allK
mBox = reshape(maskBox.',    1, []);
maxBelow = max(errAll(:, mIn),  [], 2); % 21 x 1
maxBox   = max(errAll(:, mBox), [], 2); % 21 x 1
fprintf('points below curve: %d,  points in box: %d\n', nnz(mIn), nnz(mBox));
fprintf('  K    max err below curve    max err in box\n');
for kk = 1:numel(Kvalues)
    fprintf('%3d    %12.3e         %12.3e\n', Kvalues(kk), maxBelow(kk), maxBox(kk));
end
if any(abs(log10(maxBelow) - log10(convMax)) > 1e-6)
    warning('max error below the curve differs from dQBX_Nconv.mat (convMax).');
end

%% 3. 选择 K = 6,12,20 三个 snapshot
Kshow = [6 12 20];
Z = cell(3,1);

for idx = 1:3
    K = Kshow(idx);
    kk = find(Kvalues == K, 1);
    if isempty(kk)
        error('K = %d not found in Kvalues.', K);
    end

    row = logErr_allK(kk,:);
    row(~isfinite(row)) = -16;

    Zmat = reshape(row, [ny nx]).';
    Zmat(~maskInside) = NaN;
    Z{idx} = Zmat;
end

%% 4. 画图
fig = figure('Units','inch','Position',[0 0 13.5 5.0],'Color','w');
t   = tiledlayout(1,4,'TileSpacing','Compact','Padding','loose');

clims = [-16 0];
sub = gobjects(3,1);
colormap('jet');

% 左三幅
for i = 1:3
    sub(i) = nexttile(t,i);

    surf(xx, yy, Z{i}, 'EdgeColor','none');
    view(2);
    hold on;
    plot(xb, yb, 'k-', 'LineWidth', 2);
    plot3(bx, by, ones(size(bx)), 'b-.', 'LineWidth', 2);   % box (z = 1, above the surface)
    hold off;

    axis tight;
    axis square;

    title(sprintf('$K = %d$', Kshow(i)), 'Interpreter','latex');

    % max|err| below the curve and inside the box
    kk = find(Kvalues == Kshow(i), 1);
    subtitle(sprintf('$\\log_{10}\\|err\\|_\\infty$'),'Interpreter','latex');

    caxis(clims);
    set(sub(i), 'FontSize', 14, 'Box','on');
end

ylabel(sub(1), 'QBX Method', 'Interpreter','latex');

% 收敛曲线 K=0..20
p = nexttile(t,4);
semilogy(Kvalues, maxBelow, 'o-', 'LineWidth', 2, 'MarkerSize', 6); hold on;
semilogy(Kvalues, maxBox,   's-', 'LineWidth', 2, 'MarkerSize', 6); hold off;
grid on;
legend({'below $\Gamma$', 'box'}, 'Interpreter','latex', 'Location','southwest');
xlabel('$K$', 'Interpreter','latex');
ylabel('$\|err\|_{\infty}$', 'Interpreter','latex');
title('QBX Convergence', 'Interpreter','latex');
set(p, 'FontSize', 14, 'Box','on');

axis square;

% 底部长 colorbar
rr = colorbar(sub(end));
rr.Location = 'southoutside';
rr.Label.Interpreter = 'latex';
rr.Label.String      = '';
rr.FontSize          = 14;

drawnow;
pos = rr.Position;
pos(1) = 0.10;
pos(2) = 0.055;   % 避免挡住 xlabel，同时不裁掉 colorbar tick labels
pos(3) = 0.80;
rr.Position = pos;

exportgraphics(fig, 'LaplaceDLP_QBX_K0_20_Kshow_6_12_20.pdf', 'Resolution', 300);

end