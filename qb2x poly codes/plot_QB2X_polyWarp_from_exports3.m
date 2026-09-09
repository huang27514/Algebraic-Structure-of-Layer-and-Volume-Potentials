function plot_QB2X_polyWarp_from_exports3()
% Robust plotter for QB2X exports with suffix=2.
% Prefer MAT; fallback to CSV with robust numeric parsing.
%
% 左三幅：n = 6,8,12 的 log10|err| 场（带 boundary）
% 右一幅：收敛曲线
% 底部 colorbar 覆盖整个横向长度
%
% 修改：
%   1) 左三幅误差图可设为正方形
%   2) 右侧 convergence 图也设为正方形
%   3) colorbar 再往下放，避免挡住 xlabel
%   4) figure 高度增大，避免 colorbar tick labels 被裁掉

clear; close all; clc;

% ---------------- user options ----------------
nShow           = [6 8 12];   % must be in nList
errFloor        = 1e-16;      % avoid -Inf in log10
clims           = [-16 0];    % [] for auto, or e.g. [-16 0]
bndLW           = 2.0;        % boundary linewidth
useJet          = true;
useSquarePanels = true;       % true: square panels; false: equal aspect
useSquareConv   = true;       % true: square convergence panel
suffix          = '2';
% ---------------------------------------------

f_data_mat = ['zdataxy' suffix '.mat'];
f_err_mat  = ['zerror'  suffix '.mat'];
f_con_mat  = ['zNconvergence' suffix '.mat'];

f_data_csv = ['zdataxy' suffix '.csv'];
f_err_csv  = ['zerror'  suffix '.csv'];
f_con_csv  = ['zNconvergence' suffix '.csv'];

fprintf('[INFO] Using files with suffix="%s"\n', suffix);

% ---------------- load error data ----------------
[logErrAll, nList, nN, nP] = load_zerror_any(f_err_mat, f_err_csv);

fprintf('[zerror] size(logErr) = %dx%d\n', size(logErrAll,1), size(logErrAll,2));
fprintf('[zerror] finite entries = %d\n', sum(isfinite(logErrAll(:))));

% ---------------- load grid points ----------------
[Nx, Ny, xCoord, yCoord] = load_zdataxy_any(f_data_mat, f_data_csv, nP);
fprintf('[zdataxy] Nx=%d Ny=%d nP=%d\n', Nx, Ny, numel(xCoord));

% reshape to grid (flatten order consistent with Mathematica construction)
X = reshape(xCoord, [Ny, Nx]).';
Y = reshape(yCoord, [Ny, Nx]).';

% boundary
sc = @(x) 5*x.^4 + 2*x.^2;
xb = linspace(min(X(:)), max(X(:)), 3000);
yb = sc(xb);

% inside mask
maskInside = (Y < sc(X));
maskInsideFlat = reshape(maskInside.', [], 1);

% ---------------- convergence (optional) ----------------
[convMax, convOk] = load_conv_any(f_con_mat, f_con_csv, nList);

if ~convOk
    convMax = nan(nN,1);
    for ii = 1:nN
        row = logErrAll(ii,:).';
        row(~isfinite(row)) = log10(errFloor);
        row = max(row, log10(errFloor));
        maxLog = max(row(maskInsideFlat));
        convMax(ii) = 10.^maxLog;
    end
end

% ---------------- prepare Z for shown n ----------------
Z = cell(numel(nShow),1);
maxErrShown = nan(numel(nShow),1);

for k = 1:numel(nShow)
    n = nShow(k);
    idx = find(abs(nList - n) < 1e-12, 1);
    if isempty(idx)
        error('nShow contains n=%g not in nList.', n);
    end

    row = logErrAll(idx,:);
    row(~isfinite(row)) = log10(errFloor);
    row = max(row, log10(errFloor));

    Zmat = reshape(row, [Ny Nx]).';
    Zmat(~maskInside) = NaN;
    Z{k} = Zmat;

    if idx <= numel(convMax) && isfinite(convMax(idx))
        maxErrShown(k) = convMax(idx);
    else
        tmp = row(:);
        maxErrShown(k) = 10.^max(tmp(maskInsideFlat));
    end
end

% ---------------- figure/layout ----------------
fig = figure('Units','inch','Position',[0 0 13.5 5.0],'Color','w');
t = tiledlayout(1,4,'TileSpacing','Compact','Padding','loose');

if useJet
    colormap('jet');
end

% ---------------- left three panels ----------------
ax = gobjects(3,1);
for k = 1:numel(nShow)
    ax(k) = nexttile(t,k);

    surf(X, Y, Z{k}, 'EdgeColor','none');
    view(2);
    hold on;
    plot(xb, yb, 'k-', 'LineWidth', bndLW);
    hold off;

    if useSquarePanels
        axis tight;
        axis square;
    else
        axis equal tight;
    end

    title(sprintf('$n = %d$', nShow(k)), 'Interpreter','latex');
    subtitle(sprintf('$\\log_{10}(\\|err\\|_\\infty) = %.3f$', log10(maxErrShown(k))), ...
        'Interpreter','latex');

    if ~isempty(clims)
        caxis(clims);
    end

    set(ax(k), 'FontSize', 14, 'Box','on');
end

ylabel(ax(1), 'Recursive Polynomial Method', 'Interpreter','latex');

% ---------------- right convergence panel ----------------
ax4 = nexttile(t,4);
semilogy(nList, convMax, 'o-', 'LineWidth', 2, 'MarkerSize', 6);
grid on;
xlabel('$n$', 'Interpreter','latex');
ylabel('$\|err\|_{\infty}$', 'Interpreter','latex');
title('Convergence (Laplace DLP)', 'Interpreter','latex');
set(ax4, 'FontSize', 14, 'Box','on');

if useSquareConv
    axis square;
end

% ---------------- colorbar ----------------
cb = colorbar(ax(end));
cb.Location = 'southoutside';
cb.Label.Interpreter = 'latex';
cb.Label.String = '';
cb.FontSize = 14;

drawnow;
pos = cb.Position;
pos(1) = 0.10;
pos(2) = 0.055;   % 再往下，但仍保留 tick labels 空间
pos(3) = 0.80;
cb.Position = pos;

fprintf('Done.\n');

end

% ======================================================================
% helpers
% ======================================================================

function [logErrAll, nList, nN, nP] = load_zerror_any(matFile, csvFile)
% Prefer MAT: fields nList, logErrAllK, logErrMask, sentinel.
if isfile(matFile)
    S = load(matFile);
    if isfield(S,'nList') && isfield(S,'logErrAllK')
        nList = double(S.nList(:));
        logErrAll = double(S.logErrAllK);

        if isfield(S,'logErrMask')
            mask = double(S.logErrMask);
            if isequal(size(mask), size(logErrAll))
                logErrAll(mask==1) = NaN;
            end
        end
        if isfield(S,'sentinel')
            sentinel = double(S.sentinel);
            logErrAll(logErrAll==sentinel) = NaN;
        end

        [nN,nP] = size(logErrAll);
        return;
    end
end

% CSV fallback (ROBUST): parse cell-by-cell to double
if ~isfile(csvFile)
    error('Cannot find zerror: neither %s nor %s', matFile, csvFile);
end

T = readcell(csvFile);

if isempty(T) || size(T,1) < 2
    error('zerror CSV seems empty.');
end

% parse nList from first column
nList = nan(size(T,1)-1,1);
for i = 2:size(T,1)
    nList(i-1) = cell_to_double(T{i,1});
end

% parse data block
nN = size(T,1)-1;
nP = size(T,2)-1;
logErrAll = nan(nN, nP);

for i = 1:nN
    for j = 1:nP
        logErrAll(i,j) = cell_to_double(T{i+1,j+1});
    end
end

% trim trailing all-NaN columns
while nP > 1 && all(isnan(logErrAll(:,nP)))
    logErrAll(:,nP) = [];
    nP = nP - 1;
end

end

function [Nx, Ny, xRe, yIm] = load_zdataxy_any(matFile, csvFile, nP_expected)
if isfile(matFile)
    D = load(matFile);

    if isfield(D,'Nx') && isfield(D,'Ny') && isfield(D,'wevalRe') && isfield(D,'wevalIm')
        Nx = double(D.Nx);
        Ny = double(D.Ny);
        xRe = double(D.wevalRe(:));
        yIm = double(D.wevalIm(:));
        sanity_nP(Nx,Ny,xRe,yIm,nP_expected);
        return;
    end

    if isfield(D,'Expression1')
        E = D.Expression1;
        if isnumeric(E) && ~isreal(E)
            pts = E(:);
            xRe = double(real(pts));
            yIm = double(imag(pts));
            [Nx,Ny] = infer_dims(numel(pts));
            sanity_nP(Nx,Ny,xRe,yIm,nP_expected);
            return;
        end
        if iscell(E) && numel(E)==1 && isnumeric(E{1}) && ~isreal(E{1})
            pts = E{1}(:);
            xRe = double(real(pts));
            yIm = double(imag(pts));
            [Nx,Ny] = infer_dims(numel(pts));
            sanity_nP(Nx,Ny,xRe,yIm,nP_expected);
            return;
        end
    end
end

% CSV fallback robust
if ~isfile(csvFile)
    error('Cannot find zdataxy: neither %s nor %s', matFile, csvFile);
end

T = readcell(csvFile);
if size(T,2) < 2
    error('zdataxy CSV must have at least 2 columns (x,y).');
end

xRe = nan(size(T,1)-1,1);
yIm = nan(size(T,1)-1,1);
for i = 2:size(T,1)
    xRe(i-1) = cell_to_double(T{i,1});
    yIm(i-1) = cell_to_double(T{i,2});
end

[Nx,Ny] = infer_dims(numel(xRe));
sanity_nP(Nx,Ny,xRe,yIm,nP_expected);
end

function [convMax, ok] = load_conv_any(matFile, csvFile, nList)
ok = false;
convMax = [];

if isfile(matFile)
    C = load(matFile);
    if isfield(C,'convMax')
        convMax = double(C.convMax(:));
        if isfield(C,'convMask')
            m = double(C.convMask(:));
            if numel(m)==numel(convMax)
                convMax(m==1)=NaN;
            end
        end
        if isfield(C,'sentinel')
            s = double(C.sentinel);
            convMax(convMax==s) = NaN;
        end
        ok = true;
        return;
    end
end

if isfile(csvFile)
    T = readcell(csvFile);
    ncol = nan(size(T,1)-1,1);
    vcol = nan(size(T,1)-1,1);
    for i = 2:size(T,1)
        ncol(i-1) = cell_to_double(T{i,1});
        vcol(i-1) = cell_to_double(T{i,2});
    end

    convMax = nan(numel(nList),1);
    for i = 1:numel(nList)
        j = find(abs(ncol - nList(i)) < 1e-12, 1);
        if ~isempty(j)
            convMax(i) = vcol(j);
        end
    end
    ok = true;
end
end

function v = cell_to_double(c)
% robust cell -> double
if isnumeric(c)
    v = double(c);
    return;
end

if isstring(c) || ischar(c)
    s = strtrim(string(c));
    if s == "" || lower(s) == "nan" || lower(s) == "missing"
        v = NaN;
        return;
    end
    v = str2double(s);
    if isnan(v)
        s2 = replace(s, ",", "");
        v = str2double(s2);
    end
    return;
end

v = NaN;
end

function [Nx,Ny] = infer_dims(nP)
if nP == 10100
    Nx = 101; Ny = 100; return;
end
if nP == 300
    Nx = 20; Ny = 15; return;
end

Nx = round(sqrt(nP));
Ny = nP / Nx;
if abs(Ny - round(Ny)) < 1e-12
    Ny = round(Ny);
else
    if mod(nP,100)==0
        Ny = 100;
        Nx = nP / 100;
    else
        Ny = nP;
        Nx = 1;
    end
end
end

function sanity_nP(Nx,Ny,xRe,yIm,nP_expected)
nP = numel(xRe);
if numel(yIm) ~= nP
    error('zdataxy: length mismatch.');
end
if ~isempty(nP_expected) && nP ~= nP_expected
    warning('zdataxy: nP=%d but zerror expects nP=%d. Plot may be wrong.', nP, nP_expected);
end
if Nx*Ny ~= nP
    warning('zdataxy: Nx*Ny=%d but nP=%d. reshape may be off.', Nx*Ny, nP);
end
end