% Risk aware importation rates on superspreading potential
clearvars; clc; close all; tic;

% Assumptions and notes
% - x introduction, y infections, n-x susceptibles
% - no structure (random mixing) but n-x dilutes beta
% - heterogeneity in infectiousness param by R0 and k
% - Dirichlet sample weightings for imports only
% - reproduces Fig 4 and Fig S3 (k = 0.1) and Fig S4 (k = 10)
% - for larger infections tau = d and k = 0.1

% Directory of some main code and plotting options
thisDir = cd; cd('main'); mainDir = cd;
cd(thisDir); addpath(genpath(mainDir));
% Default plotting options
[grey1, grey2, cmap, fnt] = defaultSet(10);

% All event sizes (susceptibles + x = n) and size weight
n = 5:50; m = length(n); nsize = n./sum(n);
% Domains for each that x or y can take for any n
dom = arrayfun(@(x) 0:x, n, 'UniformOutput', false);

% Population R0 and nSamps
R0 = 3; nSamps = 10000;
% Infectious period d, event time tau, R limit
d = 5; tau = d; Rlim = R0*tau/d;

% Dispersion k and prevalence
k = 0.1; rhoM = 0.01;
% Weight the rhoM over n with sorted Dirichlet
whet = [3 1 0.5]; lenw = length(whet);

%% Deterministic statistics of infections and R with k

% Dirichlet weighted import rate statistics
Rmh = cell(1, lenw); ymh = Rmh; Rvh = Rmh; yvh = Rmh;
% Fixed prevalence import rate statistics
Rm = Rmh; Rv = Rm; ym = Rm; yv = Rm; rhoW = Rm;
% Variance to mean ratios
Rvm = Rmh; yvm = Rmh; Rvmh = Rmh; yvmh = Rmh;
% Ratio of mean infections to that from prevalence and quantiles
ymrho = Rm; Rmrho = Rm; yvrho = Rm; Rvrho = Rm; quant = [0.025, 0.5, 0.975];

% For largest and smallest k and every Dirichlet weight compute R and E[y]
for j = 1:lenw
    % Specific k and weight set to be evaluated and get statistics
    [ym{j}, yv{j}, Rm{j}, Rv{j}, ymh{j}, yvh{j}, Rmh{j}, Rvh{j}, rhoW{j}] =...
        weightDirichlet(n, nsize, m, k, whet(j), nSamps, rhoM, dom, Rlim);

    % Variance to mean ratios
    Rvm{j} = Rv{j}./Rm{j}; Rvmh{j} = Rvh{j}./Rmh{j};
    yvm{j} = yv{j}./ym{j}; yvmh{j} = yvh{j}./ymh{j};
    % Ratios to fixed prevalence case
    Rmrho{j} = Rmh{j}./Rm{j}; ymrho{j} = ymh{j}./ym{j};
    Rvrho{j} = Rvh{j}./Rv{j}; yvrho{j} = yvh{j}./yv{j};
end


% Get quantiles of reproduction numbers and infection statistics
Rmhq = cellfun(@(x) quantile(x, quant), Rmh, 'UniformOutput', false);
ymhq = cellfun(@(x) quantile(x, quant), ymh, 'UniformOutput', false);
Rvhq = cellfun(@(x) quantile(x, quant), Rvh, 'UniformOutput', false);
yvhq = cellfun(@(x) quantile(x, quant), yvh, 'UniformOutput', false);

% Quantiles of sampled weights and variance to mean ratios
rhoWq = cellfun(@(x) quantile(x, quant), rhoW, 'UniformOutput', false);
Rvmhq = cellfun(@(x) quantile(x, quant), Rvmh, 'UniformOutput', false);
yvmhq = cellfun(@(x) quantile(x, quant), yvmh, 'UniformOutput', false);

% Sum total infections across the event sizes
ysum = cellfun(@(x) sum(x,2), ymh, 'UniformOutput', false);
% Ratios over fixed prevalence quantiles
Rmrhoq = cellfun(@(x) quantile(x, quant), Rmrho, 'UniformOutput', false);
ymrhoq = cellfun(@(x) quantile(x, quant), ymrho, 'UniformOutput', false);
Rvrhoq = cellfun(@(x) quantile(x, quant), Rvrho, 'UniformOutput', false);
yvrhoq = cellfun(@(x) quantile(x, quant), yvrho, 'UniformOutput', false);

% Mean statistics
ymavg = cellfun(@mean, ymh, 'UniformOutput', false);
yvavg = cellfun(@mean, yvh, 'UniformOutput', false);
Rmavg = cellfun(@mean, Rmh, 'UniformOutput', false);
Rvavg = cellfun(@mean, Rvh, 'UniformOutput', false);
rhoavg = cellfun(@mean, rhoW, 'UniformOutput', false);

% Find optimal n where mean crossover occurs
idcross = zeros(1, lenw); ncross = idcross;
for j = 1:lenw
    idcross(j) = find(abs(rhoavg{j} - rhoM) == min(abs(rhoavg{j} - rhoM)));
    ncross(j) = n(idcross(j));
end

%% Figure 4: risk awareness and ratios to null model

figure('Position', [10 10 800 800]);
cols = {'b', 'g', 'r'};

% Size biased Dirichlet weights
subplot(2, 2, 1); hold on;
for i = 1:3
    % Risk aware import rates
    plotCIRaw(n', rhoWq{i}(2, :)', rhoWq{i}(1, :)', rhoWq{i}(3, :)', cols{i});
end
% Prevalence
plot(n, rhoM*ones(1, m), 'k--', 'LineWidth', 2); h = gca;
% Crossover points
for i = 1:3
    % Crossover points
    plot([ncross(i) ncross(i)], h.YLim, '--', 'Color', cols{i}, 'LineWidth', 2);
end
ylabel('risk aware import rate $\epsilon(n)$', 'FontSize', fnt);
xlim([min(n) max(n)]); xlabel('event sizes $n$', 'FontSize', fnt);
hold off;  grid off; box off;
legend('', ['$r = $ ' num2str(whet(1))], '', ['$r = $ ' num2str(whet(2))],...
    '', ['$r = $ ' num2str(whet(3))], 'FontSize', fnt);

% Mean y infections due to risk awareness
subplot(2, 2, 2); hold on;
for i = 1:3
    plotCIRaw(n', ymhq{i}(2, :)', ymhq{i}(1, :)', ymhq{i}(3, :)', cols{i});
end
plot(n, ym{1}, 'k--', 'LineWidth', 2);
hold off;  grid off; box off;
xlim([min(n) max(n)]); xlabel('event sizes $n$', 'FontSize', fnt);
ylabel('mean infections, E$[y | \epsilon]$', 'FontSize', fnt);

% Mean ratios
axes('Position', [0.65, 0.73, 0.17, 0.17]);
hold on;
for i = 1:3
    plot(n, ymavg{i}./ym{1}, 'Color', cols{i}, 'LineWidth', 2);
end
plot(n, ones(1, m), 'k--', 'LineWidth', 2);
hold off; grid off; box off;
xlim([min(n) max(n)]); xlabel('$n$', 'FontSize', fnt);
ylabel('E$[y | \epsilon]$ ratio', 'FontSize', fnt);

% Variance of y infections due to risk awareness
subplot(2, 2, 3); hold on;
for i = 1:3
    plotCIRaw(n', yvhq{i}(2, :)', yvhq{i}(1, :)', yvhq{i}(3, :)', cols{i});
end
plot(n, yv{1}, 'k--', 'LineWidth', 2);
hold off;  grid off; box off;
ylabel('var infections, V$[y | \epsilon]$', 'FontSize', fnt);
xlim([min(n) max(n)]); xlabel('event sizes $n$', 'FontSize', fnt);

% Var ratios
axes('Position', [0.21, 0.26, 0.17, 0.17]);
hold on;
for i = 1:3
    plot(n, yvavg{i}./yv{1}, 'Color', cols{i}, 'LineWidth', 2);
end
plot(n, ones(1, m), 'k--', 'LineWidth', 2);
hold off; grid off; box off;
xlim([min(n) max(n)]); xlabel('$n$', 'FontSize', fnt);
ylabel('V$[y | \epsilon]$ ratio', 'FontSize', fnt);


% % For reviewer 2 response
% % Variance of y infections due to risk awareness
% subplot(2, 2, 3); hold on;
% for i = 1:3
%     plotCIRaw(n', Rmhq{i}(2, :)', Rmhq{i}(1, :)', Rmhq{i}(3, :)', cols{i});
% end
% plot(n, Rm{1}, 'k--', 'LineWidth', 2);
% plot(n, Rlim*ones(size(n)), 'k', 'LineWidth', 2);
% hold off;  grid off; box off;
% ylabel('mean event R, E$[R(x) | \epsilon]$', 'FontSize', fnt);
% xlim([min(n) max(n)]); xlabel('event sizes $n$', 'FontSize', fnt);
% 
% % Var ratios
% axes('Position', [0.21, 0.26, 0.17, 0.17]);
% hold on;
% for i = 1:3
%     plot(n, Rmavg{i}./Rm{1}, 'Color', cols{i}, 'LineWidth', 2);
% end
% plot(n, ones(1, m), 'k--', 'LineWidth', 2);
% hold off; grid off; box off;
% xlim([min(n) max(n)]); xlabel('$n$', 'FontSize', fnt);
% ylabel('E$[R(x) | \epsilon]$ ratio', 'FontSize', fnt);


% Histogram of infections in total
ysumH = cell2mat(ysum);
ymin = min(min(ysumH)); ymax = max(max(ysumH));

subplot(2, 2, 4); hold on;
for i = 1:3
    h = histogram(ysumH(:, i), 'Normalization', 'probability',...
        'BinEdges', linspace(ymin, ymax, 100));
    h.EdgeAlpha = 0; h.FaceAlpha = 0.4; h.FaceColor = cols{i};
end
hold off;  grid off; box off;
xlabel('total mean infections', 'FontSize', fnt);

%% Fig S3: statistics of event R and infections

figure('Position', [10 10 800 800]);

% Variance to mean of R 
subplot(3, 2, 1); hold on;
for i = 1:3
    plotCIRaw(n', Rvmhq{i}(2, :)', Rvmhq{i}(1, :)', Rvmhq{i}(3, :)', cols{i});
    % All fixed prevalence statistics will be same
    if i == 3
        plot(n, Rvm{i}, 'k--', 'LineWidth', 2);
    end
end
ylabel('VM of event $R$', 'FontSize', fnt);
hold off;  grid off; box off; xlim([min(n) max(n)]);

% Variance to mean of y 
subplot(3, 2, 2); hold on;
for i = 1:3
    plotCIRaw(n', yvmhq{i}(2, :)', yvmhq{i}(1, :)', yvmhq{i}(3, :)', cols{i});
    % All fixed prevalence statistics will be same
    if i == 3
        plot(n, yvm{i}, 'k--', 'LineWidth', 2);
    end
end
ylabel('VM of infections $y$', 'FontSize', fnt);
hold off;  grid off; box off; xlim([min(n) max(n)]);

% Variances of R 
subplot(3, 2, 3); hold on;
for i = 1:3
    plotCIRaw(n', Rvhq{i}(2, :)', Rvhq{i}(1, :)', Rvhq{i}(3, :)', cols{i});
    % All fixed prevalence statistics will be same
    if i == 3
        plot(n, Rv{i}, 'k--', 'LineWidth', 2);
    end
end
ylabel('var of event $R$', 'FontSize', fnt);
hold off;  grid off; box off; xlim([min(n) max(n)]);

% Variances of y 
subplot(3, 2, 4); hold on;
for i = 1:3
    plotCIRaw(n', yvhq{i}(2, :)', yvhq{i}(1, :)', yvhq{i}(3, :)', cols{i});
    % All fixed prevalence statistics will be same
    if i == 3
        plot(n, yv{i}, 'k--', 'LineWidth', 2);
    end
end
ylabel('var of infections $y$', 'FontSize', fnt);
hold off;  grid off; box off; xlim([min(n) max(n)]);

% Means of R 
subplot(3, 2, 5); hold on;
for i = 1:3
    plotCIRaw(n', Rmhq{i}(2, :)', Rmhq{i}(1, :)', Rmhq{i}(3, :)', cols{i});
    % All fixed prevalence statistics will be same
    if i == 3
        plot(n, Rm{i}, 'k--', 'LineWidth', 2);
    end
end
ylabel('mean event $R$', 'FontSize', fnt);
hold off;  grid off; box off; xlim([min(n) max(n)]);

% Variances of y 
subplot(3, 2, 6); hold on;
for i = 1:3
    plotCIRaw(n', ymhq{i}(2, :)', ymhq{i}(1, :)', ymhq{i}(3, :)', cols{i});
    % All fixed prevalence statistics will be same
    if i == 3
        plot(n, ym{i}, 'k--', 'LineWidth', 2);
    end
end
ylabel('mean infections $y$', 'FontSize', fnt);
xlabel(['$n | \rho = $ ' num2str(rhoM)], 'FontSize', fnt);
hold off;  grid off; box off; xlim([min(n) max(n)]);

legend('', ['$r = $ ' num2str(whet(1))], '', ['$r = $ ' num2str(whet(2))],...
    '', ['$r = $ ' num2str(whet(3))], 'FontSize', fnt);
xlabel(['event sizes $n \, | \, $prevalence $\rho = $ ' num2str(rhoM)...
    ' at dispersion $k = $ ' num2str(k)], 'FontSize', fnt);

% Running time of script
tsim = toc/60; disp(['Run time = ' num2str(tsim)]);