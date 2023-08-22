% Risk aware importation rates on superspreading potential
clearvars; clc; close all; tic;

% Assumptions and notes
% - x introduction, y infections, n-x susceptibles
% - no structure (random mixing) but n-x dilutes beta 
% - heterogeneity in infectiousness param by R0 and k
% - Dirichlet sample weightings for imports only
% - reproduces Fig 3 of paper and supplement Fig S2

% Directory of some main code and plotting options
thisDir = cd; cd('main'); mainDir = cd;
cd(thisDir); addpath(genpath(mainDir));
% Default plotting options
[grey1, grey2, cmap, fnt] = defaultSet(10);

% All event sizes (susceptibles + x = n) and size weight
n = 5:50; m = length(n); nsize = n./sum(n);
%n = [10*ones(1, 25) 25*ones(1, 10) 50]; m = length(n); nsize = n./sum(n);
% Domains for each that x or y can take for any n
dom = arrayfun(@(x) 0:x, n, 'UniformOutput', false);

% Mean prevalence and population R0 and nSamps
rhoM = 0.01; R0 = 3; nSamps = 10000;
% Infectious period d, event time tau, R limit
d = 5; tau = d/10; Rlim = R0*tau/d;

% Dispersion k controlling heterogeneity
k = [10 0.1]; lenk = length(k);
% Weight the rhoM over n with sorted Dirichlet 
whet = [5 1 0.5]; lenw = length(whet); 

%% Deterministic statistics of infections and R with k

% Dirichlet weighted import rate statistics
Rmh = cell(lenk, lenw); ymh = Rmh; Rvh = Rmh; yvh = Rmh;
% Fixed prevalence import rate statistics 
Rm = Rmh; Rv = Rm; ym = Rm; yv = Rm; rhoW = Rm;
% Variance to mean ratios
Rvm = Rmh; yvm = Rmh; Rvmh = Rmh; yvmh = Rmh;
% Ratio of mean infections to that from prevalence and quantiles
ymrho = Rm; Rmrho = Rm; quant = [0.025, 0.5, 0.975];

% For largest and smallest k and every Dirichlet weight compute R and E[y]
for i = 1:lenk
    for j = 1:lenw
        % Specific k and weight set to be evaluated and get statistics
        [ym{i,j}, yv{i,j}, Rm{i,j}, Rv{i,j}, ymh{i,j}, yvh{i,j}, Rmh{i,j}, Rvh{i,j},...
            rhoW{i,j}] = weightDirichlet(n, nsize, m, k(i), whet(j), nSamps, rhoM, dom, Rlim);

        % Variance to mean ratios
        Rvm{i,j} = Rv{i,j}./Rm{i,j}; Rvmh{i,j} = Rvh{i,j}./Rmh{i,j};
        yvm{i,j} = yv{i,j}./ym{i,j}; yvmh{i,j} = yvh{i,j}./ymh{i,j};
        % Ratios to fixed prevalence case
        Rmrho{i,j} = Rmh{i,j}./Rm{i,j}; ymrho{i,j} = ymh{i,j}./ym{i,j};
    end
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

%% Figure 3: E_het[E_imp[R]] and E[y|x,n] for risk awareness 

% Importation rate and mean risk statistics and ratios
figure('Renderer', 'painters', 'Position', [10 10 800 800]);
cols = {'b', 'g', 'r'}; 
% Size biased Dirichlet weights (used second k but results same)
subplot(4, 2, [1 3]);
for i = 1:3
    hold on;
    plotCIRaw(n', rhoWq{2,i}(2, :)', rhoWq{2,i}(1, :)', rhoWq{2,i}(3, :)', cols{i});
    % All fixed prevalence statistics will be same
    if i == 3
        plot(n, rhoM*ones(1, m), 'k--', 'LineWidth', 2);
    end
    hold off;  grid off; box off;
    ylabel('$\epsilon(n) | r$', 'FontSize', fnt);
    xlim([min(n) max(n)]);
end
legend('', ['$r = $ ' num2str(whet(1))], '', ['$r = $ ' num2str(whet(2))],...
            '', ['$r = $ ' num2str(whet(3))], ['$\rho = $ ' num2str(rhoM)], 'FontSize', fnt-1);

% Ratio of y infections on average due to alpha
subplot(4, 2, [2 4]);
for i = 1:3
    hold on;
    plotCIRaw(n', ymrhoq{2,i}(2, :)', ymrhoq{2,i}(1, :)', ymrhoq{2,i}(3, :)', cols{i});
    hold off;  grid off; box off;
    ylabel('E$[y | \epsilon]$ E$[y | \rho]^{-1}$', 'FontSize', fnt);
    xlim([min(n) max(n)]);
end

% Means of R under most and least heterogenous cases
for j = 1:2
    subplot(4, 2, 4+2*j-1); hold on;
    for i = 1:3
        plotCIRaw(n', Rmhq{j,i}(2, :)', Rmhq{j,i}(1, :)', Rmhq{j,i}(3, :)', cols{i});
        % All fixed prevalence statistics will be same 
        if i == 3
            plot(n, Rm{j,i}, 'k--', 'LineWidth', 2);
        end
    end
    ylabel(['E$[R] | k = $ ' num2str(k(j))], 'FontSize', fnt);
    if j == 2
        xlabel(['$n | \rho = $ ' num2str(rhoM)], 'FontSize', fnt);
    end
     hold off; grid off; box off; xlim([min(n) max(n)]);
end
% Means of y under most and least heterogenous cases
for j = 1:2
    subplot(4, 2, 4+2*j); hold on;
    for i = 1:3
        plotCIRaw(n', ymhq{j,i}(2, :)', ymhq{j,i}(1, :)', ymhq{j,i}(3, :)', cols{i});
        % All fixed prevalence statistics will be same 
        if i == 3
            plot(n, ym{j,i}, 'k--', 'LineWidth', 2);
        end
    end
    ylabel(['E$[y|n] | k = $ ' num2str(k(j))], 'FontSize', fnt);
    if j == 2
        xlabel(['$n | \rho = $ ' num2str(rhoM)], 'FontSize', fnt);
    end
    hold off;  grid off; box off; xlim([min(n) max(n)]);
end

%% Fig S2: VM ratios and variances with Dirichlet

figure('Renderer', 'painters', 'Position', [10 10 800 800]);

% Variance to mean of R under most and least heterogenous cases
for j = 1:2
    subplot(4, 2, j); hold on;
    for i = 1:3
        plotCIRaw(n', Rvmhq{j,i}(2, :)', Rvmhq{j,i}(1, :)', Rvmhq{j,i}(3, :)', cols{i});
        % All fixed prevalence statistics will be same 
        if i == 3
            plot(n, Rvm{j,i}, 'k--', 'LineWidth', 2);
        end
    end
    ylabel('VM$[R]$', 'FontSize', fnt);
    hold off;  grid off; box off; xlim([min(n) max(n)]);
    title(['$k = $ ' num2str(k(j))])
end
% Means of y under most and least heterogenous cases
for j = 1:2
    subplot(4, 2, 2+j); hold on;
    for i = 1:3
        plotCIRaw(n', yvmhq{j,i}(2, :)', yvmhq{j,i}(1, :)', yvmhq{j,i}(3, :)', cols{i});
        % All fixed prevalence statistics will be same 
        if i == 3
            plot(n, yvm{j,i}, 'k--', 'LineWidth', 2);
        end
    end
    ylabel('VM$[y]$', 'FontSize', fnt);
    hold off;  grid off; box off; xlim([min(n) max(n)]);
end

% Variance to mean of R under most and least heterogenous cases
for j = 1:2
    subplot(4, 2, 4+j); hold on;
    for i = 1:3
        plotCIRaw(n', Rvhq{j,i}(2, :)', Rvhq{j,i}(1, :)', Rvhq{j,i}(3, :)', cols{i});
        % All fixed prevalence statistics will be same 
        if i == 3
            plot(n, Rv{j,i}, 'k--', 'LineWidth', 2);
        end
    end
    ylabel('V$[R]$', 'FontSize', fnt);
    hold off;  grid off; box off; xlim([min(n) max(n)]);
end
% Variances of y under most and least heterogenous cases
for j = 1:2
    subplot(4, 2, 6+j); hold on;
    for i = 1:3
        plotCIRaw(n', yvhq{j,i}(2, :)', yvhq{j,i}(1, :)', yvhq{j,i}(3, :)', cols{i});
        % All fixed prevalence statistics will be same 
        if i == 3
            plot(n, yv{j,i}, 'k--', 'LineWidth', 2);
        end
    end
    ylabel('V$[y]$', 'FontSize', fnt);
    xlabel(['$n | \rho = $ ' num2str(rhoM)], 'FontSize', fnt);
    hold off;  grid off; box off; xlim([min(n) max(n)]);
end
legend('', ['$r = $ ' num2str(whet(1))], '', ['$r = $ ' num2str(whet(2))],...
            '', ['$r = $ ' num2str(whet(3))], ['$\rho = $ ' num2str(rhoM)], 'FontSize', fnt-1);


% Running time of script
tsim = toc/60; disp(['Run time = ' num2str(tsim)]);