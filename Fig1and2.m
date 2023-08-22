% Event reproduction numbers and clusters for a single event
clearvars; clc; close all; tic;

% Assumptions and notes
% - x introduction, y infections, n-x susceptibles
% - no structure (random mixing) but n-x dilutes beta 
% - heterogeneity in infectiousness param by R0 and k
% - reproduces Fig 1 and Fig 2 of paper (Fig S1 at n = 100)

% Directory of some main code and plotting options
thisDir = cd; cd('main'); mainDir = cd;
cd(thisDir); addpath(genpath(mainDir));
% Default plotting options
[grey1, grey2, cmap, fnt] = defaultSet(10);

% Event size, R0 and domain x or y takes
n = 30; dom = 0:n; ldom = length(dom); R0 = 3;
% Infectious period d, event time tau, R limit
d = 5; tau = d/10; Rlim = R0*tau/d;
% Dispersion k controlling heterogeneity
k = [0.1 0.25 0.5 1 5 10]; lenk = length(k);

% Mean prevalence and samples
rhoM = 0.05; nSamps = 100000;
% Domain for sum of R0 values 
R0dom = 0.01:0.01:200; lRdom = length(R0dom);

%% Transmission with heterogeneity as a function of imports

% Prob of y given x and of transmission for heterogeneous R0
Pyx = cell(1, lenk); ptransRx = Pyx; pinfx = Pyx;
% Event R based on import number and x that maximises y 
Rx = Pyx; Vx = Pyx; R1 = zeros(1, lenk); V1 = R1; 
% Samples of R and y, avg y for each x and max x
Rsamp = Rx; ysamp = Rx; yavg = Pyx; xmaxAvg = R1; 

% Generate P(y|x,n) for all possible x with heterogeneity
for j = 1:lenk
    % Event reproduction numbers R(x) as function of x
    [Rx{j}, pinfx{j}, Vx{j}, ysamp{j}, Rsamp{j}] = getRxSSEsamp(n, R0, tau, d, k(j), nSamps);

    % Standard event R (1 import) and its variance
    R1(j) = Rx{j}(2); V1(j) = Vx{j}(2);

    % Probability of y new infections for all possible x
    Pyx{j} = getygivenxSSE(dom, n, pinfx{j});
    % Check that average num infections matches formulations
    yavg{j} = Pyx{j}*dom'; yavgR = Rx{j}.*dom;

    % Import count with maximum y on average
    [~, yi] = max(yavg{j}); xmaxAvg(j) = dom(yi);
end
% Find prevalence with max average x (same as xmaxAvg/n)
rhoMax = binofit(xmaxAvg, n);

% Statistics from infection clusters
ysampM = cellfun(@mean, ysamp, 'UniformOutput', false);
ysampV = cellfun(@var, ysamp, 'UniformOutput', false);
ysampM = cell2mat(ysampM')'; ysampV = cell2mat(ysampV')';

% Convert to matrices
yavg = cell2mat(yavg); Rx = cell2mat(Rx')'; Vx = cell2mat(Vx')'; 

% Statistics from event R
RsampM = cellfun(@mean, Rsamp, 'UniformOutput', false);
RsampV = cellfun(@var, Rsamp, 'UniformOutput', false);
RsampM = cell2mat(RsampM')'; RsampV = cell2mat(RsampV')';
% Variance to mean ratios
RsampVM = RsampV./RsampM; ysampVM = ysampV./ysampM; 

%% Figure 1: R(x) and E[y|x,n] and VM ratios 

% Basic statistics of small R and y
figure('Renderer', 'painters', 'Position', [10 10 800 600]);
% Plot how Rev(x) and avg y|x depends on every x at different k
subplot(3, 2, [1 3]); hold on;
stairs(dom, RsampM(:, 2:end-1), 'LineWidth', 2, 'Color', grey1);
stairs(dom, RsampM(:, 1), 'LineWidth', 2, 'Color', 'b');
stairs(dom, RsampM(:, end), 'LineWidth', 2, 'Color', 'r');
stairs(dom, Rlim*ones(size(dom)), 'k--', 'LineWidth', 1);
hold off; box off; grid off; ylim([0 0.01+Rlim]);
ylabel('E$_{het}[R(x)]$', 'FontSize', fnt);

subplot(3, 2, [2 4]); hold on;
stairs(dom, ysampM(:, 2:end-1), 'LineWidth', 2, 'Color', grey1);
stairs(dom, ysampM(:, 1), 'LineWidth', 2, 'Color', 'b');
stairs(dom, ysampM(:, end), 'LineWidth', 2, 'Color', 'r');
stairs(dom, n-dom, 'k--', 'LineWidth', 1);
hold off; box off; grid off; ylim([0 6]);
ylabel('E$_{het}$[E$[y|x, n]$]', 'FontSize', fnt);

% Plot VMs of Rev(x) and avg y|x depends on every x at different k
subplot(3, 2, 5); hold on;
stairs(dom, RsampVM(:, 2:end-1), 'LineWidth', 2, 'Color', grey1);
stairs(dom, RsampVM(:, 1), 'LineWidth', 2, 'Color', 'b');
stairs(dom, RsampVM(:, end), 'LineWidth', 2, 'Color', 'r');
hold off; box off; grid off;
ylabel('VM$_{het}[R(x)]$', 'FontSize', fnt);
xlabel(['$x | n = $ ' num2str(n)], 'FontSize', fnt); 
legend('', '', '', '', ['$k = $ ' num2str(k(1))], ['$k = $ ' num2str(k(end))])

subplot(3, 2, 6); hold on;
stairs(dom, ysampVM(:, 2:end-1), 'LineWidth', 2, 'Color', grey1);
stairs(dom, ysampVM(:, 1), 'LineWidth', 2, 'Color', 'b');
stairs(dom, ysampVM(:, end), 'LineWidth', 2, 'Color', 'r');
hold off; box off; grid off;
ylabel('VM$_{het}$[E$[y|x]$]', 'FontSize', fnt);
xlabel(['$x | n = $ ' num2str(n)], 'FontSize', fnt);


%% Influence of import models on clustering statistics

% Import prob function for a given mean
Pxrhofn = @(alpha) binopdf(dom, n, alpha);
% Null scenario of sampling from prevalence or uniform
Pxfix = Pxrhofn(rhoM); Pxunif = (1/(n+1))*ones(1, n+1);

% Mean infections for this event size 
yfix = Pxfix*yavg; Rfix = Pxfix*Rx;
% Only focus on most heterogeneous and most deterministic case
yHet = ysamp{1}; yDet = ysamp{end}; RHet = Rsamp{1}; RDet = Rsamp{end};

% Curve if increase prevalence of population
rho = 0.01:0.01:0.1; lenrho = length(rho);
Pxrho = zeros(lenrho, n+1);
for i = 1:lenrho
    Pxrho(i, :) = Pxrhofn(rho(i));
end

% Means of y for these alpha from samples
yHmeans = zeros(1, lenrho); yDmeans = yHmeans; yHvars = yHmeans; yDvars = yHmeans;
RHmeans = yHmeans; RDmeans = yHmeans; RHvars = yHmeans; RDvars = yHmeans;
% Survival probabilities indicating superspreading
syD = cell(1, lenrho); syH = syD; syDx = syD; syHx = syD; yHi = syD;
sRD = cell(1, lenrho); sRH = sRD; sRDx = sRD; sRHx = sRD; yDi = syD;

for i = 1:lenrho
    % Thinned by probability of imports
    [yHi{i}, RHi] = impWeight(n, yHet, RHet, Pxrho(i, :), nSamps);
    [yDi{i}, RDi] = impWeight(n, yDet, RDet, Pxrho(i, :), nSamps);
    % Statistics and survival probabilities
    yHmeans(i) = mean(yHi{i}); yDmeans(i) = mean(yDi{i});
    yHvars(i) = var(yHi{i}); yDvars(i) = var(yDi{i});
    [syD{i}, syDx{i}] = ecdf(yDi{i}); [syH{i}, syHx{i}] = ecdf(yHi{i});
    RHmeans(i) = mean(RHi); RDmeans(i) = mean(RDi);
    RHvars(i) = var(RHi); RDvars(i) = var(RDi);
    [sRD{i}, sRDx{i}] = ecdf(RDi); [sRH{i}, sRHx{i}] = ecdf(RHi);
end

% Define bin edges based on domain
domEdg = sort(unique([dom-0.5, dom+0.5]));
% Shorter domain for testing
domSh = 0:round(0.66*n); lSh = length(domSh);
id = [1 2 round(lenrho/2) lenrho];

% Integer probabilities P(y = c) and cumulative P(y <= c) for clusters
pH = zeros(lenrho, ldom); pD = pH; cpH = pH; cpD = pH;
for i = 1:lenrho
    % Heterogeneous case
    pH(i, :) = histcounts(yHi{i}, 'Normalization', 'probability', 'BinEdges', domEdg);
    cpH(i, :) = cumsum(pH(i, :));
    % Homoegeneous case
    pD(i, :) = histcounts(yDi{i}, 'Normalization', 'probability', 'BinEdges', domEdg);
    cpD(i, :) = cumsum(pD(i, :));
end

%% Figure 2: cluster size histograms and survival curves

figure('Renderer', 'painters', 'Position', [10 10 800 800]);
id1 = [1 round(lenrho/2) lenrho];
% Survival tail probabilities with k and alpha
cols = {'b', grey1, grey1, 'r'};
subplot(3, 6, [1:3 7:9]);  hold on;
for i = 1:4
    stairs(syHx{id(i)}, log(1-syH{id(i)}), '-', 'LineWidth', 2, 'Color', cols{i});
    stairs(syDx{id(i)}, log(1-syD{id(i)}), '--', 'LineWidth', 2, 'Color', cols{i});
end
hold off; box off; grid off;
ylabel('$\log$ P$(y \geq c)$', 'FontSize', fnt);
xlabel(['$c | n = $ ' num2str(n)], 'FontSize', fnt);
% Mean and variances as functions of alpha
subplot(3, 6, [4:6 10:12]);  hold on;
for i = 1:4
    stairs(sRHx{id(i)}, log(1-sRH{id(i)}), '-', 'LineWidth', 2, 'Color', cols{i});
    stairs(sRDx{id(i)}, log(1-sRD{id(i)}), '--', 'LineWidth', 2, 'Color', cols{i});
end
hold off; box off; grid off;
ylabel('$\log$ P$(R \geq c)$', 'FontSize', fnt);
xlabel(['$c | n = $ ' num2str(n)], 'FontSize', fnt);
legend([num2str(rho(id(1))) ',' num2str(k(1))], [num2str(rho(id(1))) ',' ...
    num2str(k(end))], '', '', '', '', [num2str(rho(id(end))) ',' ... 
    num2str(k(1))], [num2str(rho(id(end))) ',' num2str(k(end))], FontSize = fnt);

% Histogram of infections at lowest, middle and peak prevalence
for i = 1:3
    subplot(3, 6, [12+(2*i-1) 12+2*i]); hold on;
    h = histogram(yHi{id1(i)}, 'Normalization', 'probability', 'BinMethod', 'integers');
    h.EdgeAlpha = 0; h.FaceAlpha = 0.4;
    h = histogram(yDi{id1(i)}, 'Normalization', 'probability', 'BinMethod', 'integers');
    h.EdgeAlpha = 0; h.FaceAlpha = 0.4; h = gca;
    plot(yHvars(id1(i))*ones(1, 2), h.YLim, 'b--', 'LineWidth', 2);
    plot(yDvars(id1(i))*ones(1, 2), h.YLim, 'r--', 'LineWidth', 2);
    hold off; box off; grid off;
    xlabel(['$y|\rho = $ ' num2str(rho(id1(i)))], 'FontSize', fnt); 
    ylabel('P$(y)$', 'FontSize', fnt); 
    xlim([-1 9]); ylim([0 1.05]);
    if i == 3
        legend(['$k = $ ' num2str(k(1))], ['$k = $ ' num2str(k(end))],...
            FontSize = fnt-1, Location="best");
    elseif i == 1
        legend('', '', ['V$[y] | k = $ ' num2str(k(1))], ['V$[y] | k = $ ' ...
            num2str(k(end))], FontSize = fnt-1, Location = "best");
    end
end

% Running time of script
tsim = toc/60; disp(['Run time = ' num2str(tsim)]);
