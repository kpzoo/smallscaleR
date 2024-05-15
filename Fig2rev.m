% Event reproduction numbers and clusters for a single event
clearvars; clc; close all; tic;

% Assumptions and notes
% - x introduction, y infections, n-x susceptibles
% - no structure (random mixing) but n-x dilutes beta 
% - heterogeneity in infectiousness param by R0 and k
% - reproduces Fig 2 at R0 = 3, tau = d/10, rho = 0.05 
% - reproduces Fig S1 at R0 = 6, tau = d/2, rho = 0.1

% Directory of some main code and plotting options
thisDir = cd; cd('main'); mainDir = cd;
cd(thisDir); addpath(genpath(mainDir));
% Default plotting options
[grey1, grey2, cmap, fnt] = defaultSet(10);

% Event size, R0 and domain x or y takes
n = 30; dom = 0:n; ldom = length(dom); R0 = 6;
% Infectious period d, event time tau, R limit
d = 5; tau = d/2; Rlim = R0*tau/d;
% Dispersion k controlling heterogeneity
k = [0.1 0.25 0.5 1 5 10]; lenk = length(k);

% Mean prevalence and samples
rhoM = 0.1; nSamps = 100000;
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

%% Figure 2: R(x) and E[y|x,n] and VM ratios 

% Basic statistics of small R and y
figure('Position', [10 10 800 600]);
% Plot how Rev(x) and avg y|x depends on every x at different k
subplot(3, 2, [1 3]); hold on;
stairs(dom, RsampM(:, 2:end-1), 'LineWidth', 2, 'Color', grey1);
stairs(dom, RsampM(:, 1), 'LineWidth', 2, 'Color', 'b');
stairs(dom, RsampM(:, end), 'LineWidth', 2, 'Color', 'r');
stairs(dom, Rlim*ones(size(dom)), 'k-', 'LineWidth', 2);
hold off; box off; grid off; ylim([0 0.01+Rlim]);
ylabel('mean event R, $R(x)$', 'FontSize', fnt);

subplot(3, 2, [2 4]); hold on;
stairs(dom, ysampM(:, 2:end-1), 'LineWidth', 2, 'Color', grey1);
stairs(dom, ysampM(:, 1), 'LineWidth', 2, 'Color', 'b');
stairs(dom, ysampM(:, end), 'LineWidth', 2, 'Color', 'r');
stairs(dom, n-dom, 'k-', 'LineWidth', 2);
hold off; box off; grid off; %ylim([0 6]);
ylabel('mean exp infections, E$[y | x]$', 'FontSize', fnt);

% Plot VMs of Rev(x) and avg y|x depends on every x at different k
subplot(3, 2, 5); hold on;
stairs(dom, RsampVM(:, 2:end-1), 'LineWidth', 2, 'Color', grey1);
stairs(dom, RsampVM(:, 1), 'LineWidth', 2, 'Color', 'b');
stairs(dom, RsampVM(:, end), 'LineWidth', 2, 'Color', 'r');
hold off; box off; grid off;
ylabel('VM of $R(x)$', 'FontSize', fnt);
xlabel(['imports $x |$event size $n = $ ' num2str(n)], 'FontSize', fnt); 
legend('', '', '', '', ['$k = $ ' num2str(k(1))], ['$k = $ ' num2str(k(end))]);

subplot(3, 2, 6); hold on;
stairs(dom, ysampVM(:, 2:end-1), 'LineWidth', 2, 'Color', grey1);
stairs(dom, ysampVM(:, 1), 'LineWidth', 2, 'Color', 'b');
stairs(dom, ysampVM(:, end), 'LineWidth', 2, 'Color', 'r');
hold off; box off; grid off;
ylabel('VM of E$[y|x]$', 'FontSize', fnt);
xlabel('imports $x$', 'FontSize', fnt);
%xlabel(['imports $x | n = $ ' num2str(n)], 'FontSize', fnt);

% Running time of script
tsim = toc/60; disp(['Run time = ' num2str(tsim)]);
