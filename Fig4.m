% Superspreading potential from risk awareness strength and prevalence
clearvars; clc; close all; tic;

% Assumptions and notes
% - x introduction, y infections, n-x susceptibles
% - no structure (random mixing) but n-x dilutes beta 
% - heterogeneity in infectiousness param by R0 and k
% - Dirichlet sample weightings for imports only
% - variations of mean prevalence and rik awareness
% - reproduces Fig 4 of paper and Fig S3 of supplement

% Directory of some main code and plotting options
thisDir = cd; cd('main'); mainDir = cd;
cd(thisDir); addpath(genpath(mainDir));
% Default plotting options
[grey1, grey2, cmap, fnt] = defaultSet(10);

% All event sizes (susceptibles + x = n) and size weight
n = 5:50; m = length(n); nsize = n./sum(n);
% Domains for each that x or y can take for any n
dom = arrayfun(@(x) 0:x, n, 'UniformOutput', false);

% Fixed R0, k and num of samples
R0 = 3; nSamps = 10000; k = 0.5;
% Infectious period d, event time tau, R limit
d = 5; tau = d/10; Rlim = R0*tau/d;

% Mean prevalence range
rho = 0.01:0.01:0.1; lenrho = length(rho);
% Risk awareness strengths (r in paper)
lenw = 10; whet = sort(linspace(0.1, 1, lenw), 'descend');
% Specific whet and rho to take slices and median
whet0 = 0.5; rho0 = 0.05; quant = 0.5;

%% Infections and R at various prevalence for given Dirichlet shape

% Dirichlet weighted import rate statistics
Rmh1 = cell(1, lenrho); ymh1 = Rmh1; Rvh1 = Rmh1; yvh1 = Rmh1;
% Fixed prevalence import rate statistics and weights
Rm1 = Rmh1; Rv1 = Rm1; ym1 = Rm1; yv1 = Rm1; rhoW1 = Rm1;
% Ratio of mean infections to that from prevalence
ymrho1 = Rm1; yvrho1 = Rm1; Rmrho1 = Rm1; Rvrho1 = Rm1; 

% Preallocate function to prevent overhead
weightFn = @(whetj, rhoi) weightDirichlet(n, nsize, m, k, whetj, nSamps,...
    rhoi, dom, Rlim);

% Slice of outputs at a given weight choice
parfor i = 1:lenrho
    % Specific k and weight set to be evaluated and get statistics
    [ym1{i}, yv1{i}, Rm1{i}, Rv1{i}, ymh1{i}, yvh1{i}, Rmh1{i}, Rvh1{i},...
        rhoW1{i}] = weightFn(whet0, rho(i));

    % Ratios to fixed prevalence case
    Rmrho1{i} = Rmh1{i}./Rm1{i}; Rvrho1{i} = Rvh1{i}./Rv1{i};
    ymrho1{i} = ymh1{i}./ym1{i}; yvrho1{i} = yvh1{i}./yv1{i};
    disp(['Completed ' num2str(i) ' of ' num2str(lenrho)]);
end

% Get quantiles of reproduction numbers and infection statistics
Rmh1q = cellfun(@(x) quantile(x, quant), Rmh1, 'UniformOutput', false);
ymh1q = cellfun(@(x) quantile(x, quant), ymh1, 'UniformOutput', false);
Rvh1q = cellfun(@(x) quantile(x, quant), Rvh1, 'UniformOutput', false);
yvh1q = cellfun(@(x) quantile(x, quant), yvh1, 'UniformOutput', false);

% Quantiles of sampled weights and sum of mean infections
rhoW1q = cellfun(@(x) quantile(x, quant), rhoW1, 'UniformOutput', false);
ysum1 = cellfun(@(x) sum(x,2), ymh1, 'UniformOutput', false);
% Ratios over fixed prevalence quantiles
Rmrho1q = cellfun(@(x) quantile(x, quant), Rmrho1, 'UniformOutput', false);
ymrho1q = cellfun(@(x) quantile(x, quant), ymrho1, 'UniformOutput', false);
yvrho1q = cellfun(@(x) quantile(x, quant), yvrho1, 'UniformOutput', false);

%% %% Infections and R at various Dirichlet shapes for given prevalence

% Dirichlet weighted import rate statistics
Rmh2 = cell(1, lenw); ymh2 = Rmh2; Rvh2 = Rmh2; yvh1 = Rmh2;
% Fixed prevalence import rate statistics and weights
Rm2 = Rmh2; Rv2 = Rm2; ym2 = Rm2; yv2 = Rm2; rhoW2 = Rm2;
% Ratio of mean infections to that from prevalence 
ymrho2 = Rm2; yvrho2 = Rm2; Rmrho2 = Rm2; Rvrho2 = Rm2; 

% Slice of outputs at a given weight choice
parfor i = 1:lenw
    % Specific k and weight set to be evaluated and get statistics
    [ym2{i}, yv2{i}, Rm2{i}, Rv2{i}, ymh2{i}, yvh2{i}, Rmh2{i}, Rvh2{i},...
        rhoW2{i}] = weightFn(whet(i), rho0);

    % Ratios to fixed prevalence case
    Rmrho2{i} = Rmh2{i}./Rm2{i}; Rvrho2{i} = Rvh2{i}./Rv2{i};
    ymrho2{i} = ymh2{i}./ym2{i}; yvrho2{i} = yvh2{i}./yv2{i};
    disp(['Completed ' num2str(i) ' of ' num2str(lenw)]);
end

% Get quantiles of reproduction numbers and infection statistics
Rmh2q = cellfun(@(x) quantile(x, quant), Rmh2, 'UniformOutput', false);
ymh2q = cellfun(@(x) quantile(x, quant), ymh2, 'UniformOutput', false);
Rvh2q = cellfun(@(x) quantile(x, quant), Rvh2, 'UniformOutput', false);
yvh2q = cellfun(@(x) quantile(x, quant), yvh2, 'UniformOutput', false);

% Quantiles of sampled weights and sum of mean infections
rhoW2q = cellfun(@(x) quantile(x, quant), rhoW2, 'UniformOutput', false);
ysum2 = cellfun(@(x) sum(x,2), ymh2, 'UniformOutput', false);
% Ratios over fixed prevalence quantiles
Rmrho2q = cellfun(@(x) quantile(x, quant), Rmrho2, 'UniformOutput', false);
Rvrho2q = cellfun(@(x) quantile(x, quant), Rvrho2, 'UniformOutput', false);
ymrho2q = cellfun(@(x) quantile(x, quant), ymrho2, 'UniformOutput', false);
yvrho2q = cellfun(@(x) quantile(x, quant), yvrho2, 'UniformOutput', false);

%% At small and large n from range show tail probabilities

% Define two extreme event sizes (high H and low L)
nval = [24 48]; domL = 0:nval(1); domH = 0:nval(2);
idL = find(n == nval(1), 1, 'first'); idH = find(n == nval(2), 1, 'first'); 

% Get heterogeneous samples of y at fixed k
[~, ~, ~, yL, RL] = getRxSSEsamp(nval(1), R0, tau, d, k, 10^4);
[~, ~, ~, yH, RH] = getRxSSEsamp(nval(2), R0, tau, d, k, 10^4);

% For given shape show how prevalence influences tails
PxrhoL1 = zeros(lenrho, nval(1)+1); PxrhoH1 = zeros(lenrho, nval(2)+1);
for i = 1:lenrho
    % Import probabilities at different overall prevalence
    PxrhoL1(i, :) = binopdf(domL, nval(1), rho(i));
    PxrhoH1(i, :) = binopdf(domH, nval(2), rho(i));
end

% Survival probabilities indicating superspreading at fixed r
syL1 = cell(1, lenrho); syH1 = syL1; syLx1 = syL1; syHx1 = syL1; yHi1 = syL1;
sRL1 = cell(1, lenrho); sRH1 = sRL1; sRLx1 = sRL1; sRHx1 = sRL1; yLi1 = syL1;
for i = 1:lenrho
    % Thinned by probability of imports
    [yLi1{i}, ~] = impWeight(nval(1), yL, RL, PxrhoL1(i, :), 10^4);
    [yHi1{i}, ~] = impWeight(nval(2), yH, RH, PxrhoH1(i, :), 10^4);
    % Empirical CDFs
    [syL1{i}, syLx1{i}] = ecdf(yLi1{i}); [syH1{i}, syHx1{i}] = ecdf(yHi1{i});
end

% For given overall rho but various shapes changing event based rho
rhoWM = cellfun(@mean, rhoW2, 'UniformOutput', false);
rhoL = zeros(1, lenw); rhoH = rhoL;
PxrhoL2 = zeros(lenw, nval(1)+1); PxrhoH2 = zeros(lenw, nval(2)+1);
for i = 1:lenw
    % Actual mean prevalence from risk awareness
    rhoL(i) = rhoWM{i}(idL); rhoH(i) = rhoWM{i}(idH);
    % Resulting import probabilities
    PxrhoL2(i, :) = binopdf(domL, nval(1), rhoL(i));
    PxrhoH2(i, :) = binopdf(domH, nval(2), rhoH(i));
end

% Survival probabilities indicating superspreading at fixed rho
syL2 = cell(1, lenw); syH2 = syL2; syLx2 = syL2; syHx2 = syL2; yHi2 = syL2;
sRL2 = cell(1, lenw); sRH2 = sRL2; sRLx2 = sRL2; sRHx2 = sRL2; yLi2 = syL2;
for i = 1:lenrho
    % Thinned by probability of imports
    [yLi2{i}, ~] = impWeight(nval(1), yL, RL, PxrhoL2(i, :), 10^4);
    [yHi2{i}, ~] = impWeight(nval(2), yH, RH, PxrhoH2(i, :), 10^4);
    % Empirical CDFs
    [syL2{i}, syLx2{i}] = ecdf(yLi2{i}); [syH2{i}, syHx2{i}] = ecdf(yHi2{i});
end

%% Figure S3: Strength of risk awareness and prevalence statistics

figure('Renderer', 'painters', 'Position', [10 10 800 800]);
% Panels on how rho at fixed r (whet) controls y and R
subplot(4, 2, 1); hold on;
for i = 2:lenrho-1
    plot(n, Rmh1q{i}, 'LineWidth', 2, 'Color', grey1);
end
plot(n, Rmh1q{1}, 'LineWidth', 2, 'Color', 'b');
plot(n, Rmh1q{end}, 'LineWidth', 2, 'Color', 'r');
plot(n, Rlim*ones(size(n)), 'k--', 'LineWidth', 2);
hold off; box off; grid off;
ylabel('E$[R]$', 'FontSize', fnt); 
xlim([min(n) max(n)]); %ylim([0 Rlim+0.01]);

subplot(4, 2, 3); hold on;
for i = 2:lenrho-1
    plot(n, ymh1q{i}, 'LineWidth', 2, 'Color', grey1);
end
plot(n, ymh1q{1}, 'LineWidth', 2, 'Color', 'b');
plot(n, ymh1q{end}, 'LineWidth', 2, 'Color', 'r');
hold off; box off; grid off;
ylabel('E$[y]$', 'FontSize', fnt); xlim([min(n) max(n)]);
legend('', '', '', '', '', '', '', '', ['$\rho$ = ' num2str(rho(1))],...
    ['$\rho$ = ' num2str(rho(end))], fontsize=fnt-1);

subplot(4, 2, 5); hold on;
for i = 2:lenrho-1
    plot(n, Rvh1q{i}, 'LineWidth', 2, 'Color', grey1);
end
plot(n, Rvh1q{1}, 'LineWidth', 2, 'Color', 'b');
plot(n, Rvh1q{end}, 'LineWidth', 2, 'Color', 'r');
hold off; box off; grid off;
ylabel('V$[R]$', 'FontSize', fnt); xlim([min(n) max(n)]);

subplot(4, 2, 7); hold on;
for i = 2:lenrho-1
    plot(n, yvh1q{i}, 'LineWidth', 2, 'Color', grey1);
end
plot(n, yvh1q{1}, 'LineWidth', 2, 'Color', 'b');
plot(n, yvh1q{end}, 'LineWidth', 2, 'Color', 'r');
hold off; box off; grid off;
ylabel('V$[y]$', 'FontSize', fnt); xlim([min(n) max(n)]);
xlabel(['$n | r = $' num2str(whet0)], 'FontSize', fnt);

% Panels on r (whet) at fixed rho controls y and R
subplot(4, 2, 2); hold on;
for i = 2:lenw-1
    plot(n, Rmh2q{i}, 'LineWidth', 2, 'Color', grey1);
end
plot(n, Rmh2q{1}, 'LineWidth', 2, 'Color', 'b');
plot(n, Rmh2q{end}, 'LineWidth', 2, 'Color', 'r');
plot(n, Rlim*ones(size(n)), 'k--', 'LineWidth', 2);
hold off; box off; grid off;
xlim([min(n) max(n)]); ylim([0 Rlim+0.01]);

subplot(4, 2, 4); hold on;
for i = 2:lenrho-1
    plot(n, ymh2q{i}, 'LineWidth', 2, 'Color', grey1);
end
plot(n, ymh2q{1}, 'LineWidth', 2, 'Color', 'b');
plot(n, ymh2q{end}, 'LineWidth', 2, 'Color', 'r');
hold off; box off; grid off;
xlim([min(n) max(n)]); 
legend('', '', '', '', '', '', '', '', ['$r$ = ' num2str(whet(1))],...
    ['$r$ = ' num2str(whet(end))], fontsize=fnt-1);

subplot(4, 2, 6); hold on;
for i = 2:lenrho-1
    plot(n, Rvh2q{i}, 'LineWidth', 2, 'Color', grey1);
end
plot(n, Rvh2q{1}, 'LineWidth', 2, 'Color', 'b');
plot(n, Rvh2q{end}, 'LineWidth', 2, 'Color', 'r');
hold off; box off; grid off;
xlim([min(n) max(n)]);

subplot(4, 2, 8); hold on;
for i = 2:lenrho-1
    plot(n, yvh2q{i}, 'LineWidth', 2, 'Color', grey1);
end
plot(n, yvh2q{1}, 'LineWidth', 2, 'Color', 'b');
plot(n, yvh2q{end}, 'LineWidth', 2, 'Color', 'r');
hold off; box off; grid off;
xlim([min(n) max(n)]); 
xlabel(['$n | \rho = $' num2str(rho0)], 'FontSize', fnt);

%% Figure 4: Risk awareness and prevalence and superspreading

figure('Renderer', 'painters', 'Position', [10 10 800 800]);
% Panels on E[y] and V[y] ratios across rho ar given r
subplot(3, 2, 1); hold on;
for i = 2:lenrho-1
    plot(n, ymrho1q{i}, 'LineWidth', 2, 'Color', grey1);
end
plot(n, ymrho1q{1}, 'LineWidth', 2, 'Color', 'b');
plot(n, ymrho1q{end}, 'LineWidth', 2, 'Color', 'r');
plot(n, ones(size(n)), 'k--', 'LineWidth', 2);
h = gca; ylim = h.YLim;
plot([nval; nval], [ylim; ylim]', '--', 'LineWidth', 2, 'Color', grey2);
hold off; box off; grid off;
ylabel('E$[y | \epsilon]$ E$[y | \rho]^{-1}$', 'FontSize', fnt);
xlim([min(n) max(n)]); 
xlabel(['$n | r = $' num2str(whet0)], 'FontSize', fnt);

subplot(3, 2, 2); hold on;
for i = 2:lenrho-1
    plot(n, ymrho2q{i}, 'LineWidth', 2, 'Color', grey1);
end
plot(n, ymrho2q{1}, 'LineWidth', 2, 'Color', 'b');
plot(n, ymrho2q{end}, 'LineWidth', 2, 'Color', 'r');
plot(n, ones(size(n)), 'k--', 'LineWidth', 2);
h = gca; ylim = h.YLim;
plot([nval; nval], [ylim; ylim]', '--', 'LineWidth', 2, 'Color', grey2);
hold off; box off; grid off;
ylabel('E$[y | \epsilon]$ E$[y | \rho]^{-1}$', 'FontSize', fnt);
xlim([min(n) max(n)]); 
xlabel(['$n | \rho = $' num2str(rho0)], 'FontSize', fnt);

% Panels on survival probabilities for low and high n
subplot(3, 2, 3); hold on;
for i = 2:lenrho-1
    stairs(syLx1{i}, log(1-syL1{i}), '-', 'LineWidth', 2, 'Color', grey1);
end
stairs(syLx1{1}, log(1-syL1{1}), '-', 'LineWidth', 2, 'Color', 'b');
stairs(syLx1{end}, log(1-syL1{end}), '-', 'LineWidth', 2, 'Color', 'r');
hold off; box off; grid off;
ylabel('$\log$ P$(y \geq c)$', 'FontSize', fnt);
xlabel(['$c | n = $ ' num2str(nval(1))], 'FontSize', fnt);
xlim([0 8]); h = gca; h.YLim = [-8 0];
legend('', '', '', '', '', '', '', '', ['$\rho$ = ' num2str(rho(1))],...
    ['$\rho$ = ' num2str(rho(end))], fontsize=fnt-1);


subplot(3, 2, 5); hold on;
for i = 2:lenrho-1
    stairs(syHx1{i}, log(1-syH1{i}), '-', 'LineWidth', 2, 'Color', grey1);
end
stairs(syHx1{1}, log(1-syH1{1}), '-', 'LineWidth', 2, 'Color', 'b');
stairs(syHx1{end}, log(1-syH1{end}), '-', 'LineWidth', 2, 'Color', 'r');
hold off; box off; grid off;
ylabel('$\log$ P$(y \geq c)$', 'FontSize', fnt);
xlabel(['$c | n = $ ' num2str(nval(2))], 'FontSize', fnt);
xlim([0 8]); h = gca; h.YLim = [-8 0];

subplot(3, 2, 4); hold on;
for i = 2:lenw-1
    stairs(syLx2{i}, log(1-syL2{i}), '-', 'LineWidth', 2, 'Color', grey1);
end
stairs(syLx2{1}, log(1-syL2{1}), '-', 'LineWidth', 2, 'Color', 'b');
stairs(syLx2{end}, log(1-syL2{end}), '-', 'LineWidth', 2, 'Color', 'r');
hold off; box off; grid off;
ylabel('$\log$ P$(y \geq c)$', 'FontSize', fnt);
xlabel(['$c | n = $ ' num2str(nval(1))], 'FontSize', fnt);
xlim([0 8]); h = gca; h.YLim = [-8 0];
legend('', '', '', '', '', '', '', '', ['$r$ = ' num2str(whet(1))],...
    ['$r$ = ' num2str(whet(end))], fontsize=fnt-1);

subplot(3, 2, 6); hold on;
for i = 2:lenw-1
    stairs(syHx2{i}, log(1-syH2{i}), '-', 'LineWidth', 2, 'Color', grey1);
end
stairs(syHx2{1}, log(1-syH2{1}), '-', 'LineWidth', 2, 'Color', 'b');
stairs(syHx2{end}, log(1-syH2{end}), '-', 'LineWidth', 2, 'Color', 'r');
hold off; box off; grid off;
ylabel('$\log$ P$(y \geq c)$', 'FontSize', fnt);
xlabel(['$c | n = $ ' num2str(nval(2))], 'FontSize', fnt);
xlim([0 8]); h = gca; h.YLim = [-8 0];

% Running time of script
tsim = toc/60; disp(['Run time = ' num2str(tsim)]);