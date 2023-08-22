% Compute heterogenous statistics for Dirichlet weights
function [ym, yv, Rm, Rv, ymh, yvh, Rmh, Rvh, rhoW] = weightDirichlet(n, nsize,...
    m, k, whet, nSamps, rhoM, dom, Rlim)

% Assumptions and notes
% - x introduction, y infections, n-x susceptibles
% - apply a set of Dirichlet sized weights alphaW
% - considers one k value and one Dirichlet weight setting
% - total of m events and whet is shape of Dirichlet

% Fixed importation rate statistics
Rm = zeros(1, m); ym = Rm; Rv = Rm; yv = Rm;
% Size biased importation rate statistics
Rmh = zeros(nSamps, m); ymh = Rmh; Rvh = Rmh; yvh = Rmh;

% Weighting of Dirichlet distribution
wn = whet*ones(1, m); pnrho = rhoM*drchrnd(wn, nSamps);
% Sort so larger n has larger weights and normalise by size
rhoW = sort(pnrho, 2, 'ascend')./nsize;

% Import prob function for a given mean and with weights
Pxrhofn = cell(1, m); Pxfix = Pxrhofn; Pxhet = cell(nSamps, m);
for i = 1:m
    %Pxrhofn{i} = @(rho) binopdf(dom{i}, n(i), rho);
    Pxrhofn{i} = @(rho) getBinPMF(n(i), rho);

    % Fixed and constant alpha (prevalence)
    Pxfix{i} = Pxrhofn{i}(rhoM);
    % Get P(x) under the heterogeneous weightings on alphaM
    for j = 1:nSamps
        % Fixed and constant rho (prevalence)
        Pxhet{j, i} = Pxrhofn{i}(rhoW(j, i));
    end
end

% For all n generate P(y|x,n) for all possible x
for i = 1:m
    % Variables for a given n and k
    RM = zeros(1, n(i)+1); RV = RM; yM = RM; yV = RM;

    % For every possible import value compute y and R stats
    for ii = 2:n(i)+1
        % Number of imports x from 1 to n
        x = ii - 1;

        % Elements inside summation over x for means
        RM(ii) = (n(i)/x - 1)*(1 - (1 + Rlim/(k*n(i)))^(-k*x));
        yM(ii) = (n(i) - x)*(1 - (1 + Rlim/(k*n(i)))^(-k*x));

        % Elements inside summation over x for variances
        RV(ii) = (1 + 2*Rlim/(k*n(i)))^(-k*x);
        RV(ii) = RV(ii) - (1 + Rlim/(k*n(i)))^(-2*k*x);
        yV(ii) = ((n(i) - x)^2)*RV(ii);
        RV(ii) = ((n(i)/x - 1)^2)*RV(ii);
    end

    % Statistics weighted properly by P(x) fixed
    Rm(i) = RM*Pxfix{i}'; ym(i) = yM*Pxfix{i}';
    Rv(i) = RV*(Pxfix{i}.^2)'; yv(i) = yV*(Pxfix{i}.^2)';

    % Weight by size aware x distributions
    for j = 1:nSamps
        Rmh(j,i) = RM*Pxhet{j, i}';
        ymh(j,i) = yM*Pxhet{j, i}';
        Rvh(j,i) = RV*(Pxhet{j, i}.^2)';
        yvh(j,i) = yV*(Pxhet{j, i}.^2)';
    end

end