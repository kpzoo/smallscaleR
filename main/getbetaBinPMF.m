% Beta-binomial prob mass function
function Pdom = getbetaBinPMF(n, a, b)

% Assumptions and notes
% - input trials (n) and beta distrib hyperparams (a, b)
% - the domain is implicit as 0:n

% BetaBin prob mass across domain
Pdom = zeros(1, n+1); B0 = beta(a, b);

% Evaluation using beta functions
for i = 1:n+1
    % For x successes
    x = i-1; nx = nchoosek(n, x);
    Pdom(i) = nx*beta(x+a, n-x+b)/B0;
end

% Check normalisation
if(abs(sum(Pdom) - 1) > 10^-9)
    assignin('base', 'Pdom', Pdom);
    error('Issues in pmf normalisation');
end

