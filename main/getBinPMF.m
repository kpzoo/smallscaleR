% Binomial prob mass function
function Pdom = getBinPMF(n, p)

% Assumptions and notes
% - input trials (n) and prob of success p
% - the domain is implicit as 0:n

% BetaBin prob mass across domain
Pdom = zeros(1, n+1); 

% Evaluation using beta functions
for i = 1:n+1
    % For x successes
    x = i-1; nx = nchoosek(n, x);
    Pdom(i) = nx*(p^x)*(1-p)^(n-x);
end

% Check normalisation
if(abs(sum(Pdom) - 1) > 10^-9)
    assignin('base', 'Pdom', Pdom);
    error('Issues in pmf normalisation');
end
