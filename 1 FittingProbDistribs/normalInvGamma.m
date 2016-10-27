function res = normalInvGamma( alpha, beta, delta, gammaVar, sigma, mu )
% Compute the probability that a given normal was generated using this
% normal inverse gamma funcion.

% TODO fill out this function
res = sqrt(gammaVar) / (sigma * sqrt(2*pi)) * beta^alpha / gamma(alpha) * (1 / sigma^2)^(alpha + 1) * exp (- (2*beta + (gamma(delta-mu))^2)/(2*sigma^2));
        
