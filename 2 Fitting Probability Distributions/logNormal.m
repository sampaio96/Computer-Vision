function res = logNormal( X, sigma, mu )
% Computes the likelihood that the data X have been generated using the
% given parameters of a one-dimensional normal, but in log space, so the 
% equation should be derived and implemented differently here than in
% normal().
N = length(X);
% TODO fill out this function
res = -0.5 * N * log(2*pi) - 0.5 * N * log (sigma^2) - 0.5 * sum (((X-mu).^2)./sigma^2);