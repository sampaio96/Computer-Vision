function plotGaussian( mu, sigma, col )
% plot Gaussian
normcurv = @(x) exp( -(x-mu).^2./(2*sigma^2) )/(sigma*sqrt(2*pi));
xPlot = [mu - 4*sigma : 0.01 : mu + 4*sigma];
yPlot = normcurv(xPlot);
plot(xPlot, yPlot, col, 'LineWidth', 1.5)
