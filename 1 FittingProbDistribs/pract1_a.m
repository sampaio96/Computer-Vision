%% 
%
% Lab 1 Part a
%
% We will use maximum likelihood (ML) to fit both parameters (mu and sigma) of a
% univariate normal distribution to the data.
%
%%
 
clear all
close all
 
% generate some data points from a normal distribution
mu = 1; % mean
sigma = 2; %standard deviation = sqrt(variance)
N = 20; % number of datapoints
X = mu + sigma.*randn(N,1);
 
%% ML Parameter Estimation of normal
% TODO 1 - On paper derive the ML parameter estimates for the normal
% distribution and then insert the equations here.
% Note: For today, it's safe to assume that the data points were generated
% independently of each other.
muML = 0;    % So put your answer here in place of zero.
sigmaML = 0; % So put your answer here in place of zero.
 
% Plotting univariate normal
% normalizing histograms - assuming bins are evenly sized
[histCnt, histPos] = hist(X, 20);
binWidth = histPos(2)-histPos(1);
bar(histPos, histCnt./(sum(histCnt)*binWidth), 1);
hold on;plotNormal(mu, sigma, 'g')
hold on;plotNormal(muML, sigmaML, 'r')
legend('data', 'input', 'ML estimate')
title('ML parameter fited distribution')
hold off

% TODO remove this break when you have completed the first section
break
 
%% Likelihood function
% 
offSet = 1;
sigmaRange = [sigma-offSet : 0.01 : sigma+offSet];
muRange = [mu-offSet : 0.01 : mu+offSet];
 
% Note: in the code below, there may be nothing missing, but functions like 
% normal() and logNormal() are being called, yet are potentially blank. To 
% accomplish these TODO's, fill them in!
% TODO 2 - Compute the likelihood function for the normal distribution
% TODO 3 - Empirically verify that the maximum is at the same position as
% your ML solution
% TODO 4 - Compute the log-likelihood function and verify that it has the
% same maximum
lfun = zeros(length(sigmaRange), length(muRange));
llfun = zeros(length(sigmaRange), length(muRange));
for s=1:length(sigmaRange)
    for m=1:length(muRange)
        lfun(s,m) = normal( X, sigmaRange(s), muRange(m) );
        llfun(s,m) = logNormal( X, sigmaRange(s), muRange(m) );
    end
end
 
% Plotting log-likelihood function
figure
subplot(1,2,1)
imagesc(muRange,sigmaRange,llfun);colormap hot
[val ind] = max(llfun(:)); % find max
[maxSigma maxMu] = ind2sub(size(llfun), ind);
hold on
plot(muRange(maxMu), sigmaRange(maxSigma), 'bo')% peak
plot(muML, sigmaML, 'rx') % ML parameters
plot(mu, sigma, 'go') % input
title('Log-likelihood function');
ylabel('sigma');xlabel('mu')
legend('peak', 'ml estimate', 'input');
 
% Plotting likelihood function
subplot(1,2,2)
imagesc(muRange,sigmaRange, lfun);colormap hot
[val ind] = max(lfun(:)); % find max
[maxSigma maxMu] = ind2sub(size(lfun), ind);
hold on
plot(muRange(maxMu), sigmaRange(maxSigma), 'bo')% peak
plot(muML, sigmaML, 'rx') % ML parameters
plot(mu, sigma, 'go') % input
title('Likelihood function');
ylabel('sigma');xlabel('mu')
legend('peak', 'ml estimate', 'input');
 
 
% TODO 5 - Observe what happens to the likelihood function as you increase the number of datapoints
 
 
