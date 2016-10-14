%% 
%
% Lab 1 Part b
%
% In this part of the lab we will be using maximum a posteriori (MAP). This
% allows us to introduce prior information we may have about the
% parameters. We shall use the Normal Inverse Gamma as our conjugate prior.
%
%
%  PLEASE limit your use of built-in functions, so use built-in gamma that 
% we mentioned in class if you like, but not more complicated or 
% conveneince functions like normpdf, normlike etc.
%%
 
clear all
close all
 
% generate some data points from a normal distribution
mu = 1; % mean
sigma = 1; % standard deviation = sqrt(variance)
N = 5; % number of datapoints
X = mu + sigma.*randn(N,1);
 
% parameters for prior - normal inverse gamma
alpha = 1;
beta = 1;
gamma = 1;
delta = 0;
 
 
%% MAP Parameter Estimation of univariate normal
% TODO 1 - On paper derive the MAP parameter estimates for the normal
% distribution using the normal inverse gamma prior and then insert the
% equations below
muML = 0;
sigmaML = 0;
muMAP = 0;
sigmaMAP = 0;
 
% TODO remove this break when you have completeed the first section
break
 
%% Likelihood function
sigmaRange = [0 : 0.01 : 2];
muRange = [ -2: 0.01 : 2 ];
% Note: again the TODO's are calling for you to dig into the functions
% called within the for-loops.
% TODO 2 - Compute the likelihood function for the normal distribution
% using the same code from Part a
% TODO 3 - Compute the likelihood function for the prior
lfun = zeros(length(sigmaRange), length(muRange));
prior = zeros(length(sigmaRange), length(muRange));
for s=1:length(sigmaRange)
    for m=1:length(muRange)        
        lfun(s,m) = normal( X, sigmaRange(s), muRange(m) );                
        prior(s,m) = normalInvGamma( alpha, beta, delta, gamma, sigmaRange(s), muRange(m) );                   
    end
end
% TODO 4 - Compute the posterior
% TODO 5 - Empirically verify that the maximum of the posterior is at the same position as
% your MAP solution
posterior = zeros(size(prior));
 
% Plotting the likelihood, prior and posterior
figure
subplot(2, 2, 1)
imagesc(muRange,sigmaRange, lfun);colormap hot
hold on
plot(muML, sigmaML, 'bx') % ML parameters
ylabel('sigma');xlabel('mu')
title('likelihood')
 
subplot(2, 2, 2)
imagesc(muRange,sigmaRange, prior);colormap hot
ylabel('sigma');xlabel('mu')
title('prior')
 
subplot(2, 2, 3)
imagesc(muRange,sigmaRange, posterior);colormap hot
hold on
plot(muML, sigmaML, 'bx') % ML parameters
plot(muMAP, sigmaMAP, 'mx') % MAP parameters
ylabel('sigma');xlabel('mu')
title('posterior')
 
subplot(2, 2, 4)
imagesc(zeros(size(posterior)));
hold on
plot(muML, sigmaML, 'bx') % ML parameters
plot(muMAP, sigmaMAP, 'mx') % MAP parameters
legend('ML', 'MAP')
axis off
 
% TODO 6 - Comment on the effects of introducing the prior as the number of
% datapoints is low as compared to high
 
 


