%% 
%
% Lab 1 Part c
%
% In this final part of the lab we will explore Bayesian parameter estimation. 
% The ML and MAP solutions only give a point estimate of the parameters where as
% the Bayesian approach gives a full distribution over the parameter space. 
% Additionally, by using a conjugate prior we are guaranteed a closed form expression for
% this posterior distribution over the parameters. 
%
%%
 
clear all
close all
 
mu = 1;
sigma = 1; %standard deviation = sqrt(variance)
N = 5;
X = mu + sigma.*randn(N,1);
 
% parameters for prior - normal inverse gamma
alpha = 1;
beta = 1;
gamma = 1;
delta = 0;
 
% parameters for posterior
% TODO 1 - Define the parameters for the Bayesian posterior
alphaP = alpha + N/2;
gammaP = gamma + N;
deltaP = (gamma*delta + sum(X)) / (gamma + N);
betaP = sum(X.^2)/2 + beta + gamma*delta^2 / 2 - (gamma*delta + sum(X))^2/(2*(gamma+N));
 
 
%% ML + MAP Parameter Estimation of univariate normal
% TODO - fill these equations from the previous section
muML = sum(X) / N;
sigmaML = sqrt(sum((X - muML).^2) / N);
muMAP = (N * muML + gamma*delta)/(N+gamma);
sigmaMAP = sqrt((sum((X-muMAP).^2) + 2*beta + gamma((delta-mu)^2)) / (N + 3 + 2*alpha));
 
 
%% Likelihood function
sigmaRange = [ 0 : 0.01 : 2 ];
muRange = [ -2 : 0.01 : 2 ];
conjPosterior = zeros(length(sigmaRange), length(muRange));
prior = zeros(length(sigmaRange), length(muRange));
lfun = zeros(length(sigmaRange), length(muRange));
for s=1:length(sigmaRange)
    for m=1:length(muRange)
        lfun(s,m) = normal( X, sigmaRange(s), muRange(m) );                
        prior(s,m) = normalInvGamma( alpha, beta, delta, gamma, sigmaRange(s), muRange(m) ); 
        
        % TODO 2 - Compute the posterior given the new closed form expression
        conjPosterior(s,m) = normalInvGamma( alphaP, betaP, deltaP, gammaP, sigmaRange(s), muRange(m) );
    end
end
% TODO - Estimate the posterior
posterior = zeros(size(prior));
posterior = lfun.*prior;
% TODO 3 - Empirically compare the closed form posterior to the product of
% likelihood and prior from Part b
% TODO 4 - Show that peak of this distribution again corresponds to the MAP
% solution
 
figure
subplot(1,2,1)
imagesc(muRange,sigmaRange, conjPosterior);colormap hot
hold on
plot(muML, sigmaML, 'bx') % ML parameters
plot(muMAP, sigmaMAP, 'mx') % MAP parameters
ylabel('sigma');xlabel('mu')
legend('ML', 'MAP')
title('closed form estimation of posterior')
 
subplot(1,2,2)
imagesc(muRange,sigmaRange, posterior);colormap hot
hold on
plot(muML, sigmaML, 'bx') % ML parameters
plot(muMAP, sigmaMAP, 'mx') % MAP parameters
ylabel('sigma');xlabel('mu')
legend('ML', 'MAP')
title('posterior = likelihood*prior')
 
% BONUS TODO 5 - Write code to estimate the probability that a new data point
% belongs to the same model. Compare the ML, MAP and fully Bayesian
% methods to do this. What are the advantages of using the Bayesian
% approach?
 


