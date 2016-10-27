function r=practicalRegress2

%The goal of this part is to investigate the Bayesian approach to
%regression

%Work your way through these Matlab examples filling in code where it says
%"TO DO"


%close all figures
close all;

%set seeds for random number generators 
%so we get the same random numbers each time
rand('seed',5);
randn('seed',5);


%define actual parameters
offsetActual = 1.5;  % this is phi_0 in the notes
slopeActual = -0.8;  % this is phi_1 in the notes
noiseActual = 0.01;  % this is sigma^2 in the notes


%generate some data 
nData = 5; % I = 5
x = rand(1,nData);
w = offsetActual + slopeActual * x + sqrt(noiseActual) * randn(1,nData);
%Transpose w so that it is a column vector like in the notes
w = w'; % w : (Ix1) = (5x1)
%Add a row of ones to the start of x
x = [ones(1,size(x,2));x]; % x : (DxI) = (2x5)
nDim = size(x,1); % nDim : D = 2


%display this data
figure; set(gcf,'Color',[1 1 1]);
plot(x(2,:),w','r.');
xlabel('Data, x'); ylabel('World, w');
set(gca,'Box','Off'); 
xlim([0 1]);ylim([0 2]);

%define a normal prior on the slope and offset parameters
sigmaPSq = 1;
priorPhiMean = zeros(nDim,1);
priorPhiCov =  eye(nDim)*sigmaPSq;

%fit the noise parameter to the model 
%I've done this for you (it's at the bottom), but you should check out how
%it works
sigmaSq= estNoiseParameter(x,w,sigmaPSq);

%now the problem is to fit the data - we find a normal  probability distribution
%over the parameter vector phi
%TODO  - fill in this routine (template at the bottom of this file)
[posteriorPhiMean posteriorPhiCov] = fitLinearRegressionBayes(x,w,sigmaSq,sigmaPSq)


%********************************
%display the prior and posterior over the parameters
%create arrays of x and w positions
xFigRange = -2.0:0.01:2.0;
wFigRange = -2.0:0.01:2.0;
xFig = repmat(xFigRange,length(wFigRange),1);
wFig = repmat(wFigRange',1,length(xFigRange));
%evaluate the Gaussians before and after
xwDrawData = [xFig(:) wFig(:)]';
prior =  normalProb(xwDrawData,priorPhiMean,priorPhiCov);
prior = reshape(prior,size(xFig));
posterior =  normalProb(xwDrawData,posteriorPhiMean,posteriorPhiCov);
posterior = reshape(posterior,size(xFig));
figure; set(gcf,'Color',[1 1 1]);
subplot(1,2,1); imagesc(prior);
xlabel('Offset'); ylabel('Slope'); title('Prior');axis image;
set(gca,'XTick',[1 length(xFigRange)]); set(gca,'XTickLabel',{'-2.0', '2.0'});
set(gca,'YTick',[1 length(wFigRange)]); set(gca,'YTickLabel',{'-2.0', '2.0'});
subplot(1,2,2); imagesc(posterior);
xlabel('Offset'); ylabel('Slope'); title('Posterior'); axis image;
set(gca,'XTick',[1 length(xFigRange)]); set(gca,'XTickLabel',{'-2.0', '2.0'});
set(gca,'YTick',[1 length(wFigRange)]); set(gca,'YTickLabel',{'-2.0', '2.0'});
%********************************


%Let's draw a figure to show the posterior probability 
%of the world as a function of the data
%For every xStar, the model predicts a normal distribution over w
xStar = 0:0.005:1;
%Add ones to the start of xStar
xStar = [ones(1,size(xStar,2));xStar];
[nDim nData] = size(xStar); % 2x201

%TO DO  - Compute A matrix - (using either equation 8.11 or 8.13 from
%notes)
%Replace this:
A = 1/sigmaSq * x * x.' + 1/sigmaPSq * eye(nDim); % A : DxD = 2x2

if ~(size(A, 1) == nDim && size(A, 2) == nDim)
    error('A should be a square matrix, size nDim x nDim');
end

% TO DO Compute the mean and variance of the prediction for each element of
% x2.  Replace this:

for (cX = 1:size(xStar,2)) % xStar(cX) : 2x1
    predMean(cX) = 1/sigmaSq * xStar(:,cX).' * inv(A) * x * w; % 1x1 -> predMean : 1x201
    predVar(cX) = xStar(:,cX).' * inv(A) * xStar(:,cX) + sigmaSq; % 1x1 -> predVar : 1x201
end;

% checking predMean and predVar are of the correct dimensions
if ~(isvector(predMean) && size(predMean, 2)==length(xStar))
    error('predMean should be a vector of same length as xStar');
elseif ~(isvector(predVar) && size(predVar, 2)==length(xStar))
    error('predVar should be a vector of same length as xStar');
end

%compute the pixel colours for the image - one Gaussian in each column
wFig = 0:0.005:2;
postFun = zeros(length(wFig),length(xStar));
for (cX= 1:length(xStar))
    postFun(:,cX)= (1./sqrt(2*pi*predVar(cX))).*exp(-0.5 *((wFig-predMean(cX)).^2)./predVar(cX));
end;

%draw the figure
figure; set(gcf,'Color',[1 1 1]);
imagesc(postFun); colormap(hot);
set(gca,'YDir','Normal');
hold on;
%draw the points on top - they need to be rescaled so that they are in the 
%units of the pixels used to draw the picture
plot((x(2,:)-xStar(2,1))/(xStar(2,2)-xStar(2,1)), (w'-wFig(1))/(wFig(2)-wFig(1)),'bo');
%draw the units on the graph
set(gca,'XTick',[1 length(xStar)]);
set(gca,'XTickLabel',{'0','1'});
set(gca,'YTick',[1 length(wFig)]);
set(gca,'YTickLabel',{'0','1'});
xlabel('x');
ylabel('w');
title('Pr(w|x)');

%TO DO (AT HOME IF YOU ARE KEEN): 
%Convert this to non-linear regression by passing the data
%through a non-linear transformation.  When you do this you'll have to 
%comment out the section between the ***** ***** as the prior and posterior
%are no longer 2 dimensional.
%You'll also have to load in the non-linear data and make sure that the
%final plot is drawn with sensible axes.



%The goal of this routine is to take data x and w and the prior variance
%scale sigmaPSq and return the mean and covariance of the posterior
%distribution of the vector phi = [offset;slope];
function [posteriorPhiMean posteriorPhiCov] = fitLinearRegressionBayes(x,w,sigmaSq,sigmaPSq);

%retrieve dimensionality and number of data points 
[nDim nData] = size(x);


%TO DO compute A matrix (using either equation 8.11 or 8.13 from notes)
%REPLACE THIS
A = 1/sigmaSq * x * x.' + 1/sigmaPSq * eye(nDim);

%TODO compute posterior mean and variance of phi vector (contains offset and slope)
%REPLACE THIS:
posteriorPhiMean = 1/sigmaSq * inv(A) * x * w;
posteriorPhiCov = inv(A);

% checking sizes of output values are correct
if ~(size(posteriorPhiCov, 1) == nDim && size(posteriorPhiCov, 2) == nDim)
    error('posteriorPhiCov should be a nDim x nDim matrix');
elseif ~(isvector(posteriorPhiMean) && length(posteriorPhiMean) == nDim)
    error('posteriorPhiMean should be a nDim x 1 vector');
end

%The goal of this routine is to take data x and w and the prior variance
%scale sigmaPSq and return the estimated variance of the linear regression
%model.  It does this by optimizing the marginal likelihood
function sigmaSq = estNoiseParameter(x,w,sigmaPSq)

%initial estimate for variance
sigmaSqInit = std(w);
%take logarithm (so optimization routine doesn`t have to constrain to be positive)
logSigmaSqInit = log(sigmaSqInit);
%minimizes the negative log marginal likelihood
logSigmaSq = fminsearch(@(logSigmaSq) negLogMarginalLike(logSigmaSq, x,w,sigmaPSq),logSigmaSqInit);
%exponentiate to return to positive value
sigmaSq = exp(logSigmaSq);


%This function returns the negative log marginal likelihood
function L = negLogMarginalLike(logSigmaSq,x,w,sigmaPSq)

sigmaSq =exp(logSigmaSq);
L = -sum(log(normalProb(w,zeros(size(w)),x'*x*sigmaPSq+eye(size(x,2))*sigmaSq)));


%evaluates data against the multivariate normal distribution
function l = normalProb(x,normalMean,normCov);

[nDim nData] = size(x);

x = x-repmat(normalMean,1,nData);
l = sum((inv(normCov)*x).*x);
l = (1/sqrt(det(normCov)*((2*pi)^nDim)))*exp(-0.5*l);
